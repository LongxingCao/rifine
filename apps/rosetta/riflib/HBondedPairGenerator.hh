// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:

#ifndef INCLUDED_riflib_HBondedPairGenerator_hh
#define INCLUDED_riflib_HBondedPairGenerator_hh

// headers

	// #include <core/pack/rotamer_set/RotamerSet.hh>
	// #include <core/pack/rotamer_set/RotamerSetFactory.hh>
	// #include <core/pack/dunbrack/RotamerLibrary.hh>
	#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
	// #include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>
	// #include <core/pack/task/TaskFactory.hh>
	// #include <core/graph/Graph.hh>
	// #include <core/pack/packer_neighbors.hh>

	#include <ObjexxFCL/format.hh>
	#include <ObjexxFCL/string.functions.hh>

	#include <core/chemical/AtomType.hh>
	#include <core/chemical/orbitals/OrbitalType.hh>
	#include <core/chemical/ChemicalManager.hh>
	#include <core/chemical/ResidueTypeSet.hh>
	// #include <core/conformation/symmetry/util.hh>
	#include <core/conformation/ResidueFactory.hh>
	// #include <core/conformation/util.hh>
	// #include <core/import_pose/import_pose.hh>
	// #include <core/io/silent/SilentFileData.hh>
	// #include <core/pose/PDBInfo.hh>
	#include <core/pose/Pose.hh>
	// #include <core/pose/annotated_sequence.hh>
	#include <core/pose/motif/reference_frames.hh>
	// #include <core/pose/util.hh>
	// #include <core/pose/symmetry/util.hh>
	// #include <core/import_pose/import_pose	.hh>
	// #include <core/kinematics/MoveMap.hh>
	#include <core/scoring/Energies.hh>
	// #include <core/scoring/EnergyGraph.hh>
	#include <core/scoring/ScoreFunction.hh>
	#include <core/scoring/ScoreFunctionFactory.hh>
	// #include <core/scoring/ScoreTypeManager.hh>
	// #include <core/scoring/dssp/Dssp.hh>
	// #include <core/scoring/etable/Etable.hh>
	#include <core/scoring/hbonds/HBondOptions.hh>
	// #include <core/scoring/hbonds/HBondSet.hh>
	// #include <core/scoring/hbonds/hbonds.hh>
	#include <core/scoring/motif/util.hh>

	#include <core/scoring/methods/EnergyMethodOptions.hh>
	// #include <core/scoring/packing/compute_holes_score.hh>
	// #include <core/scoring/rms_util.hh>
	// #include <core/scoring/sasa.hh>
	// #include <core/scoring/symmetry/SymmetricScoreFunction.hh>
	#include <numeric/conversions.hh>
	#include <numeric/model_quality/rms.hh>
	// #include <numeric/random/random.hh>
	#include <numeric/xyz.functions.hh>
	#include <numeric/xyz.io.hh>
	#include <numeric/xyzVector.hh>
	// #include <protocols/idealize/IdealizeMover.hh>
	// #include <protocols/sicdock/Assay.hh>
	// #include <protocols/sic_dock/SICFast.hh>
	// #include <protocols/sic_dock/util.hh>
	// #include <protocols/sic_dock/read_biounit.hh>
	// #include <protocols/simple_moves/MinMover.hh>
	// #include <protocols/simple_moves/AddConstraintsToCurrentConformationMover.hh>
	// #include <utility/io/izstream.hh>
	// #include <utility/io/ozstream.hh>
	#include <utility/tools/make_vector1.hh>
	// #include <utility/fixedsizearray1.hh>
	// #include <utility/file/file_sys_util.hh>
	// #include <numeric/geometry/hashing/SixDHasher.hh>
	// #include <numeric/HomogeneousTransform.hh>

	#include <riflib/RotamerGenerator.hh>

	// #include <apps/pilot/will/will_util.ihh>

	#include <boost/foreach.hpp>
	#include <boost/iterator/iterator_facade.hpp>


namespace devel {
namespace scheme {

core::Real
get_rot_score(
	core::pose::Pose & pose,
	core::Size ir,
	core::pack::dunbrack::RotamerLibraryScratchSpaceOP scratch
);

void dump_pdb_atom(
	std::ostream & out,
	double x, double y, double z
);


///@brief generate BCC lattice "hyper-ring" to sample rotationr around
/// the hbond axis as well as tilting the axis a up to short_tol
utility::vector1<numeric::xyzMatrix<double> >
get_rotation_samples( double short_tol, double long_tol, double reslang );


class HBondedPairGenerator
    : public boost::iterator_facade< HBondedPairGenerator,
                                     core::pose::Pose,
                                     boost::forward_traversal_tag,
                                     core::pose::Pose const &
             						>
{
public:
	core::pose::Pose pose_,pose0_;
	core::pack::dunbrack::RotamerLibraryScratchSpaceOP scratch_;
	// core::pack::rotamer_set::RotamerSetOP rotset1_,rotset2_;
	utility::vector1< utility::vector1< float > > rotset1_, rotset2_;
	core::Size irot1_,irot2_,idon_,iacc_,iorb_,ihbr_,irotindex_start1_,irotindex_start2_;
	core::Size nrots1_,nrots2_;
	core::Size clash_atom_start1_, clash_atom_start2_;
	utility::vector1<core::Size> donor_atoms_;
	utility::vector1<core::Size> donor_bases_;
	utility::vector1<core::Size> acceptor_atoms_;
	utility::vector1<utility::vector1<core::Size> > acceptor_orbitals_;
	utility::vector1<numeric::xyzMatrix<double> > rot_samples_;
	core::scoring::ScoreFunctionOP score_func_;
	core::id::AtomID align_atom1_, align_atom2_, align_atom3_;
	core::Real score_;
	bool fix_donor_, fix_acceptor_;
	HBondedPairGenerator(){ irot1_=0; irot2_=0; idon_=0; iacc_=0; iorb_=0; ihbr_=0; }
	core::chemical::ResidueTypeOP rtype1op_, rtype2op_;
	virtual ~HBondedPairGenerator(){}
	void init(
		RotamerIndex const & rot_index,
		std::string resn1,
		std::string resn2,
		double tip_tol=20.0,
		double rot_resl=5.0,
		double rot_range=360.0,
		bool fix_donor=false,
		bool fix_acceptor=false,
		core::id::AtomID align_atom1 = core::id::AtomID(1,1),
		core::id::AtomID align_atom2 = core::id::AtomID(2,1),
		core::id::AtomID align_atom3 = core::id::AtomID(3,1),
		std::map<std::string,core::pose::Pose> const & exemplars = std::map<std::string,core::pose::Pose>()
	);
	bool has_more_samples() const { return irot1_!=0; }
	core::Size raw_num_samples() const {
		core::Size nacc_orbs = 0;
		for(core::Size i = 1; i <= acceptor_orbitals_.size(); ++i)
			nacc_orbs += acceptor_orbitals_[i].size();
		return rot_samples_.size()              *
			   nacc_orbs                        *
		       nrots2_                          *
		       donor_atoms_.size()              *
		       nrots1_                          ;
	}
	size_t rot1_of_current() const { return irot1_ + irotindex_start1_ - 1 ; }
	size_t rot2_of_current() const { return irot2_ + irotindex_start2_ - 1 ; }
	core::Real score_of_current() const { return score_; }
private:
    friend class boost::iterator_core_access;
    core::pose::Pose const & dereference() const {
    	return pose_;
    }
    void update_chi(){
    	if( ! fix_donor_ )
			for(core::Size ichi=1; ichi <= pose0_.residue(1).nchi(); ++ichi)
				pose0_.set_chi(ichi,1, rotset1_[irot1_][ichi] );
		if( ! fix_acceptor_ )
			for(core::Size ichi=1; ichi <= pose0_.residue(2).nchi(); ++ichi)
				pose0_.set_chi(ichi,2, rotset2_[irot2_][ichi] );
    }

    virtual bool update_hb_and_score(bool=true) { utility_exit_with_message("Abstract base class"); }

    void increment(){
    	bool clash = true;
    	while(clash){
			// std::cout << "INCR " << irot1_ << " " << idon_ << " " << irot2_ << " " << iacc_ << " " << iorb_ << " " << ihbr_ << std::endl;
    		++ihbr_;
	    	if( ihbr_  > rot_samples_.size()              ){ ihbr_  = 1; ++iorb_ ; }
    		if( iorb_  > acceptor_orbitals_[iacc_].size() ){ iorb_  = 1; ++iacc_ ; }
    		if( iacc_  > acceptor_atoms_.size()           ){ iacc_  = 1; ++irot2_; }
			if( irot2_ > nrots2_                          ){ irot2_ = 1; ++idon_ ; }
			if( idon_  > donor_atoms_.size()              ){ idon_  = 1; ++irot1_;  }
			if( irot1_ > nrots1_                          ){ irot1_ = irot2_ = idon_ = iacc_ = iorb_ = ihbr_ = 0;  return; }
			this->update_chi();
			clash = ! this->update_hb_and_score(ihbr_==1);
		}
    }
    bool equal(HBondedPairGenerator const& o) const {
    	return o.irot1_==irot1_ && o.irot2_==irot2_ && o.iacc_==iacc_ && o.idon_==idon_ && iorb_==o.iorb_ && ihbr_==o.ihbr_;
    }

};

// typedef utility::pointer::shared_ptr<HBondedPairGenerator> HBondedPairGeneratorOP;

class SingleHbondedPairGenerator : public HBondedPairGenerator {

    virtual bool update_hb_and_score(bool realign=true);

};

class BidentateHbondedPairGenerator : public HBondedPairGenerator {

   virtual bool update_hb_and_score(bool realign=true);

};

}
}

#endif

