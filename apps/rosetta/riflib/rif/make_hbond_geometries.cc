// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.



#include <riflib/rif/make_hbond_geometries.hh>


	#include <ObjexxFCL/FArray3D.hh>
	#include <ObjexxFCL/format.hh>

	#include <basic/Tracer.hh>
	#include <basic/options/keys/in.OptionKeys.gen.hh>
	#include <basic/options/keys/mh.OptionKeys.gen.hh>
	#include <basic/options/keys/out.OptionKeys.gen.hh>
	#include <basic/options/keys/packing.OptionKeys.gen.hh>
	#include <basic/options/option_macros.hh>

	#include <boost/assign/std/vector.hpp>
	#include <boost/foreach.hpp>
	#include <boost/lexical_cast.hpp>
	#include <boost/math/special_functions/sign.hpp>
	#include <boost/random/mersenne_twister.hpp>
	#include <boost/random/uniform_real.hpp>

	#include <core/id/AtomID.hh>
	#include <core/import_pose/import_pose.hh>
	#include <core/pose/Pose.hh>
	#include <core/pose/motif/reference_frames.hh>
	#include <core/scoring/ScoreFunctionFactory.hh>
	#include <core/scoring/motif/motif_hash_stuff.hh>
	#include <core/scoring/motif/util.hh>

	#include <devel/init.hh>
	#include <riflib/EtableParams_init.hh>
	#include <riflib/HBondedPairGenerator.hh>
	#include <riflib/RotamerGenerator.hh>
	#include <riflib/rosetta_field.hh>
	#include <riflib/util.hh>

	#include <map>
	#include <parallel/algorithm>

	#include <scheme/actor/Atom.hh>
	#include <scheme/actor/BackboneActor.hh>
	#include <scheme/actor/VoxelActor.hh>
	#include <scheme/chemical/RotamerIndex.hh>
	#include <scheme/kinematics/Director.hh>
	#include <scheme/kinematics/Scene.hh>
	#include <scheme/nest/pmap/OriTransMap.hh>
	#include <scheme/objective/ObjectiveFunction.hh>
	#include <scheme/objective/hash/XformMap.hh>
	#include <scheme/objective/storage/RotamerScores.hh>
	#include <scheme/objective/voxel/FieldCache.hh>
	// #include <scheme/objective/voxel/VoxelArray.hh>
	#include <scheme/rosetta/score/RosettaField.hh>
	#include <scheme/util/StoragePolicy.hh>

	#include <utility/file/file_sys_util.hh>
	#include <utility/io/izstream.hh>
	#include <utility/io/ozstream.hh>
	#include <utility/tools/make_vector1.hh>

using namespace core;

namespace devel {
namespace scheme {
namespace rif {

void make_hbond_geometries(
	RotamerIndex const & rot_index,
	std::string resn1,
	std::string resn2,
	bool fix_donor,
	bool fix_acceptor,
	std::map<std::string,core::pose::Pose> const & exemplars,
	utility::vector1< RelRotPos > & hbonders,
	MakeHbondGeomOpts const & opts
 ){
	using basic::options::option;
	using namespace basic::options::OptionKeys;
	using namespace core::scoring;
	using namespace core::scoring::motif;
		using namespace devel::scheme;
		typedef numeric::xyzVector<core::Real> Vec;
		typedef numeric::xyzVector<float> Vecf;
		typedef numeric::xyzMatrix<core::Real> Mat;
		typedef numeric::xyzTransform<core::Real> Xform;
		using core::id::AtomID;
		using std::cout;
		using std::endl;

		core::Size totcount = 0;

		omp_lock_t init_lock, cout_lock;
		omp_init_lock(&init_lock);
		omp_init_lock(&cout_lock);		

	devel::scheme::HBondedPairGenerator * gen_ptr = new devel::scheme::SingleHbondedPairGenerator;
	devel::scheme::HBondedPairGenerator & gen(*gen_ptr);

	omp_set_lock( &cout_lock );
	cout << "==================================================================================================================" << endl;
	cout << "===================== Thread " << omp_thread_num_1() << " HBondedPairGenerator "
		  << resn1 << " " << resn2 << " ================================================" << endl;
	cout << "==================================================================================================================" << endl;
	omp_unset_lock( &cout_lock );

	int frame_atomno_1 = 1;
	int frame_atomno_2 = 2;
	int frame_atomno_3 = 3;
	bool fixed_is_nonstandard;
	if( fix_donor ){
		core::chemical::ResidueTypeSetCAP rts = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");
		core::chemical::ResidueType const & rtype = rts.lock()->name_map(resn1);
		frame_atomno_1 = rtype.nheavyatoms()-2;
		frame_atomno_2 = rtype.nheavyatoms()-1;
		frame_atomno_3 = rtype.nheavyatoms()-0;
		fixed_is_nonstandard = ( exemplars.find(resn1)!=exemplars.end() );
		cout << "alignment atomnos " << resn1 << " " << frame_atomno_1 << " " << frame_atomno_2 << " " << frame_atomno_3 << endl;
		if( resn1=="SER" ){
			frame_atomno_1 = 11; // SER is special case because so small... use HG instead of bb O
		}
		if( resn1=="GLY" ){
			frame_atomno_1 = 5; // H
			frame_atomno_2 = 1; // N
			frame_atomno_3 = 2; // CA
		}
		// runtime_assert( frame_atomno_1 > 4 );
	}
	if( fix_acceptor ){
		core::chemical::ResidueTypeSetCAP rts = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");
		core::chemical::ResidueType const & rtype = rts.lock()->name_map(resn2);
		frame_atomno_1 = rtype.nheavyatoms()-2;
		frame_atomno_2 = rtype.nheavyatoms()-1;
		frame_atomno_3 = rtype.nheavyatoms()-0;
		fixed_is_nonstandard = ( exemplars.find(resn2)!=exemplars.end() );
		cout << "alignment atomnos " << resn2 << " " << frame_atomno_1 << " " << frame_atomno_2 << " " << frame_atomno_3 << endl;
		if( resn2=="SER" ){
			frame_atomno_1 = 11; // SER is special case because so small... use HG instead of bb O
		}
		if( resn2=="GLY" ){
			frame_atomno_1 = 4; // O
			frame_atomno_2 = 3; // C
			frame_atomno_3 = 2; // CA
		}
		if( resn2=="ADX" || resn2=="CYX" || resn2=="GUX" || resn2=="THX" ){
			frame_atomno_1 = rtype.atom_index("OP1");
			frame_atomno_2 = rtype.atom_index("P");
			frame_atomno_3 = rtype.atom_index("OP2");
		}
		// runtime_assert( frame_atomno_1 > 4 );
	}
	if( fix_donor && fix_acceptor ) utility_exit_with_message("can't fix donor and acceptor");
	core::Size ifixed = fix_donor ? 1 : 2;
	core::Size imove  = fix_donor ? 2 : 1;



	omp_set_lock( &init_lock );
	gen.init(
		rot_index,
		resn1,
		resn2,
		opts.tip_tol_deg,
		opts.rot_samp_resl,
		opts.rot_samp_range,
		fix_donor,
		fix_acceptor,
		core::id::AtomID(frame_atomno_1, ifixed),
		core::id::AtomID(frame_atomno_2, ifixed),
		core::id::AtomID(frame_atomno_3, ifixed),
		exemplars
	);
	omp_unset_lock( &init_lock );

	omp_set_lock( &cout_lock );
	{
		cout << "==================================================================================================================" << endl;
		cout << "============ Thread " << omp_thread_num_1() << " generating hotspot geometries "
				<< resn1 << " " << resn2 << " =======================================" << endl;
		cout << "==================================================================================================================" << endl;
	}
	omp_unset_lock( &cout_lock );

	while( gen.has_more_samples() ){
		++totcount;
		if(totcount%100000==99999){
			omp_set_lock(&cout_lock);
			std::cout << "Thread " << omp_thread_num_1() << " Samples " << resn1 << " -> " << resn2 << " "
						 << (fix_donor?"fix_donor":"fix_acceptor") << " " << KMGT(totcount+1) << " of at most " << KMGT(gen.raw_num_samples())
						 << " progress:  " << 100.0*(double)totcount/(double)gen.raw_num_samples() << "%" << std::endl;
			omp_unset_lock(&cout_lock);
		}
		core::pose::Pose const & pose(*gen);

		// static int count = 0;
		// pose.dump_pdb( "hbgen_test_" + boost::lexical_cast<std::string>(++count) + ".pdb" );
		// if( count > 1000 ) utility_exit_with_message("check hbond gen");
		// #ifdef DEBUG
		// Xform hbonder_frame( pose.xyz(AtomID(frame_atomno_1,ifixed)), pose.xyz(AtomID(frame_atomno_2,ifixed)), pose.xyz(AtomID(frame_atomno_3,ifixed)) );
		// runtime_assert( hbonder_frame == Xform::identity() );
		// #endif
		Vec N  = pose.residue(imove).xyz(1);
		Vec CA = pose.residue(imove).xyz(2);
		Vec C  = pose.residue(imove).xyz(3);

		Vec S1 = pose.residue(ifixed).xyz(frame_atomno_1);
		Vec S2 = pose.residue(ifixed).xyz(frame_atomno_2);
		Vec S3 = pose.residue(ifixed).xyz(frame_atomno_3);

		// cout << "stored STUB " << S1 << " " << S2 << " " << S3 << endl;

		float score = gen.score_of_current();

		// if( score > 0.0 ){
		// 	pose.dump_pdb("hbond_bad_score.pdb");
		// 	utility_exit_with_message("bad score hbond generator");
		// }

		int rotamer_num = imove==1 ? gen.rot1_of_current() : gen.rot2_of_current();

		RelRotPos rel_rot_pos;
		rel_rot_pos.score = score;
		rel_rot_pos.rotamer = rotamer_num;
		for(int i = 0; i < 3; ++i){
			rel_rot_pos.n [i] =  N[i];
			rel_rot_pos.ca[i] = CA[i];
			rel_rot_pos.c [i] =  C[i];
			rel_rot_pos.stub1[i] = S1[i];
			rel_rot_pos.stub2[i] = S2[i];
			rel_rot_pos.stub3[i] = S3[i];
		}
		if( fix_acceptor ){
			rel_rot_pos.hbonding_atom_num       = gen.donor_bases_   .at(gen.idon_) - 2 ; // 0-index and skip O
			if( resn2=="GLY" )
				rel_rot_pos.hbonding_atom_num_other = 3;
			else
				rel_rot_pos.hbonding_atom_num_other = gen.acceptor_atoms_.at(gen.iacc_) - 1;//(fixed_is_nonstandard?1:2); // skip O iff std res
		}
		if( fix_donor ){
			rel_rot_pos.hbonding_atom_num       = gen.acceptor_atoms_.at(gen.iacc_) - 2 ; // 0-index and skip O
			if( resn1=="GLY" )
				rel_rot_pos.hbonding_atom_num_other = 0;
			else
				rel_rot_pos.hbonding_atom_num_other = gen.donor_bases_   .at(gen.idon_) - 1;//(fixed_is_nonstandard?1:2); // skip O iff std res
		}


		// std::string aname = rot_index.rotamers_.at(rotamer_num).atoms_.at( rel_rot_pos.hbonding_atom_num  ).data().atomname;
		// // cout << (fix_donor?resn2+"_ACC ":resn1+"_DON ") << rot_index.rotamers_[rotamer_num].atoms_[
		//rel_rot_pos.hbonding_atom_num-1].data().atomname << " " << rel_rot_pos.hbonding_atom_num-1 << endl;
		// cout << (fix_donor?resn2+"_ACC ":resn1+"_DON ") << rot_index.rotamers_[rotamer_num].atoms_[
		//rel_rot_pos.hbonding_atom_num  ].data().atomname << " " << rel_rot_pos.hbonding_atom_num   << endl;
		// // cout << (fix_donor?resn2+"_ACC ":resn1+"_DON ") << rot_index.rotamers_[rotamer_num].atoms_[
		//rel_rot_pos.hbonding_atom_num+1].data().atomname << " " << rel_rot_pos.hbonding_atom_num+1 << endl;
		// // cout << endl;

		hbonders.push_back( rel_rot_pos );
		++gen;
	}

	if( hbonders.size() == 0 ){
		utility_exit_with_message("no hbonders generated for " + resn1 + " / " + resn2 );
	}

	delete gen_ptr; // total hack , use shared_ptr!

	omp_set_lock( &cout_lock );	
	std::cout << "total " << resn1 << " " << resn2 << " " << KMGT(totcount) << std::endl;
	omp_unset_lock( &cout_lock );

	// utility_exit_with_message("debug hbgeom");
	omp_destroy_lock( &init_lock );
	omp_destroy_lock( &cout_lock );		


 }



}
}
}