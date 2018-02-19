// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:

#include <riflib/HBondedPairGenerator.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibrary.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibraryFactory.hh>

#include <Eigen/Geometry>

namespace devel {
namespace scheme {

typedef numeric::xyzVector<core::Real> Vec;
typedef numeric::xyzMatrix<core::Real> Mat;
typedef numeric::xyzTransform<core::Real> Xform;





Xform
my_align( Vec to, Vec from )
{
	double const angle = std::acos(to.normalized().dot(from.normalized()));
	if( fabs( angle ) < 0.0001 ){
		return Xform::identity();
	} else if( fabs(angle) > 2*M_PI-0.0001 ){
		// arbitrary perp. vector
		Xform(numeric::rotation_matrix(to.cross(Vec(1,2,3)),-angle));
	}
	return Xform(numeric::rotation_matrix(to.cross(from),-angle));
}
// // broken!
// Xform
// align_fast( Vec to, Vec from )
// {
// 	// rot 180 around avg. axis
// 	return Xform( 2.0*numeric::projection_matrix<double>(to+from)-Mat::identity() );
// }




core::Real
get_rot_score(
	core::pose::Pose & pose,
	core::Size ir,
	core::pack::dunbrack::RotamerLibraryScratchSpaceOP scratch
){
	core::pack::rotamers::SingleResidueRotamerLibraryCOP rotlib =
		core::pack::rotamers::SingleResidueRotamerLibraryFactory::get_instance()->get( pose.residue(ir).type() );
	return rotlib->rotamer_energy( pose.residue(ir), pose, *scratch );
}

void dump_pdb_atom(
	std::ostream & out,
	double x, double y, double z
){
	// std::string atomname,resname,elem;
	// int atomnum,resnum;
	// char chain;
	// bool ishet;
	// float occ,bfac;
	BOOST_VERIFY( x<10000 && x > -1000 );
	BOOST_VERIFY( y<10000 && y > -1000 );
	BOOST_VERIFY( z<10000 && z > -1000 );
	// std::cout << "ATOM   1604  C   GLU A 220       5.010  12.933   1.553  1.00 41.10           C" << std::endl;
	char buf[128];
	// std::string aname = a.atomname;
	// if( aname.size() == 1 ) aname = aname+"  ";
	// if( aname.size() == 2 ) aname = aname+" ";
	snprintf(buf,128,"%s%5i %4s %3s %c%4i    %8.3f%8.3f%8.3f%6.2f%6.2f %11s\n",
		true ? "HETATM":"ATOM  ",
		1,
		"DOTS",
		"DOT",
		'A',
		1,
		x,y,z,
		1.0,
		1.0,
		"D"
	);
	out << buf;
}



///@brief generate BCC lattice "hyper-ring" to sample rotationr around
/// the hbond axis as well as tilting the axis a up to short_tol
utility::vector1<numeric::xyzMatrix<double> >
get_rotation_samples( double short_tol, double long_tol, double reslang ){
	using namespace Eigen;
	typedef Map<Matrix3d> MMap;

	short_tol = numeric::conversions::radians(short_tol);
	long_tol = numeric::conversions::radians(long_tol);
	reslang = numeric::conversions::radians(reslang);

	if( reslang == 0.0 ) return utility::vector1<numeric::xyzMatrix<double> >(1,numeric::xyzMatrix<double>::identity());

	double spacing = sqrt(3.0)*reslang;
	int const n1 = std::ceil(long_tol/spacing);
	spacing = long_tol/n1;
	int const n2 = std::ceil(short_tol/spacing);
	std::cout << "gen_hbond_rotations: spacing " << spacing << " n1 " << n1 << " n2 " << n2 << std::endl;
	std::cout << "   short_tol " << short_tol << std::endl;
	std::cout << "   long_tol  " << long_tol << std::endl;
	std::cout << "   reslang   " << reslang << std::endl;

	// utility::io::ozstream out("test.pdb");
	// generate bcc "square rod" wrapped around in 3sphere wx plane
	utility::vector1<numeric::xyzMatrix<double> > rots;
	for(int i =  0 ; i <  n1; ++i){
	for(int j = -n2; j <= n2; ++j){
	for(int k = -n2; k <= n2; ++k){
	for(int o =  0;  o <=  1; ++o){
		// if(o && ( j==n2 || k==n2 ) ) continue;
		double const wx = i*spacing + (o?spacing/2.0:0.0) - long_tol/2.0;
		double const  w = cos( wx/2.0 );
		double const  x = sin( wx/2.0 );
		double const  y = ( j*spacing + (o?spacing/2.0:0.0) ) / 2.0;
		double const  z = ( k*spacing + (o?spacing/2.0:0.0) ) / 2.0;
		if( y*y + z*z > short_tol*short_tol/4.0 ) continue;
		// double const r = 1.0+y;
		// dump_pdb_atom(out,r*w*50.0,r*x*50.0,r*z*50.0);
		Quaterniond q(w,x,y,z);
		// printf("%4d %4d %4d %8.3f %8.3f %8.3f %8.3f\n",i,j,k,wx,y,z,q.norm());
		q.normalize();
		// dump_pdb_atom(out,wx*50.0,y*50.0,z*50.0);
		numeric::xyzMatrix<double> rot;
		MMap mm((double*)&rot);
		mm = q.matrix();
		// std::cout << rot << std::endl;
		// std::cout << q.matrix() << std::endl;
		// std::cout << std::endl;
		rots.push_back(rot);
	}}}}

	// out.close();

// std::exit(0);

	return rots;
}




void HBondedPairGenerator::init(
	RotamerIndex const & rot_index,
	std::string resn1,
	std::string resn2,
	double tip_tol,
	double rot_resl,
	double rot_range,
	bool fix_donor,
	bool fix_acceptor,
	core::id::AtomID align_atom1,
	core::id::AtomID align_atom2,
	core::id::AtomID align_atom3,
	std::map<std::string,core::pose::Pose> const & exemplars
){
	std::cout << "ROTATION_SAMPLING tip_tol: " << tip_tol << " rot_range: " << rot_range << " resl: " << rot_resl << std::endl;
	rot_samples_ = get_rotation_samples(tip_tol,rot_range,rot_resl);
	fix_donor_ = fix_donor;
	fix_acceptor_ = fix_acceptor;

	scratch_ = core::pack::dunbrack::RotamerLibraryScratchSpaceOP( new core::pack::dunbrack::RotamerLibraryScratchSpace );

	// #ifdef USE_OPENMP
	// #pragma omp critical
	// #endif
	{
		// only hb score
		#pragma omp critical
		{
			// make sure this has been done in a critical block
			core::pack::rotamers::SingleResidueRotamerLibraryFactory::get_instance();

			score_func_ = core::scoring::ScoreFunctionOP( new core::scoring::ScoreFunction );
			// score_func_ = core::scoring::get_score_function();
			if( resn1=="GLY" || resn2=="GLY" ){
				score_func_->set_weight( core::scoring::hbond_bb_sc, 1.1 );
			} else {
				score_func_->set_weight( core::scoring::hbond_sc, 1.1 );
			}
			// score_func_->set_weight( core::scoring::fa_atr, 0.8 );
			// score_func_->set_weight( core::scoring::fa_rep, 0.44 );
			// score_func_->set_weight( core::scoring::fa_sol, 0.75 );
			// set use_hb_env_dep false
				core::scoring::methods::EnergyMethodOptions opts = score_func_->energy_method_options();
				core::scoring::hbonds::HBondOptions hopts = opts.hbond_options();
				// std::cout << "use_hb_env_dep " <<  hopts.use_hb_env_dep() << std::endl;
				hopts.use_hb_env_dep( false );
				opts.hbond_options( hopts );
				// std::cout << score_func_->energy_method_options().hbond_options().use_hb_env_dep() << std::endl;
				score_func_->set_energy_method_options( opts );
				// std::cout << score_func_->energy_method_options().hbond_options().use_hb_env_dep() << std::endl;
				// std::cout << core::scoring::get_score_function()->energy_method_options().hbond_options().use_hb_env_dep() << std::endl;
				// std::cout << "SCORE_FUNCTION: " << std::endl << *score_func_ << std::endl;
				// utility_exit_with_message("roitsn");
		}
		core::chemical::ResidueTypeSetCAP rts;
		#pragma omp critical
		{		
			rts = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");
		}
 		rtype1op_ = rts.lock()->name_map(resn1).clone();
		rtype2op_ = rts.lock()->name_map(resn2).clone();
		if( rtype1op_->has( "H" ) && resn1 != "GLY" ) rtype1op_->set_atom_type( "H", "VIRT" );
		if( rtype1op_->has( "O" ) && resn1 != "GLY" ) rtype1op_->set_atom_type( "O", "VIRT" );
		if( rtype2op_->has( "H" ) && resn2 != "GLY" ) rtype2op_->set_atom_type( "H", "VIRT" );
		if( rtype2op_->has( "O" ) && resn2 != "GLY" ) rtype2op_->set_atom_type( "O", "VIRT" );
 		// core::chemical::ResidueType const & rtype1 = rts.lock()->name_map(resn1);
		// core::chemical::ResidueType const & rtype2 = rts.lock()->name_map(resn2);
	}

	core::chemical::ResidueType const & rtype1 = *rtype1op_;
	core::chemical::ResidueType const & rtype2 = *rtype2op_;
	core::conformation::ResidueOP res1op, res2op;
	#pragma omp critical
	{
		res1op = core::conformation::ResidueFactory::create_residue(rtype1);
		res2op = core::conformation::ResidueFactory::create_residue(rtype2);
	}

	if( exemplars.find(resn1) == exemplars.end() ){
		pose0_.append_residue_by_jump( *res1op, 1 );
	} else {
		// #ifdef USE_OPENMP
		// #pragma omp critical
		// #endif
		std::cout << "found " << resn1 << " in exemplars"  << std::endl;
		pose0_.append_residue_by_jump( exemplars.find(resn1)->second.residue(1), 1 );
	}
	if( exemplars.find(resn2) == exemplars.end() ){
		pose0_.append_residue_by_jump( *res2op, 1 );
	} else {
		// #ifdef USE_OPENMP
		// #pragma omp critical
		// #endif
		std::cout << "found " << resn2 << " in exemplars"  << std::endl;
		pose0_.append_residue_by_jump( exemplars.find(resn2)->second.residue(1), 1 );
	}

	// core::pose::add_lower_terminus_type_to_pose_residue(pose0_,1);
	// core::pose::add_lower_terminus_type_to_pose_residue(pose0_,2);
	// core::pose::add_upper_terminus_type_to_pose_residue(pose0_,1);
	// core::pose::add_upper_terminus_type_to_pose_residue(pose0_,2);
	// pose0_.set_phi(1,64.0);
	// pose0_.set_phi(2,64.0);
	// pose0_.set_psi(1,41.0);
	// pose0_.set_psi(2,41.0);

	pose_ = pose0_;

	{
		std::pair<size_t,size_t> b1 = rot_index.index_bounds(resn1.substr(0,3)), b2 = rot_index.index_bounds(resn2.substr(0,3));
		irotindex_start1_ = b1.first;
		irotindex_start2_ = b2.first;
		if( b1.first < b1.second ){
			std::cout << "adding rotamers res1 " << resn1 << " " << b1.first << "-" << b1.second-1;
			for(size_t i = b1.first; i < b1.second; ++i){
				// irotindex_start1_=98; // for testing asn BTN issue
				// if( i != 98 ) continue;
				std::cout << " " << i;
				utility::vector1<float> chis;
				for(size_t j = 0; j < rot_index.rotamers_[i].chi_.size(); ++j)
					chis.push_back( rot_index.rotamers_[i].chi_[j] );
				rotset1_.push_back( chis );
				runtime_assert( pose_.residue(1).nchi()==rotset1_.back().size() );
				// runtime_assert( i == rotset1_.size()+irotindex_start1_-1 );
			}
			std::cout << std::endl;
		} else {
			std::cout << "no rotamers for res1 " << resn1 << std::endl;
			utility::vector1<float> chis;
			for(size_t j = 0; j < pose0_.residue(1).nchi(); ++j)
				chis.push_back( pose0_.residue(1).chi(j+1) );
			rotset1_.push_back( chis );
		}
		if( b2.first < b2.second ){
			std::cout << "adding rotamers res2 " << resn2 << " " << b2.first << "-" << b2.second-1;
			for(size_t i = b2.first; i < b2.second; ++i){
				std::cout << " " << i;
				utility::vector1<float> chis;
				for(size_t j = 0; j < rot_index.rotamers_[i].chi_.size(); ++j)
					chis.push_back( rot_index.rotamers_[i].chi_[j] );
				rotset2_.push_back( chis );
				runtime_assert( pose_.residue(2).nchi()==rotset2_.back().size() );
				runtime_assert( i == rotset2_.size()+irotindex_start2_-1 );
			}
			std::cout << std::endl;
		} else {
			std::cout << "no rotamers for res2 " << resn2 << std::endl;
			utility::vector1<float> chis;
			for(size_t j = 0; j < pose0_.residue(2).nchi(); ++j)
				chis.push_back( pose0_.residue(2).chi(j+1) );
			rotset2_.push_back( chis );
		}
	}

	nrots1_ = rotset1_.size();
	nrots2_ = rotset2_.size();
	runtime_assert( nrots1_ );
	runtime_assert( nrots2_ );

	align_atom1_ = align_atom1;
	align_atom2_ = align_atom2;
	align_atom3_ = align_atom3;

	{
		std::ostringstream oss;
		oss << "HBondedPairGenerator: align atoms: "
			<< pose_.residue(align_atom1_.rsd()).atom_name(align_atom1_.atomno()) << " "
			<< pose_.residue(align_atom2_.rsd()).atom_name(align_atom2_.atomno()) << " "
			<< pose_.residue(align_atom3_.rsd()).atom_name(align_atom3_.atomno()) << std::endl;
		for(size_t i = 1; i <= pose_.residue(align_atom1_.rsd()).natoms(); ++i){
			oss << i << " " << pose_.residue(align_atom1_.rsd()).atom_name(i) << "  ";
		}
		std::cout << oss.str() << std::endl;
	}

	clash_atom_start1_ = clash_atom_start2_ = 1;
	if( fix_donor    ){
		nrots1_ = 1;
		// std::cout << "WARNING: clash checking bb atoms for donor with -fix_donor" << std::endl;
		clash_atom_start1_ = std::min( align_atom1_.atomno(), (core::Size)6 ); // TODO: check last controlling chi or something instead
		runtime_assert( align_atom1_.rsd() == 1 );
		runtime_assert( align_atom2_.rsd() == 1 );
		runtime_assert( align_atom3_.rsd() == 1 );
	}
	if( fix_acceptor ){
		nrots2_ = 1;
		// std::cout << "WARNING: clash checking bb atoms for acceptor with -fix_acceptor" << std::endl;
		clash_atom_start2_ = std::min( align_atom1_.atomno(), (core::Size)6 ); // TODO: check last controlling chi or something instead
		runtime_assert( align_atom1_.rsd() == 2 );
		runtime_assert( align_atom2_.rsd() == 2 );
		runtime_assert( align_atom3_.rsd() == 2 );
	}

	if( fix_donor && fix_acceptor ) std::cout << "WARNING: doesn't make sense to fix both donor and acceptor" << std::endl;
	donor_atoms_ = rtype1.Hpos_polar_sc();
	if( resn1=="GLY" ) donor_atoms_ = rtype1.Hpos_polar();
	for( auto i : donor_atoms_) donor_bases_.push_back( rtype1.atom_base(i) );
	acceptor_atoms_ = rtype2.accpt_pos_sc();
	if( resn2=="GLY" || resn2=="PRO" ||
	    resn2=="ADX" || resn2=="CYX" || resn2=="GUX" || resn2=="THX" // bbonly dna stuff
	){
		acceptor_atoms_ = rtype2.accpt_pos();
	}
	acceptor_orbitals_.resize(acceptor_atoms_.size());
	std::vector< Vec > orbseenit;
	std::vector<int> acceptors_to_remove;
	for(core::Size iacc = 1; iacc <= acceptor_atoms_.size(); ++iacc){
		// #ifdef USE_OPENMP
		// #pragma omp critical
		// #endif
		std::cout << "consider orbitals for atom " << iacc << " " << rtype2.atom_name(acceptor_atoms_[iacc]) << std::endl;
		for( auto j :rtype2.bonded_orbitals(acceptor_atoms_[iacc])){
			if(
				// rtype2.orbital_type(j).orbital_enum() == core::chemical::orbitals::C_pi_sp2    ||
			    // rtype2.orbital_type(j).orbital_enum() == core::chemical::orbitals::N_pi_sp2    ||
			    rtype2.orbital_type(j).orbital_enum() == core::chemical::orbitals::N_p_sp2     ||
			    // rtype2.orbital_type(j).orbital_enum() == core::chemical::orbitals::O_pi_sp2    ||
			    rtype2.orbital_type(j).orbital_enum() == core::chemical::orbitals::O_p_sp2     ||
			    rtype2.orbital_type(j).orbital_enum() == core::chemical::orbitals::S_p_sp3     ||
			    // rtype2.orbital_type(j).orbital_enum() == core::chemical::orbitals::O_pi_sp2_bb ||
			    // rtype2.orbital_type(j).orbital_enum() == core::chemical::orbitals::O_p_sp2_bb  ||
			    rtype2.orbital_type(j).orbital_enum() == core::chemical::orbitals::O_p_sp3
			){
				bool isnew = true;
				for( auto seenorb : orbseenit ){
					if( pose_.residue(2).orbital_xyz(j).distance_squared( seenorb ) < 0.001 ){
						isnew = false;
					}
				}
				if( ! isnew ){
					// #ifdef USE_OPENMP
					// #pragma omp critical
					// #endif
					std::cout << "HBPairGen WARNING: duplicate orbital found! "
					          << rtype2.name() << " "
					          << rtype2.atom_name(acceptor_atoms_[iacc]) << " "
					          << j << " "
					          << rtype2.orbital_type(j).name()
					          << std::endl;
			        // utility_exit_with_message("dup orb.");
				} else {
					if( acceptor_orbitals_[iacc].size() >= 2 ){
						// #ifdef USE_OPENMP
						// #pragma omp critical
						// #endif
						std::cout << "HBPairGen WARNING: skipping orbitals beyond first two on " << rtype2.name() << " "
						          << rtype2.atom_name(acceptor_atoms_[iacc]) << " orbname: "
	    					      << rtype2.orbital_type(j).name()
						          <<  std::endl;
			          utility_exit_with_message("extra orbs");
					} else {
						orbseenit.push_back( pose_.residue(2).orbital_xyz(j) );
						acceptor_orbitals_[iacc].push_back(j);
					}
				}
			}
		}
		if( acceptor_orbitals_[iacc].size()==0 ){
			std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
			std::cout << "WARNING!!! no acceptor orbitals "+resn2+" atom: "+ pose_.residue(2).atom_name(acceptor_atoms_[iacc]) << std::endl;
			std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
			acceptors_to_remove.push_back(iacc);
		}
	}

	if( acceptors_to_remove.size() ) { // remove "acceptors" without any applicable orbitals
		utility::vector1<core::Size> acceptor_atoms_tmp;
		utility::vector1<utility::vector1<core::Size> > acceptor_orbitals_tmp;
		for( size_t i = 1; i <= acceptor_atoms_.size(); ++i ){
			if( std::find(acceptors_to_remove.begin(),acceptors_to_remove.end(),i) != acceptors_to_remove.end() ) continue;
			acceptor_atoms_tmp.push_back( acceptor_atoms_.at(i) );
			acceptor_orbitals_tmp.push_back( acceptor_orbitals_.at(i) );
		}
		acceptor_atoms_ = acceptor_atoms_tmp;
		acceptor_orbitals_ = acceptor_orbitals_tmp;
	}

	if( donor_atoms_.size()==0 ) utility_exit_with_message("no donor atoms");
	if( acceptor_atoms_.size()==0 ) utility_exit_with_message("no acceptor atoms");
	// #ifdef DEBUG
		// for(core::Size i = 1; i <= nrots1_; ++i){
		// 	if( pose0_.residue(1).nchi()==1 ) printf("donor chi %7.2f \n",rotset1_->rotamer(i)->chi(1));
		// 	if( pose0_.residue(1).nchi()==2 ) printf("donor chi %7.2f %7.2f \n",rotset1_->rotamer(i)->chi(1),rotset1_->rotamer(i)->chi(2));
		// 	if( pose0_.residue(1).nchi()==3 ) printf("donor chi %7.2f %7.2f %7.2f \n",
		// 		rotset1_->rotamer(i)->chi(1),rotset1_->rotamer(i)->chi(2),rotset1_->rotamer(i)->chi(3));
		// 	if( pose0_.residue(1).nchi()==4 ) printf("donor chi %7.2f %7.2f %7.2f %7.2f \n",
		// 		rotset1_->rotamer(i)->chi(1),rotset1_->rotamer(i)->chi(2),rotset1_->rotamer(i)->chi(3),rotset1_->rotamer(i)->chi(4));
		// }
		// std::cout << "acceptor nchi: " << pose0_.residue(2).nchi() << std::endl;
		// for(core::Size i = 1; i <= nrots2_; ++i){
		// 	if( pose0_.residue(2).nchi()==1 ) printf("acceptor chi %7.2f \n",rotset2_->rotamer(i)->chi(1));
		// 	if( pose0_.residue(2).nchi()==2 ) printf("acceptor chi %7.2f %7.2f \n",rotset2_->rotamer(i)->chi(1),rotset2_->rotamer(i)->chi(2));
		// 	if( pose0_.residue(2).nchi()==3 ) printf("acceptor chi %7.2f %7.2f %7.2f \n",
		// 		rotset2_->rotamer(i)->chi(1),rotset2_->rotamer(i)->chi(2),rotset2_->rotamer(i)->chi(3));
		// 	if( pose0_.residue(2).nchi()==4 ) printf("acceptor chi %7.2f %7.2f %7.2f %7.2f \n",
		// 		rotset2_->rotamer(i)->chi(1),rotset2_->rotamer(i)->chi(2),rotset2_->rotamer(i)->chi(3),rotset2_->rotamer(i)->chi(4));
		// }
		// #ifdef USE_OPENMP
		// #pragma omp critical
		// #endif
		{
			for(core::Size iacc=1; iacc<=acceptor_atoms_.size(); ++iacc)
				std::cout << "acceptor orbitals, index: " << iacc << " atom: " << rtype2.atom_name(acceptor_atoms_[iacc]) << " orbitals: " << acceptor_orbitals_[iacc] << std::endl;
			std::cout << "donor atoms: ";
			for( auto idonatom : donor_atoms_ ) std::cout << rtype1.atom_name(idonatom) << ", ";
			std::cout << std::endl;
			std::cout << "donor nchi: " << pose0_.residue(1).nchi() << std::endl;
			std::cout << "donor num rots: " << nrots1_ << std::endl;
			std::cout << "acceptor num rots: " << nrots2_ << std::endl;
			std::cout << "N hbond geom samples: " << rot_samples_.size() << std::endl;
			std::cout << "total raw samples: " << raw_num_samples() << std::endl;
		}
	// #endif
	irot1_=1; irot2_=1; idon_=1; iacc_=1; iorb_=1; ihbr_=1;
	update_chi();
	if( !update_hb_and_score() ) increment();
}




	bool
	SingleHbondedPairGenerator::update_hb_and_score(
		bool realign
	){
    	typedef numeric::xyzVector<core::Real> Vec;
    	typedef numeric::xyzTransform<core::Real> Xform;
		// if(realign) std::cout << "REALIGN_HB " << irot1_ << " " << idon_ << " " << irot2_ << " " << iacc_ << " " << iorb_ << " " << ihbr_ << std::endl;
		// else        std::cout << "UPDATE_HB  " << irot1_ << " " << idon_ << " " << irot2_ << " " << iacc_ << " " << iorb_ << " " << ihbr_ << std::endl;
		// check if already "aligned"
		// pose0_.dump_pdb("before_update_hb_and_score_test0.pdb");
		// pose_ .dump_pdb("before_update_hb_and_score_test.pdb");

		if( realign ){
			// this does not fix the strange rotamer scoring issue, where successive rotamers don't get correct scores
			// #ifdef USE_OPENMP
			// #pragma omp critical
			// #endif
			// scratch_ = core::pack::dunbrack::RotamerLibraryScratchSpaceOP( new core::pack::dunbrack::RotamerLibraryScratchSpace );
			Vec don  = Vec( pose0_.residue(1).xyz(donor_atoms_[idon_]) );
			Vec donb = Vec( pose0_.residue(1).xyz(donor_bases_[idon_]) );
			Vec acc  = Vec( pose0_.residue(2).xyz(acceptor_atoms_[iacc_]));
			Vec acco = Vec( pose0_.residue(2).orbital_xyz(acceptor_orbitals_.at(iacc_).at(iorb_)));
			// Vec d = (don-donb).normalized();
			// Vec a = (acc-acco).normalized();
			Xform xd = my_align( Vec(1,0,0), don-donb );
			Xform xa = my_align( Vec(1,0,0), acc-acco );


			// std::cout << xd*(don-donb).normalized() << std::endl;
			// std::cout << xa*(acc-acco).normalized() << std::endl;
			runtime_assert( fabs( (xd*((don-donb).normalized())).dot(Vec(1,0,0)) - 1.0 ) < 0.001 );
			runtime_assert( fabs( (xa*((acc-acco).normalized())).dot(Vec(1,0,0)) - 1.0 ) < 0.001 );			

			if( fabs((don-donb).y()) < 0.00001 && fabs((don-donb).z()) < 0.00001 ) xd = Xform::identity();
			if( fabs((acc-acco).y()) < 0.00001 && fabs((acc-acco).z()) < 0.00001 ) xa = Xform::identity();
			// std::cout << "ddon " << don-donb << std::endl;
			// std::cout << "dacc " << acc-acco << std::endl;
			// std::cout << "don:  " << don  << std::endl;
			// std::cout << "donb: " << donb << std::endl;
			// std::cout << "acc:  " << acc  << std::endl;
			// std::cout << "acco: " << acco << std::endl;
			// std::cout << "xd:   " << xd   << std::endl;
			// std::cout << "xa:   " << xa   << std::endl;
			core::scoring::motif::xform_pose( pose0_, xd, 1, 1 );
			core::scoring::motif::xform_pose( pose0_, xa, 2, 2 );
			core::scoring::motif::xform_pose( pose0_, Vec(-1.15,0,0)-pose0_.residue(1).xyz(donor_atoms_   [idon_]), 1, 1 );
			core::scoring::motif::xform_pose( pose0_, Vec( 0.65,0,0)-pose0_.residue(2).xyz(acceptor_atoms_[iacc_]), 2, 2 );
			pose_ = pose0_;


		} else {
			for(core::Size ia = 1; ia <= pose0_.residue(1).natoms(); ++ia)
				pose_.set_xyz( core::id::AtomID(ia,1) , pose0_.residue(1).xyz(ia) );
			for(core::Size ia = 1; ia <= pose0_.residue(2).natoms(); ++ia)
				pose_.set_xyz( core::id::AtomID(ia,2) , rot_samples_[ihbr_] * pose0_.residue(2).xyz(ia) );
		}

		// only carbon/carbon
		for(core::Size ia = clash_atom_start1_; ia <= pose_.residue(1).nheavyatoms(); ++ia){
			if( ia > 1 && pose_.residue(1).atom_type_index(ia) > 6 ) continue;
			// if(ia==4) continue; // skip BB O, don't know psi angle
			core::Real const lj1 = pose_.residue(1).atom_type(ia).lj_radius();
			for(core::Size ja = clash_atom_start2_; ja <= pose_.residue(2).nheavyatoms(); ++ja){
				if( ja > 1 && pose_.residue(2).atom_type_index(ja) > 6 ) continue;
				// if( ia== donor_bases_[idon_] && ja==acceptor_atoms_[iacc_] ) continue;
				core::Real const lj2 = pose_.residue(2).atom_type(ja).lj_radius();
				core::Real const clash_dis = (lj1+lj2-0.5);
				core::Real const clash_dis2 = clash_dis*clash_dis;
				if( clash_dis2 > pose_.residue(1).xyz(ia).distance_squared( pose_.residue(2).xyz(ja) ) ){
					// std::cout << sqrt(clash_dis2) << " " << pose_.residue(1).xyz(ia).distance( pose_.residue(2).xyz(ja)) << std::endl;
					return false;
				}
			}
		}

		using core::id::AtomID;

		// utility_exit_with_message("oarintseiransteionraeti");
		core::scoring::motif::xform_pose( pose_, ~Xform( pose_.xyz(align_atom1_), pose_.xyz(align_atom2_), pose_.xyz(align_atom3_) ) );

		// if( realign ){
		// 	static int count = 0;
		// 	++count;
		// 	pose0_ .dump_pdb("after_realign_rotamer_"+boost::lexical_cast<std::string>(count)+".pdb");
		// 	pose_  .dump_pdb("after_realign_to_frame_"+boost::lexical_cast<std::string>(count)+".pdb");
		// }

		try {
			score_ = score_func_->score(pose_);
			// if( nrots1_ > 1 ) score_ += 0.56/2.0 * get_rot_score(pose_,1,scratch_);  // this is bugged, successive rotamers don't get proper scores
			// if( nrots2_ > 1 ) score_ += 0.56/2.0 * get_rot_score(pose_,2,scratch_);
			// std::cout << score_func_->score(pose_) << " " << get_rot_score(pose_,1,scratch_) << " " << get_rot_score(pose_,2,scratch_) << " " << score_ << std::endl;
		} catch(...) {
			{
				std::cout << "WARNING: scoring failed!" << std::endl;
				if(realign) std::cout << "REALIGN_HB " << irot1_ << " " << idon_ << " " << irot2_ << " " << iacc_ << " " << iorb_ << " " << ihbr_ << std::endl;
				else        std::cout << "UPDATE_HB  " << irot1_ << " " << idon_ << " " << irot2_ << " " << iacc_ << " " << iorb_ << " " << ihbr_ << std::endl;
				pose_.dump_pdb("score_fail.pdb");
				utility_exit_with_message("scoring failed!!");
			}
		}
		// std::cout << "DUN 1 " << get_rot_score(pose_,1,scratch_) << std::endl;
		// std::cout << "DUN 2 " << get_rot_score(pose_,2,scratch_) << std::endl;

		return true;
    }






	bool
	BidentateHbondedPairGenerator::update_hb_and_score(
		bool realign
	){
    	typedef numeric::xyzVector<core::Real> Vec;
    	typedef numeric::xyzMatrix<core::Real> Mat;
    	typedef numeric::xyzTransform<core::Real> Xform;

   		if( iacc_!=1 || iorb_!=2 ) return false;
   		if( idon_!=4 && idon_!=5 ) return false;
   		core::Size iacc_other = 2;
   		core::Size idon_other = idon_==4 ? 1 : 3;
		// #ifdef USE_OPENMP
		// #pragma omp critical
		// #endif
   		std::cout << idon_ << " " << idon_other << " " << iacc_ << " " << iacc_other << std::endl;

		// if(realign) std::cout << "REALIGN_HB " << irot1_ << " " << idon_ << " " << irot2_ << " " << iacc_ << " " << iorb_ << " " << ihbr_ << std::endl;
		// else        std::cout << "UPDATE_HB  " << irot1_ << " " << idon_ << " " << irot2_ << " " << iacc_ << " " << iorb_ << " " << ihbr_ << std::endl;

		// check if already "aligned"
		if( realign ){
			scratch_ = core::pack::dunbrack::RotamerLibraryScratchSpaceOP(new core::pack::dunbrack::RotamerLibraryScratchSpace);
			Vec don1  = Vec( pose0_.residue(1).xyz(donor_atoms_[idon_]) );
			Vec don2  = Vec( pose0_.residue(1).xyz(donor_atoms_[idon_other]) );
			Vec donb1 = Vec( pose0_.residue(1).xyz(donor_bases_[idon_]) );
			Vec donb2 = Vec( pose0_.residue(1).xyz(donor_bases_[idon_other]) );
			Vec acc1  = Vec( pose0_.residue(2).xyz(acceptor_atoms_[iacc_]));
			Vec acc2  = Vec( pose0_.residue(2).xyz(acceptor_atoms_[iacc_other]));
			Vec acco1 = Vec( pose0_.residue(2).orbital_xyz(acceptor_orbitals_.at(iacc_     ).at(iorb_)) );
			Vec acco2 = Vec( pose0_.residue(2).orbital_xyz(acceptor_orbitals_.at(iacc_other).at(iorb_)) );

			// #ifdef USE_OPENMP
			// #pragma omp critical
			// #endif
			{
				std::cout << "DON " << don1-donb1 << std::endl;
				std::cout << "DON " << don2-donb2 << std::endl;
			}

			// std::exit(0);
			Vec const dax1 = (don1+don2-donb1-donb2).normalized();
			Vec const aax1 = (acc1+acc2-acco1-acco2).normalized();
			Mat Rd = numeric::alignVectorSets( dax1, (don1-don2).cross(dax1).normalized(), Vec(0,0,1), Vec(0,1,0) );
			Mat Ra = numeric::alignVectorSets( aax1, (acc1-acc2).cross(aax1).normalized(), Vec(0,0,1), Vec(0,1,0) );
			Xform xd = Xform( Rd, Rd*(-( don1+ don2)/2.0) + Vec(0,0,-1.050) );
			Xform xa = Xform( Ra, Ra*(-(acco1+acco2)/2.0) + Vec(0,0, 0.065) );
			core::scoring::motif::xform_pose( pose0_, xd, 1, 1 );
			core::scoring::motif::xform_pose( pose0_, xa, 2, 2 );
			// core::scoring::motif::xform_pose( pose0_, Vec(-1.05,0,0)-pose0_.residue(1).xyz(donor_atoms_   [idon_]), 1, 1 );
			// core::scoring::motif::xform_pose( pose0_, Vec( 0.65,0,0)-pose0_.residue(2).xyz(acceptor_atoms_[iacc_]), 2, 2 );
			pose_ = pose0_;
		}

		for(core::Size ia = 1; ia <= pose0_.residue(2).natoms(); ++ia){
			pose_.set_xyz( core::id::AtomID(ia,2) , rot_samples_[ihbr_] * pose0_.residue(2).xyz(ia) );
		}

		// for(core::Size ia = clash_atom_start1_; ia <= pose0_.residue(1).nheavyatoms(); ++ia){
		// 	if(ia==4) continue;
		// 	core::Real const lj1 = pose0_.residue(1).atom_type(ia).lj_radius();
		// 	for(core::Size ja = clash_atom_start2_; ja <= pose0_.residue(2).nheavyatoms(); ++ja){
		// 		if( ja==4 || ( ia== donor_bases_[idon_] && ja==acceptor_atoms_[iacc_] ) ) continue;
		// 		core::Real const lj2 = pose0_.residue(2).atom_type(ja).lj_radius();
		// 		core::Real const clash_dis = (lj1+lj2-0.5);
		// 		core::Real const clash_dis2 = clash_dis*clash_dis;
		// 		if( clash_dis2 > pose0_.residue(1).xyz(ia).distance_squared( pose0_.residue(2).xyz(ja) ) ){
		// 			// std::cout << sqrt(clash_dis2) << " " << pose0_.residue(1).xyz(ia).distance( pose0_.residue(2).xyz(ja)) << std::endl;
		// 			return false;
		// 		}
		// 	}
		// }

		score_ = score_func_->score(pose_);
		// if( nrots1_ > 1 ) score_ += get_rot_score(pose_,1,scratch_); // this is bugged, successive rotamers don't get proper scores
		// if( nrots2_ > 1 ) score_ += get_rot_score(pose_,2,scratch_);	//

		return true;
    }





}
}

