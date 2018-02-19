// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.



#ifndef INCLUDED_riflib_util_hh
#define INCLUDED_riflib_util_hh

#include <riflib/types.hh>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
#include <core/pose/Pose.hh>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <numeric/xyzTransform.hh>
#include <exception>

#include <utility/io/ozstream.fwd.hh>
#include <utility/io/izstream.fwd.hh>

#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/id/AtomID.hh>

namespace devel {
namespace scheme {



#ifdef USE_OPENMP
 #include <omp.h>
 #endif
 static core::Size omp_max_threads_1(){
	#ifdef USE_OPENMP
		return omp_get_max_threads();
	#else
		return 1;
	#endif
 }
 static core::Size omp_thread_num_1(){
	#ifdef USE_OPENMP
		return omp_get_thread_num() + 1;
	#else
		return 1;
	#endif
 }
 static core::Size omp_max_threads(){
	#ifdef USE_OPENMP
		return omp_get_max_threads();
	#else
		return 1;
	#endif
 }
 static core::Size omp_thread_num(){
	#ifdef USE_OPENMP
		return omp_get_thread_num();
	#else
		return 0;
	#endif
 }

utility::vector1<core::Size> get_res(
	std::string fname,
	core::pose::Pose const & pose,
	bool nocgp = true
 );

utility::vector1<core::Size> get_designable_positions_best_guess(
	  core::pose::Pose pose
	, bool noloops
	, bool nocpg = true
 );


std::string
get_res_list_hash( utility::vector1<core::Size> const & reslist);


template<class T>
std::string str(T const & t, core::Size N=0){
	// std::ostringstream oss;
	// oss << t;
	std::string s = boost::lexical_cast<std::string>(t);
	while(s.size()<N) s = "0"+s;
	return s;
}




void
append_pose_to_pose(
	core::pose::Pose & pose1,
	core::pose::Pose const & pose2,
	bool new_chain = true
);

/// @brief Append specified residues of pose2 to pose1.
void
append_subpose_to_pose(
	core::pose::Pose & pose1,
	core::pose::Pose const & pose2,
	core::Size start_res,
	core::Size end_res,
	bool new_chain = true
);


std::vector<int> get_rif_atype_map();



template<class Float>
Eigen::Matrix<Float,3,3>
xyz2eigen( numeric::xyzMatrix<Float> const & m ){
	Eigen::Matrix<Float,3,3> rot;
	for(int i = 0; i < 3; ++i){
		for(int j = 0; j < 3; ++j){
			rot(i,j) = m(i+1,j+1);
		}
	}
	return rot;
}

template< class Float >
numeric::xyzMatrix<Float>
eigen2xyz( Eigen::Matrix<Float,3,3> const & rot ){
	numeric::xyzMatrix<Float> m;
	for(int i = 0; i < 3; ++i){
		for(int j = 0; j < 3; ++j){
			m(i+1,j+1) = rot(i,j);
		}
	}
	return m;
}

template< class Float >
numeric::xyzTransform<Float>
eigen2xyz( Eigen::Transform<Float,3,Eigen::AffineCompact> const & xin ){
	numeric::xyzTransform<Float> x( eigen2xyz( xin.rotation() ) );
	x.t[0] = xin.translation()[0];
	x.t[1] = xin.translation()[1];
	x.t[2] = xin.translation()[2];
	return x;
}

template< class Float >
void print_eigenxform( Eigen::Transform<Float,3,Eigen::AffineCompact> const & x, std::ostream & out = std::cout ){
	out << x.rotation() << "    " << x.translation().transpose() << std::endl;
}

static void print_header( std::string s, int n=120, int n2=0 ){
	for( int j = 0; j < n2; ++j )
		for( int i = 0; i < n; ++i )
			std::cout << "="; std::cout << std::endl;
	for( int i = 0; i < n/2-(int)s.size()/2-1; ++i )
		std::cout << "=";
	std::cout << " " << s << " ";
	for( int i = 0; i < n/2-(int)s.size()+(int)s.size()/2-1; ++i )
		std::cout << "=";
	std::cout << std::endl;
	for( int j = 0; j < n2; ++j )
		for( int i = 0; i < n; ++i )
			std::cout << "="; std::cout << std::endl;
}

std::string KMGT(double const & x, int const & w=7, int const & d=3);

template<class EigenXform>
float xform_magnitude(
	EigenXform const & x,
	float rg
){
	float err_trans2 = x.translation().squaredNorm();
	float cos_theta = (x.rotation().trace()-1.0)/2.0;

	// sanity check...
	// float cos_theta_slow = cos( Eigen::AngleAxisf( x.rotation() ).angle() );
	// if( fabs( cos_theta_slow-cos_theta ) > 0.001 ){
	// 	std::cout << "ang " << Eigen::AngleAxisf( x.rotation() ).angle()*180.0/M_PI << " " << cos_theta << " " << cos_theta_slow << std::endl;
	// }

	float err_rot = sqrt( std::max( 0.0, 1.0 - cos_theta*cos_theta ) ) * rg;
	if( cos_theta < 0 ) err_rot = rg;
	float err = sqrt( err_trans2 + err_rot*err_rot );
	return err;
}

void pose_to_ala( core::pose::Pose & pose );
void pose_to_gly( core::pose::Pose & pose );
void pose_to_ala( core::pose::Pose & pose, utility::vector1<core::Size> const & res_sel );


std::string
open_for_read_on_path(
	std::vector<std::string> const & path,
	std::string fname,
	utility::io::izstream & in
);

std::string
open_for_write_on_path(
	std::vector<std::string> const & path,
	std::string fname,
	utility::io::ozstream & out,
	bool create_directorys = false
);


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// OMG! MOVE ME
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template< class Float > Float
sqr( Float const & x ) { return x*x; }




// template< class HBondRay >
// float score_hbond_rays(
// 	HBondRay const & don,
// 	HBondRay const & acc,
// 	float non_directional_fraction = 0.0 // 0.0 - 1.0
// ){
// 	Eigen::Vector3f ho = acc.horb_cen - don.horb_cen;
// 	float ho_dist = ho.norm() - 1.05 ;
// 	ho_dist = ho_dist < 0 ? ho_dist*1.5 : ho_dist; // increase dis pen if too close
// 	float const max_diff = 1.3;
// 	ho_dist = ho_dist >  max_diff ?  max_diff : ho_dist;
// 	ho_dist = ho_dist < -max_diff ? -max_diff : ho_dist;

// 	ho.normalize();
// 	float bho_dot = don.direction.dot(ho);
// 	float hoa_dot = -ho.dot(acc.direction);
	
// 	if( hoa_dot < 0 ) hoa_dot = 0;
// 	if( bho_dot < 0 ) bho_dot = 0;

// 	// if( ho_dist > max_diff ) return 0.0;
// 	// sigmoid -like shape on distance score
// 	float score = -sqr( 1.0 - sqr( ho_dist/max_diff ) );

// 	// score -= bho_dot;
// 	// score -= hoa_dot;
// 	// score /= 3.0;

// 	score *= bho_dot;
// 	score *= hoa_dot;

// 	return score;
// }



template< class HBondRay >
float score_hbond_rays(
	HBondRay const & don,
	HBondRay const & acc,
	float non_directional_fraction = 0.0 // 0.0 - 1.0
){
	float diff = ( (don.horb_cen-acc.horb_cen).norm() - 1.05 );
	diff = diff < 0 ? diff*1.5 : diff; // increase dis pen if too close
	float const max_diff = 0.8;
	diff = diff >  max_diff ?  max_diff : diff;
	diff = diff < -max_diff ? -max_diff : diff;
	// if( diff > max_diff ) return 0.0;
	// sigmoid -like shape on distance score
	float score = sqr( 1.0 - sqr( diff/max_diff ) ) * -1.0;
	assert( score <= 0.0 );
	float dirscore = -don.direction.dot( acc.direction );
	dirscore = dirscore < 0 ? 0.0 : dirscore; // is positive
	float const nds = non_directional_fraction;
	score = ( 1.0 - nds )*( score * dirscore ) + nds * score;
	return score;
}

template< class VoxelArrayPtr, class HBondRay, class RotamerIndex >
struct ScoreRotamerVsTarget {
	::scheme::shared_ptr< RotamerIndex const > rot_index_p_ = nullptr;
	std::vector<VoxelArrayPtr> target_field_by_atype_;
	std::vector< HBondRay > target_donors_, target_acceptors_;
	float hbond_weight_ = 2.0;
	float upweight_iface_ = 1.0;
	float upweight_multi_hbond_ = 0.0;
	float min_hb_quality_for_multi_ = -0.5;
	float min_hb_quality_for_satisfaction_ = -0.6;
	ScoreRotamerVsTarget(){}

	template< class Xform, class Int >
	float
	score_rotamer_v_target(
		Int const & irot,
		Xform const & rbpos,
		float bad_score_thresh = 10.0, // hbonds won't get computed if grid score is above this
		int start_atom = 0 // to score only SC, use 4... N,CA,C,CB (?)
	) const	{
		int tmp1=-12345, tmp2=-12345;
		return score_rotamer_v_target_sat( irot, rbpos, tmp1, tmp2, bad_score_thresh, start_atom );
	}

	template< class Xform, class Int >
	float
	score_rotamer_v_target_sat(
		Int const & irot,
		Xform const & rbpos,
		int & sat1,
		int & sat2,
		float bad_score_thresh = 10.0, // hbonds won't get computed if grid score is above this
		int start_atom = 0 // to score only SC, use 4... N,CA,C,CB (?)
	) const	{
		using devel::scheme::score_hbond_rays;
		assert( rot_index_p_ );
		assert( target_field_by_atype_.size() == 22 );
		float score = 0;
		typedef typename RotamerIndex::Atom Atom;
		// typedef typename RotamerIndex::Rotamer Rotamer;
		for( int iatom = start_atom; iatom < rot_index_p_->nheavyatoms(irot); ++iatom )
		{
			Atom const & atom = rot_index_p_->rotamer(irot).atoms_.at(iatom);
			typename Atom::Position pos = rbpos * atom.position();
			score += target_field_by_atype_.at(atom.type())->at( pos );
		}
		// in one test: 244m with this, 182m without... need to optimize...
		if( score < bad_score_thresh ){
			float hbscore = 0;
			int hbcount = 0;
			if( rot_index_p_->rotamer(irot).acceptors_.size() > 0 ||
				rot_index_p_->rotamer(irot).donors_   .size() > 0 )
			{
				// alloca style stack bump... dangerous... don't piss memory here...
				bool used_tgt_donor   [target_donors_   .size()];
				bool used_tgt_acceptor[target_acceptors_.size()];
				bool used_rot_donor   [rot_index_p_->rotamer(irot).donors_   .size()];
				bool used_rot_acceptor[rot_index_p_->rotamer(irot).acceptors_.size()];
				for( int i = 0; i < target_donors_   .size(); ++i ) used_tgt_donor   [i] = false;
				for( int i = 0; i < target_acceptors_.size(); ++i ) used_tgt_acceptor[i] = false;
				for( int i = 0; i < rot_index_p_->rotamer(irot).donors_   .size(); ++i ) used_rot_donor   [i] = false;
				for( int i = 0; i < rot_index_p_->rotamer(irot).acceptors_.size(); ++i ) used_rot_acceptor[i] = false;

				for( int i_hr_rot_acc = 0; i_hr_rot_acc < rot_index_p_->rotamer(irot).acceptors_.size(); ++i_hr_rot_acc )
				{
					HBondRay hr_rot_acc = rot_index_p_->rotamer(irot).acceptors_.at(i_hr_rot_acc);
					Eigen::Vector3f dirpos = hr_rot_acc.horb_cen + hr_rot_acc.direction;
					hr_rot_acc.horb_cen  = rbpos * hr_rot_acc.horb_cen;
					hr_rot_acc.direction = rbpos * dirpos - hr_rot_acc.horb_cen;
					for( int i_hr_tgt_don = 0; i_hr_tgt_don < target_donors_.size(); ++i_hr_tgt_don )
					{
						HBondRay const & hr_tgt_don = target_donors_.at(i_hr_tgt_don);
						float const thishb = score_hbond_rays( hr_tgt_don, hr_rot_acc );
						hbscore += thishb * hbond_weight_;
						if( thishb < this->min_hb_quality_for_satisfaction_ ){
							if(      sat1==-1 ) sat1 = i_hr_tgt_don;
							else if( sat2==-1 ) sat2 = i_hr_tgt_don;
						}
						if( upweight_multi_hbond_ && thishb < min_hb_quality_for_multi_ ){
							if( !used_tgt_donor[i_hr_tgt_don] && !used_rot_acceptor[i_hr_rot_acc] ) ++hbcount;
							used_tgt_donor   [i_hr_tgt_don] = true;
							used_rot_acceptor[i_hr_rot_acc] = true;
						}
					}
				}
				for( int i_hr_rot_don = 0; i_hr_rot_don < rot_index_p_->rotamer(irot).donors_.size(); ++i_hr_rot_don )
				{
					HBondRay hr_rot_don = rot_index_p_->rotamer(irot).donors_.at(i_hr_rot_don);
					Eigen::Vector3f dirpos = hr_rot_don.horb_cen + hr_rot_don.direction;
					hr_rot_don.horb_cen  = rbpos * hr_rot_don.horb_cen;
					hr_rot_don.direction = rbpos * dirpos - hr_rot_don.horb_cen;
					for( int i_hr_tgt_acc = 0; i_hr_tgt_acc < target_acceptors_.size(); ++i_hr_tgt_acc )
					{
						HBondRay const & hr_tgt_acc = target_acceptors_.at(i_hr_tgt_acc);
						float const thishb = score_hbond_rays( hr_rot_don, hr_tgt_acc );
						hbscore += thishb * hbond_weight_;
						if( thishb < this->min_hb_quality_for_satisfaction_ ){
							if(      sat1==-1 ) sat1 = i_hr_tgt_acc + target_donors_.size();
							else if( sat2==-1 ) sat2 = i_hr_tgt_acc + target_donors_.size();
						}
						if( upweight_multi_hbond_ && thishb < min_hb_quality_for_multi_ ){
							if( !used_rot_donor[i_hr_rot_don] && !used_tgt_acceptor[i_hr_tgt_acc] ) ++hbcount;
							used_tgt_acceptor[i_hr_tgt_acc] = true;
							used_rot_donor[i_hr_rot_don] = true;
						}
					}
				}

			}
			// oh god, fix me..... what should the logic be??? probably "softer" thresh on thishb to count
			if( upweight_multi_hbond_ != 0 ){
				if( rot_index_p_->resname(irot)!="TYR" && // hack to aleviate my OH problems...
					rot_index_p_->resname(irot)!="SER" &&
					rot_index_p_->resname(irot)!="THR"
				){
					float multihb = 0;
					int nchi = rot_index_p_->nchi(irot) - rot_index_p_->nprotonchi(irot);
					switch( nchi ){
						case 0:
						case 1:
							if( hbcount <= 1 ) hbscore *= 1.0;
							if( hbcount > 1 ) multihb = 0.7*(hbcount-1);
							break;
						case 2:
							if( hbcount <= 1 ) hbscore *= 0.9;
							if( hbcount > 1 ) multihb = 1.0*(hbcount-1);
							break;
						case 3:
							if( hbcount <= 1 ) hbscore *= 0.7;
							if( hbcount > 1 ) multihb = 0.8*(hbcount-1);
							break;
						default:
							if( hbcount <= 1 ) hbscore *= 0.5;
							if( hbcount > 1 ) multihb = 0.7*(hbcount-1);
							break;
					}
					multihb = std::max( 0.0f, multihb );
					hbscore += hbscore * multihb * upweight_multi_hbond_; // should multihb be additive or multiplicative?
				}
			}
			if( hbscore < 0 ) score += hbscore;
		}
		return score * upweight_iface_;
	}
};



/// @brief Returns an map mapping the CA of the atoms given in the subset
/// to themselves
std::map< core::id::AtomID, core::id::AtomID >
residue_subset_to_CA_atom_map( 
    core::select::residue_selector::ResidueSubset const & subset, 
    core::pose::Pose const & pose );


/// @brief Returns the CA rmsd of the residues in the ResidueSubset
/// with or without superposition
core::Real
subset_CA_rmsd(
    core::pose::Pose const & pose1,
    core::pose::Pose const & pose2,
    core::select::residue_selector::ResidueSubset const & subset,
    bool superimpose);





////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

}
}

#endif
