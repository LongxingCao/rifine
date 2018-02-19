// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols
#ifndef INCLUDED_riflib_rotamer_energy_tables_hh
#define INCLUDED_riflib_rotamer_energy_tables_hh


#include <riflib/RotamerGenerator.hh>
#include <scheme/objective/storage/TwoBodyTable.hh>
#include <scheme/objective/voxel/VoxelArray.hh>
namespace devel {
namespace scheme {




void get_onebody_rotamer_energies(
	core::pose::Pose const & scaffold,
	utility::vector1<core::Size> const & scaffold_res,
	devel::scheme::RotamerIndex const & rot_index,
	std::vector<std::vector<float> > & scaffold_onebody_rotamer_energies,
	std::vector<std::string> const & cachepath,
	std::string const & cachefile,
	bool replace_with_ala = true
);

void
compute_onebody_rotamer_energies(
	core::pose::Pose const & scaffold,
	utility::vector1<core::Size> const & scaffold_res,
	RotamerIndex const & rot_index,
	std::vector<std::vector< float > > & scaffold_onebody_rotamer_energies,
	bool replace_with_ala = true
);



struct RotamerRFOpts {
	int oversample;
	float field_resl, field_spread;
	float scale_atr;
	std::string data_dir;
	RotamerRFOpts()
		: oversample(1)
		, field_resl(-1)
		, field_spread(0)
		, scale_atr(1.0)
		, data_dir("./")
	{}
};

void get_per_rotamer_rf_tables_one(
	devel::scheme::RotamerIndex const & rot_index,
	int irot,
	RotamerRFOpts opts,
	int const N_ATYPE,
	std::vector< std::vector<
		::scheme::shared_ptr<
			::scheme::objective::voxel::VoxelArray< 3, float >
		>
	> >	& fields_by_atype_per_rot,
	bool verbose
);

void get_per_rotamer_rf_tables(
	devel::scheme::RotamerIndex const & rot_index,
	RotamerRFOpts opts,
	int const N_ATYPE,
	std::vector< std::vector<
		::scheme::shared_ptr<
			::scheme::objective::voxel::VoxelArray< 3, float >
		>
	> >	& fields_by_atype_per_rot
);


struct RotamerRFTablesManager {
	typedef ::scheme::objective::voxel::VoxelArray< 3, float > VoxelArray;

	#ifdef USE_OPENMP
	std::vector<omp_lock_t> locks_;
	#endif

	static int const N_ATYPE_ = 21;

	std::vector< std::vector<
		::scheme::shared_ptr< VoxelArray >
		> >	fields_by_atype_per_rot_;

	::scheme::shared_ptr< ::devel::scheme::RotamerIndex > rot_index_p_;

	RotamerRFOpts opts_;

	bool preinitdone_;

	RotamerRFTablesManager(
		::scheme::shared_ptr< ::devel::scheme::RotamerIndex > rot_index_p,
		RotamerRFOpts opts
	) :
		rot_index_p_(rot_index_p),
		opts_(opts),
		preinitdone_(false)
	{
		fields_by_atype_per_rot_.resize( rot_index_p->size() );
		#ifdef USE_OPENMP
			locks_.resize( rot_index_p->size() );
			for( auto & lock : locks_ ) omp_init_lock( &lock );
		#endif
	}
	~RotamerRFTablesManager(){
		#ifdef USE_OPENMP
			for( auto & lock: locks_ ) omp_destroy_lock( &lock );
		#endif
	}

	std::vector< ::scheme::shared_ptr< VoxelArray > > const &
	get_rotamer_rf_tables( int irot ) {
		if( !preinitdone_ ){ // this check turns out to be a tiny bit more efficient
		                     // than just "fields_by_atype_per_rot_[irot].size() != N_ATYPE_+1"
			runtime_assert_msg( 0 <= irot && irot < rot_index_p_->size(), "bad irot" );
			if( fields_by_atype_per_rot_[irot].size() != N_ATYPE_+1 ){
				// must initialize it
				#ifdef USE_OPENMP
				omp_set_lock( &locks_[irot] );
				#endif

				// check if somebody else already fixed it....
				if( fields_by_atype_per_rot_[irot].size() != N_ATYPE_+1 ){
					if( rot_index_p_->structural_parent(irot) != irot ){
						fields_by_atype_per_rot_[irot] = get_rotamer_rf_tables( rot_index_p_->structural_parent(irot) );
					} else {
						bool verbose = ( 0==irot );
						get_per_rotamer_rf_tables_one( *rot_index_p_, irot, opts_, N_ATYPE_, fields_by_atype_per_rot_, verbose );
					}
				}

				#ifdef USE_OPENMP
				omp_unset_lock( &locks_[irot] );
				#endif
			}
		}
		return fields_by_atype_per_rot_[irot];
	}

	void preinit_all(){
		get_per_rotamer_rf_tables( *rot_index_p_, opts_, N_ATYPE_, fields_by_atype_per_rot_ );
		preinitdone_ = true;
	}

};



struct MakeTwobodyOpts {
	float onebody_threshold;
	float distance_cut;
	float hbond_weight;
	float favorable_2body_multiplier;
	MakeTwobodyOpts()
		: onebody_threshold(2.0)
		, distance_cut(15.0)
		, hbond_weight(2.0)
		, favorable_2body_multiplier(1)
	{}
};

void
make_twobody_tables(
	core::pose::Pose const & scaffold,
	devel::scheme::RotamerIndex const & rot_index,
	std::vector<std::vector<float> > const & onebody_energies,
	RotamerRFTablesManager & rotrfmanager,
	MakeTwobodyOpts opts,
	::scheme::objective::storage::TwoBodyTable<float> & twob
);

void
get_twobody_tables(
	std::vector<std::string> const & cachepath,
	std::string const & cachefile,
	std::string & description, // in/out
	core::pose::Pose const & scaffold,
	devel::scheme::RotamerIndex const & rot_index,
	std::vector<std::vector<float> > const & onebody_energies,
	RotamerRFTablesManager & rotrfmanager,
	MakeTwobodyOpts opts,
	::scheme::objective::storage::TwoBodyTable<float> & twob
);


}}

#endif
