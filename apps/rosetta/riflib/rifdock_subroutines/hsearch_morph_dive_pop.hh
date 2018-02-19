// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:


#ifndef INCLUDED_riflib_rifdock_subroutines_hsearch_morph_dive_pop_hh
#define INCLUDED_riflib_rifdock_subroutines_hsearch_morph_dive_pop_hh


#include <riflib/types.hh>
#include <riflib/rifdock_typedefs.hh>
#include <riflib/rifdock_subroutines/util.hh>
#include <riflib/rifdock_subroutines/hsearch_common.hh>

#include <core/import_pose/import_pose.hh>


using ::scheme::make_shared;
using ::scheme::shared_ptr;

typedef int32_t intRot;


namespace devel {
namespace scheme {

bool
hsearch_morph_dive_pop(
    RifDockData & rdd,
    HSearchData & d,
    shared_ptr<std::vector< SearchPointWithRots > > & hsearch_results_p
    ) {


    using namespace core::scoring;
        using std::cout;
        using std::endl;
        using namespace devel::scheme;
        typedef numeric::xyzVector<core::Real> Vec;
        typedef numeric::xyzMatrix<core::Real> Mat;
        // typedef numeric::xyzTransform<core::Real> Xform;
        using ObjexxFCL::format::F;
        using ObjexxFCL::format::I;
        using devel::scheme::print_header;
        using ::devel::scheme::RotamerIndex;

    typedef ::scheme::util::SimpleArray<3,float> F3;
    typedef ::scheme::util::SimpleArray<3,int> I3;



using ::scheme::scaffold::BOGUS_INDEX;
using ::scheme::scaffold::TreeIndex;
using ::scheme::scaffold::TreeLimits;


    // BOOST_FOREACH( ScenePtr & s, rdd.scene_pt ) s = rdd.scene_minimal->clone_specific_deep(std::vector<uint64_t> {1});


    shared_ptr<MorphingScaffoldProvider> morph_provider = std::dynamic_pointer_cast<MorphingScaffoldProvider>(rdd.scaffold_provider);

    ScaffoldDataCacheOP sdc = morph_provider->get_data_cache_slow(TreeIndex(0, 0));
    sdc->setup_onebody_tables( rdd.rot_index_p, rdd.opt);


    runtime_assert( rdd.opt.dive_resl <= rdd.RESLS.size() );
    std::vector< std::vector< SearchPoint > > samples( rdd.opt.dive_resl );
    samples[0].resize( rdd.director->size(0, RifDockIndex()).nest_index );

    uint64_t index_count = 0;
    for( uint64_t i = 0; i < rdd.director->size(0, RifDockIndex()).nest_index; ++i ) {
        samples[0][index_count++] = SearchPoint( RifDockIndex( i, TreeIndex(0, 0)) );
    }
    

    d.unique_scaffolds.resize(1);
    d.unique_scaffolds[0] = TreeIndex(0, 0);

    bool success = do_an_hsearch( 0, samples, rdd, d, rdd.opt.dump_prefix + "_" + sdc->scafftag + "_dp0" );

    d.unique_scaffolds.resize(0);

    if ( ! success ) return false;


    runtime_assert( rdd.opt.pop_resl <= rdd.opt.dive_resl );
    int dropped_resls = rdd.opt.dive_resl - rdd.opt.pop_resl;
    int shift_factor = dropped_resls * 6;

    std::unordered_map<uint64_t, bool> uniq_positions;

    for ( SearchPoint sp : samples.back() ) {
        uniq_positions[sp.index.nest_index >> shift_factor] = true;
    }

    std::vector<uint64_t> usable_positions;
    if ( rdd.opt.match_this_pdb == "" ) {
        for( std::pair<uint64_t, bool> const & pair : uniq_positions ) {
            usable_positions.push_back( pair.first );
        }
    } else {
        core::pose::Pose match_this = *core::import_pose::pose_from_file( rdd.opt.match_this_pdb );
        ::devel::scheme::pose_to_ala( match_this );
        Eigen::Vector3f match_center = pose_center(match_this,*(sdc->scaffold_res_p));


        utility::vector1<core::Size> target_res {1};
        std::vector< ::scheme::actor::Atom< Eigen::Vector3f > > match_atoms;
        std::vector< ::scheme::actor::Atom< Eigen::Vector3f > > scaff_atoms;
        
        devel::scheme::get_scheme_atoms( match_this, target_res, match_atoms, true );
        devel::scheme::get_scheme_atoms( *(sdc->scaffold_centered_p), target_res, scaff_atoms, true );

        EigenXform match_x = ::scheme::chemical::make_stub<EigenXform>(
                                                                    match_atoms[0].position(),
                                                                    match_atoms[1].position(),
                                                                    match_atoms[2].position());
        EigenXform scaff_x = ::scheme::chemical::make_stub<EigenXform>(
                                                                    scaff_atoms[0].position(),
                                                                    scaff_atoms[1].position(),
                                                                    scaff_atoms[2].position());

        EigenXform scaff2match = match_x * scaff_x.inverse();
        scaff2match.translation() = match_center;// - sdc->scaffold_center; // scaffold by definition is at the origin

        double error = (scaff2match * scaff_atoms[0].position() - match_atoms[0].position()).norm();
	std::cout << "Alignment error :" << error << std::endl;
        runtime_assert( error < 1 );

////////// test
        utility::vector1<core::Size> test_target_res {10};
        std::vector< ::scheme::actor::Atom< Eigen::Vector3f > > test_match_atoms;
        std::vector< ::scheme::actor::Atom< Eigen::Vector3f > > test_scaff_atoms;

        devel::scheme::get_scheme_atoms( match_this, test_target_res, test_match_atoms, true );
        devel::scheme::get_scheme_atoms( *(sdc->scaffold_centered_p), test_target_res, test_scaff_atoms, true );


        double test_error = (scaff2match * test_scaff_atoms[0].position() - test_match_atoms[0].position()).norm();
        std::cout << "Test Alignment error :" << error << std::endl;
        runtime_assert( test_error < 1 );
////////////////////////


        float redundancy_filter_rg = sdc->get_redundancy_filter_rg( rdd.target_redundancy_filter_rg );

	int count = 0;
        for( std::pair<uint64_t, bool> const & pair : uniq_positions ) {


            rdd.director->set_scene( RifDockIndex( pair.first, TreeIndex(0, 0)), rdd.opt.pop_resl-1, *rdd.scene_minimal );
            EigenXform x = rdd.scene_minimal->position(1);
            EigenXform xdiff = scaff2match.inverse() * x;
            float xmag =  xform_magnitude( xdiff, redundancy_filter_rg );
            // if (count++ < 10000) {
            // std::cout << xmag << " " << "  Trans: " 
            //   << F(7, 1, xdiff.translation()[0]) 
            //   << F(7, 1, xdiff.translation()[1]) 
            //   << F(7, 1, xdiff.translation()[2]) << std::endl; 
            // }
            if ( xmag < rdd.opt.match_this_rmsd ) {
                usable_positions.push_back( pair.first );
            }
        }



    }


    morph_provider->test_make_children( TreeIndex(0, 0) );

    TreeLimits limits = morph_provider->get_scaffold_index_limits();
    if (limits.size() == 1) {
        std::cout << "No scaffolds generated for morph_rules selection" << std::endl;
        return false;
    }

    uint64_t num_scaffolds = limits[1];

    for ( uint64_t scaffno = 0; scaffno < num_scaffolds; scaffno++ ) {
        TreeIndex ti(1, scaffno);
        ScaffoldDataCacheOP sdc = morph_provider->get_data_cache_slow(ti);
        // some options allow one to skip generating these here
        if ( ! sdc->local_onebody_p ) {
            sdc->setup_onebody_tables( rdd.rot_index_p, rdd.opt);
        }

        d.unique_scaffolds.push_back(ti);

    }


    std::cout << "Num unique positions: " << uniq_positions.size() << std::endl;
    std::cout << "Num usable positions: " << usable_positions.size() << std::endl;


    std::vector< std::vector< SearchPoint > > samples2( rdd.RESLS.size() - rdd.opt.pop_resl + 1 );
    samples2[0].resize( usable_positions.size()*num_scaffolds );


    uint64_t index_count2 = 0;
    for ( uint64_t scaffno = 0; scaffno < num_scaffolds; scaffno++ ) {
        for( uint64_t position : usable_positions ) {
            samples2[0][index_count2++] = SearchPoint( RifDockIndex( position, TreeIndex(1, scaffno)) );
        }
    }


    success = do_an_hsearch( rdd.opt.pop_resl-1, samples2, rdd, d, rdd.opt.dump_prefix + "_" + sdc->scafftag + "_dp1", 
        std::min((double)rdd.opt.max_beam_multiplier, (double)num_scaffolds ));

    if ( ! success ) return false;



    // for ( uint64_t scaffno = 0; scaffno < num_scaffolds; scaffno++ ) {
    //     TreeIndex ti(1, scaffno);
    //     ScaffoldDataCacheOP sdc = morph_provider->get_data_cache_slow(ti);
    //     sdc->setup_onebody_tables( rdd.rot_index_p, rdd.opt);
    // }






    std::cout << "total non-0 space size was approx " << float(d.non0_space_size)*1024.0*1024.0*1024.0 << " grid points" << std::endl;
    std::cout << "total search effort " << KMGT(d.total_search_effort) << std::endl;


    hsearch_results_p = make_shared<std::vector< SearchPointWithRots >>();
    std::vector< SearchPointWithRots > & hsearch_results = *hsearch_results_p;


    hsearch_results.resize( samples2.back().size() );
    #ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic,1024)
    #endif
    for( int ipack = 0; ipack < hsearch_results.size(); ++ipack ){
        hsearch_results[ipack].score = samples2.back()[ipack].score;
        hsearch_results[ipack].index = samples2.back()[ipack].index;
    }



    return true;



}


}}

#endif
