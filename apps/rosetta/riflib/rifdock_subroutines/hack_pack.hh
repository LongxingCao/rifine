// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:


#ifndef INCLUDED_riflib_rifdock_subroutines_hack_pack_hh
#define INCLUDED_riflib_rifdock_subroutines_hack_pack_hh


#include <riflib/types.hh>
#include <riflib/rifdock_typedefs.hh>
#include <riflib/rifdock_subroutines/util.hh>

using ::scheme::make_shared;
using ::scheme::shared_ptr;

typedef int32_t intRot;


namespace devel {
namespace scheme {


void
hack_pack(
    shared_ptr< std::vector< SearchPointWithRots > > & hsearch_results_p,
    std::vector< SearchPointWithRots > & packed_results,
    RifDockData & rdd,
    int64_t total_search_effort, int64_t & npack) {


    using namespace devel::scheme;
    using std::cout;
    using std::endl;



    std::vector< SearchPointWithRots > & hsearch_results = *hsearch_results_p;


    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
    start = std::chrono::high_resolution_clock::now();

    size_t n_packsamp = 0;
    for( n_packsamp; n_packsamp < hsearch_results.size(); ++n_packsamp ){
        if( hsearch_results[n_packsamp].score > 0 ) break;
    }
    int const config = rdd.RESLS.size()-1;
    npack = std::min( n_packsamp, (size_t)(total_search_effort *
        ( rdd.opt.hack_pack_frac / (rdd.packopts.pack_n_iters*rdd.packopts.pack_iter_mult)) ) );

    packed_results.resize( npack );


    print_header( "hack-packing top " + KMGT(npack) );

    std::cout << "Building twobody tables before hack-pack" << std::endl;
    for( int ipack = 0; ipack < npack; ++ipack ) {
        ScaffoldIndex si = hsearch_results[ipack].index.scaffold_index;
        rdd.scaffold_provider->setup_twobody_tables( si );
    }

    std::cout << "packing options: " << rdd.packopts << std::endl;
    std::cout << "packing w/rif rofts ";
    int64_t const out_interval = std::max<int64_t>(1,npack/100);
    std::exception_ptr exception = nullptr;
    #ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic,64)
    #endif
    for( int ipack = 0; ipack < npack; ++ipack ){
        if( exception ) continue;
        try {

            if( ipack%out_interval==0 ){ cout << '*'; cout.flush(); }
            RifDockIndex isamp = hsearch_results[ipack].index;
            if( hsearch_results[ipack].score > rdd.opt.global_score_cut ) continue;
            packed_results[ ipack ].index = isamp;
            packed_results[ ipack ].prepack_rank = ipack;
            ScenePtr tscene = ( rdd.scene_pt[omp_get_thread_num()] );
            bool director_success = rdd.director->set_scene( isamp, rdd.RESLS.size()-1, *tscene );

            if ( ! director_success ) {
                packed_results[ ipack ].rotamers(); // this initializes it to blank
                packed_results[ ipack ].score = 9e9;
                continue;
            }

            packed_results[ ipack ].score = rdd.packing_objective->score_with_rotamers( *tscene, packed_results[ ipack ].rotamers() );

        } catch( std::exception const & ex ) {
            #ifdef USE_OPENMP
            #pragma omp critical
            #endif
            exception = std::current_exception();
        }
    }
    if( exception ) std::rethrow_exception(exception);
    end = std::chrono::high_resolution_clock::now();

    std::cout << std::endl;
    std::cout << "full sort of packed samples" << std::endl;
    __gnu_parallel::sort( packed_results.begin(), packed_results.end() );


    int to_check = std::min(1000, (int)packed_results.size());
    std::cout << "Check " << to_check << " results after hackpack" << std::endl;
    for ( int i = 0; i < to_check; i++ ) {
        SearchPointWithRots const & packed_result = packed_results[i];
        if (packed_result.rotamers().size() == 0) continue;
        ScenePtr tscene( rdd.scene_pt[omp_get_thread_num()] );
        sanity_check_hackpack( rdd, packed_result.index, packed_result.rotamers_, tscene);
    }

    std::cout << std::endl;
    std::cout << "full sort of packed samples" << std::endl;
    __gnu_parallel::sort( packed_results.begin(), packed_results.end() );

    std::chrono::duration<double> elapsed_seconds_pack = end-start;
    std::cout << "packing rate: " << (double)npack/elapsed_seconds_pack.count()                   << " iface packs per second" << std::endl;
    std::cout << "packing rate: " << (double)npack/elapsed_seconds_pack.count()/omp_max_threads() << " iface packs per second per thread" << std::endl;


}

}}


#endif