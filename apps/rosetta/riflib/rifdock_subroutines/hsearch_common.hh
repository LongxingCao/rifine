// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols

#ifndef INCLUDED_riflib_rifdock_subroutines_hsearch_common_hh
#define INCLUDED_riflib_rifdock_subroutines_hsearch_common_hh


#include <scheme/types.hh>
#include <riflib/rifdock_typedefs.hh>
#include <riflib/rifdock_subroutines/output_results.hh>

#include <unordered_map>

namespace devel {
namespace scheme {


struct HSearchData {
    int64_t & total_search_effort;
    int64_t & non0_space_size;
    std::vector< typename ScaffoldProvider::ScaffoldIndex > unique_scaffolds;
};


bool
do_an_hsearch(uint64_t start_resl, 
    std::vector< std::vector< SearchPoint > > & samples,
    RifDockData & rdd,
    HSearchData & d,
    std::string const & dump_prefix,
    double beam_multiplier = 1.00) {


    using std::cout;
    using std::endl;
    using namespace devel::scheme;
    using ObjexxFCL::format::F;
    using ObjexxFCL::format::I;







    bool search_failed = false;
    {

        for( int this_stage = 0; this_stage < samples.size(); ++this_stage )
        {
            int iresl = this_stage + start_resl;

            bool using_csts = false;
            for ( ScaffoldIndex si : d.unique_scaffolds ) {
                ScaffoldDataCacheOP sdc = rdd.scaffold_provider->get_data_cache_slow(si);
                using_csts |= sdc->prepare_contraints( rdd.target, rdd.RESLS[iresl] );
            }

            bool need_sdc = using_csts || rdd.opt.tether_to_input_position != 0;


            cout << "HSearsh stage " << iresl+1 << " resl " << F(5,2,rdd.RESLS[iresl]) << " begin threaded sampling, " << KMGT(samples[this_stage].size()) << " samples: ";
            int64_t const out_interval = samples[this_stage].size()/50;
            std::exception_ptr exception = nullptr;
            std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
            start = std::chrono::high_resolution_clock::now();
            d.total_search_effort += samples[this_stage].size();

            #ifdef USE_OPENMP
            #pragma omp parallel for schedule(dynamic,64)
            #endif
            for( int64_t i = 0; i < samples[this_stage].size(); ++i ){
                if( exception ) continue;
                try {
                    if( i%out_interval==0 ){ cout << '*'; cout.flush(); }
                    RifDockIndex const isamp = samples[this_stage][i].index;

                    ScenePtr tscene( rdd.scene_pt[omp_get_thread_num()] );
                    bool director_success = rdd.director->set_scene( isamp, iresl, *tscene );
                    if ( ! director_success ) {
                        samples[this_stage][i].score = 9e9;
                        continue;
                    }

                    if ( need_sdc ) {
                        ScaffoldIndex si = isamp.scaffold_index;
                        ScaffoldDataCacheOP sdc = rdd.scaffold_provider->get_data_cache_slow(si);

                        if( rdd.opt.tether_to_input_position ){
                            float redundancy_filter_rg = sdc->get_redundancy_filter_rg( rdd.target_redundancy_filter_rg );

                            EigenXform x = tscene->position(1);
                            x.translation() -= sdc->scaffold_center;
                            float xmag =  xform_magnitude( x, redundancy_filter_rg );
                            if( xmag > rdd.opt.tether_to_input_position_cut + rdd.RESLS[iresl] ){
                                samples[this_stage][i].score = 9e9;
                                continue;
                            } 
                        }

                        /////////////////////////////////////////////////////
                        /////// Longxing' code  ////////////////////////////
                        ////////////////////////////////////////////////////
                        if (using_csts) {
                            EigenXform x = tscene->position(1);
                            bool pass_all = true;
                            for(CstBaseOP p : sdc->csts) {
                                if (!p->apply( x )) {
                                    pass_all = false;
                                    break;
                                }
                            }
                            if (!pass_all) {
                                samples[this_stage][i].score = 9e9;
                                continue;
                            }
                        }
                    }

                    // the real rif score!!!!!!
                    samples[this_stage][i].score = rdd.objectives[iresl]->score( *tscene );// + tot_sym_score;


                } catch( std::exception const & ex ) {
                    #ifdef USE_OPENMP
                    #pragma omp critical
                    #endif
                    exception = std::current_exception();
                }
            }
            if( exception ) std::rethrow_exception(exception);
            end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed_seconds_rif = end-start;
            float rate = (double)samples[this_stage].size()/ elapsed_seconds_rif.count()/omp_max_threads();
            cout << endl;// << "done threaded sampling, partitioning data..." << endl;

            SearchPoint max_pt, min_pt;
            int64_t len = samples[this_stage].size();
            if( samples[this_stage].size() > rdd.opt.beam_size/rdd.opt.DIMPOW2 * beam_multiplier ){
                __gnu_parallel::nth_element( samples[this_stage].begin(), samples[this_stage].begin()+rdd.opt.beam_size/rdd.opt.DIMPOW2 * beam_multiplier, samples[this_stage].end() );
                len = rdd.opt.beam_size/rdd.opt.DIMPOW2 * beam_multiplier;
                min_pt = *__gnu_parallel::min_element( samples[this_stage].begin(), samples[this_stage].begin()+len );
                max_pt = *(samples[this_stage].begin()+rdd.opt.beam_size/rdd.opt.DIMPOW2* beam_multiplier);
            } else {
                min_pt = *__gnu_parallel::min_element( samples[this_stage].begin(), samples[this_stage].end() );
                max_pt = *__gnu_parallel::max_element( samples[this_stage].begin(), samples[this_stage].end() );
            }

            cout << "HSearsh stage " << iresl+1 << " complete, resl. " << F(7,3,rdd.RESLS[iresl]) << ", "
                  << " " << KMGT(samples[this_stage].size()) << ", promote: " << F(9,6,min_pt.score) << " to "
                  << F(9,6, std::min(rdd.opt.global_score_cut,max_pt.score)) << " rate " << KMGT(rate) << "/s/t " << std::endl;

            // cout << "Answer: " << ( answer_exists ? "exists" : "doesn't exist" ) << std::endl;


            bool extra_for_dump = rdd.opt.dump_x_frames_per_resl > 0 && this_stage+1 == samples.size();

            if( this_stage+1 == samples.size() && ! extra_for_dump ) break;

            uint64_t dump_every = 0;
            if (rdd.opt.dump_x_frames_per_resl > 0) {
                dump_every = std::floor( len / rdd.opt.dump_x_frames_per_resl );
                if ( rdd.opt.dump_only_best_frames ) {
                    dump_every = std::max( 1, rdd.opt.dump_only_best_stride );
                    __gnu_parallel::sort( samples[this_stage].begin(), samples[this_stage].end() );
                }
            }

            for( int64_t i = 0; i < len; ++i ){
                uint64_t isamp0 = samples[this_stage][i].index.nest_index;
                if( samples[this_stage][i].score >= rdd.opt.global_score_cut ) continue;
                if ( ! extra_for_dump ) {
                    if( iresl == 0 ) ++d.non0_space_size;
                    for( uint64_t j = 0; j < rdd.opt.DIMPOW2; ++j ){
                        uint64_t isamp = isamp0 * rdd.opt.DIMPOW2 + j;
                        samples[this_stage+1].push_back( SearchPoint(RifDockIndex(isamp, samples[this_stage][i].index.scaffold_index)) );
                    }
                }

                if ( dump_every > 0 ) {  
                    if ( (   rdd.opt.dump_only_best_frames && i < rdd.opt.dump_x_frames_per_resl) ||
                         ( ! rdd.opt.dump_only_best_frames && ( i % dump_every ) == 0 )) {
                        std::string filename = dump_prefix + boost::str( boost::format( "_resl%i_%06i.pdb.gz" ) % iresl % (i/dump_every));
                        dump_search_point( rdd, samples[this_stage][i], filename, iresl, true );
                    }
                }

            }

            if ( extra_for_dump ) break;

            if( 0 == samples[this_stage+1].size() ){
                search_failed = true;
                std::cout << "search fail, no valid samples!" << std::endl;
                break;
            }
            samples[this_stage].clear();

        }
        if( search_failed ) return false;
        std::cout << "full sort of final samples" << std::endl;
        __gnu_parallel::sort( samples.back().begin(), samples.back().end() );
    }
    if( search_failed ) return false;

    return true;


}

}}


#endif