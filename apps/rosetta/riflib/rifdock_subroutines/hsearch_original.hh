// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:


#ifndef INCLUDED_riflib_rifdock_subroutines_hsearch_original_hh
#define INCLUDED_riflib_rifdock_subroutines_hsearch_original_hh


#include <riflib/types.hh>
#include <riflib/rifdock_typedefs.hh>
#include <riflib/rifdock_subroutines/util.hh>
#include <riflib/rifdock_subroutines/hsearch_common.hh>

using ::scheme::make_shared;
using ::scheme::shared_ptr;

typedef int32_t intRot;


namespace devel {
namespace scheme {

bool
hsearch_original(
    RifDockData & rdd,
    HSearchData & d,
    shared_ptr<std::vector< SearchPointWithRots > > & hsearch_results_p
    ) {

    using std::cout;
    using std::endl;
    using namespace devel::scheme;


    ScaffoldIndex si = scaffold_index_default_value( ScaffoldIndex());
    ScaffoldDataCacheOP sdc = rdd.scaffold_provider->get_data_cache_slow(si);
    float redundancy_filter_rg = sdc->get_redundancy_filter_rg( rdd.target_redundancy_filter_rg );
    Eigen::Vector3f scaffold_center = sdc->scaffold_center;

    sdc->setup_onebody_tables( rdd.rot_index_p, rdd.opt);

    std::string dump_prefix = rdd.opt.dump_prefix + "_" + sdc->scafftag;


    std::vector< std::vector< SearchPoint > > samples( rdd.RESLS.size() );
    samples[0].resize( rdd.director->size(0, RifDockIndex()).nest_index );
    for( uint64_t i = 0; i < rdd.director->size(0, RifDockIndex()).nest_index; ++i ) samples[0][i] 
        = SearchPoint( RifDockIndex( i, scaffold_index_default_value( ScaffoldIndex())) );

    d.unique_scaffolds.resize(1);
    d.unique_scaffolds.front() = si;

    bool search_failed = !do_an_hsearch(0, samples, rdd, d, dump_prefix);

    if ( search_failed ) return false;

    std::cout << "total non-0 space size was approx " << float(d.non0_space_size)*1024.0*1024.0*1024.0 << " grid points" << std::endl;
    std::cout << "total search effort " << KMGT(d.total_search_effort) << std::endl;


    hsearch_results_p = make_shared<std::vector< SearchPointWithRots >>();
    std::vector< SearchPointWithRots > & hsearch_results = *hsearch_results_p;


    hsearch_results.resize( samples.back().size() );
    #ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic,1024)
    #endif
    for( int ipack = 0; ipack < hsearch_results.size(); ++ipack ){
        hsearch_results[ipack].score = samples.back()[ipack].score;
        hsearch_results[ipack].index = samples.back()[ipack].index;
    }



    return true;



}

}}


#endif