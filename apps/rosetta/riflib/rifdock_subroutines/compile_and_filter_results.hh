// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:


#ifndef INCLUDED_riflib_rifdock_subroutines_compile_and_filter_results_hh
#define INCLUDED_riflib_rifdock_subroutines_compile_and_filter_results_hh


#include <riflib/types.hh>
#include <riflib/rifdock_subroutines/util.hh>
#include <unordered_map>


using ::scheme::make_shared;
using ::scheme::shared_ptr;


namespace devel {
namespace scheme {

typedef int32_t intRot;

template<class EigenXform, class ScaffoldIndex>
struct tmplXRtriple {
    EigenXform xform;
    ScaffoldIndex scaffold_index;
    uint64_t result_num;
};


// how can I fix this??? make the whole prototype into a class maybe???
// what does it do?
//  set and rescore scene with nopackscore, record more score detail
//  compute dist0
//  select results with some redundancy filtering
template<
    class EigenXform,
    // class Scene,
    class ScenePtr,
    class ObjectivePtr,
    class ScaffoldIndex
>
void
awful_compile_output_helper(
    int64_t isamp,
    int resl,
    std::vector< SearchPointWithRots > const & packed_results,
    std::vector< ScenePtr > & scene_pt,
    DirectorBase director,
    float redundancy_filter_rg,
    float redundancy_filter_mag,
    Eigen::Vector3f scaffold_center,
    std::vector< std::vector< RifDockResult > > & allresults_pt,
                 std::vector< RifDockResult >   & selected_results,
    std::vector< tmplXRtriple<EigenXform, ScaffoldIndex> > & selected_xforms,
    int n_pdb_out,
    #ifdef USE_OPENMP
        omp_lock_t & dump_lock,
    #endif
    ObjectivePtr objective,
    int & nclose,
    int nclosemax,
    float nclosethresh,
    EigenXform scaffold_perturb
) {

    typedef tmplXRtriple<EigenXform, ScaffoldIndex> XRtriple;

    SearchPointWithRots const & sp = packed_results[isamp];
    if( sp.score >= 0.0f ) return;
    ScenePtr scene_minimal( scene_pt[omp_get_thread_num()] );
    director->set_scene( sp.index, resl, *scene_minimal );
    std::vector<float> sc = objective->scores(*scene_minimal);
    float const nopackscore = sc[0]+sc[1]; //result.sum();
    float const rifscore = sc[0]; //result.template get<MyScoreBBActorRIF>();
    float const stericscore = sc[1]; //result.template get<MyClashScore>();

    // dist0 is only important to the nclose* options
    float dist0; {
        EigenXform x = scene_minimal->position(1);
        x = scaffold_perturb * x;
        x.translation() -= scaffold_center;
        dist0 = ::devel::scheme::xform_magnitude( x, redundancy_filter_rg );
    }

    RifDockResult r; // float dist0, packscore, nopackscore, rifscore, stericscore;
    r.isamp = isamp;
    r.prepack_rank = sp.prepack_rank;
    r.scene_index = sp.index;
    r.packscore = sp.score;
    r.nopackscore = nopackscore;
    r.rifscore = rifscore;
    r.stericscore = stericscore;
    r.dist0 = dist0;
    r.cluster_score = 0.0;
    r.pose_ = sp.pose_;
    allresults_pt.at( omp_get_thread_num() ).push_back( r ); // recorded w/o rotamers here

    bool force_selected = ( dist0 < nclosethresh && ++nclose < nclosemax ); // not thread-safe... is this important?

    if( selected_xforms.size() < n_pdb_out || force_selected ){

        EigenXform xposition1 = scene_minimal->position(1);
        EigenXform xposition1inv = xposition1.inverse();

        float mindiff_candidate = 9e9;
        int64_t i_closest_result;
        BOOST_FOREACH( XRtriple const & xrp, selected_xforms ){
            EigenXform const & xsel = xrp.xform;
            EigenXform const xdiff = xposition1inv * xsel;
            float diff = devel::scheme::xform_magnitude( xdiff, redundancy_filter_rg );
            if( diff < mindiff_candidate ){
                mindiff_candidate = diff;
                i_closest_result = xrp.result_num;
            }
            // todo: also compare AA composition of rotamers
        }

        if( mindiff_candidate < redundancy_filter_mag ){ // redundant result
            selected_results[i_closest_result].cluster_score += 1.0; //sp.score==0.0 ? nopackscore : sp.score;
        }

        if( mindiff_candidate > redundancy_filter_mag || force_selected ){

            #ifdef USE_OPENMP
            omp_set_lock( &dump_lock );
            #endif
            {
                // std::cout << "checking again to add selected " << selected_xforms.size() << " " << omp_get_thread_num() << std::endl;
                float mindiff_actual = 9e9;
                BOOST_FOREACH( XRtriple const & xrp, selected_xforms ){
                    EigenXform const & xsel = xrp.xform;
                    EigenXform const xdiff = xposition1inv * xsel;
                    mindiff_actual = std::min( mindiff_actual, devel::scheme::xform_magnitude( xdiff, redundancy_filter_rg ) );
                }
                if( mindiff_actual > redundancy_filter_mag || force_selected ){
                    if( redundancy_filter_mag > 0.0001 ) {
                        selected_xforms.push_back( XRtriple {
                            xposition1, 
                            sp.index.scaffold_index,
                            (int64_t)selected_results.size()
                        } );
                    }
                    r.rotamers_ = sp.rotamers_;
                    selected_results.push_back( r ); // recorded with rotamers here
                } else {
                    // std::cout << " second check failed" << std::endl;
                }
            }
            #ifdef USE_OPENMP
            omp_unset_lock( &dump_lock );
            #endif

        } // end if( mindiff > redundancy_filter_mag ){

    } // end    if( selected_xforms.size() < n_pdb_out || force_selected )

}




/*
 this needs to get fixed
 the job of this code:
 (1) do redundancy filtering
 (2) build selected_results, allresults from packed_results
 could probably split these up?
       packed_results is
            struct SearchPointWithRots {
                float score;
                uint32_t prepack_rank;
                uint64_t index;
                shared_ptr< std::vector< std::pair<intRot,intRot> > > rotamers_;
                core::pose::PoseOP pose_ = nullptr;
        allresults is
            struct RifDockResult {
                float dist0, packscore, nopackscore, rifscore, stericscore;
                uint64_t isamp, scene_index;
                uint32_t prepack_rank;
                float cluster_score;
                bool operator< ( RifDockResult const & o ) const { return packscore < o.packscore; }
                shared_ptr< std::vector< std::pair<intRot,intRot> > > rotamers_;
                core::pose::PoseOP pose_ = nullptr;
                size_t numrots() const { if(rotamers_==nullptr) return 0; return rotamers_->size(); }
                std::vector< std::pair<intRot,intRot> > const & rotamers() const { assert(rotamers_!=nullptr); return *rotamers_; }


*/

// template<class DirectorBase, class ScaffoldProvider>
// struct CompileAndFilterResultsData {
//     RifDockOpt & opt;
//     std::vector<float> & RESLS;
//     std::vector< devel::scheme::ScenePtr > & scene_pt;
//     DirectorBase & director;
//     float & target_redundancy_filter_rg;
//     // Eigen::Vector3f & scaffold_center;
//     omp_lock_t & dump_lock;
//     std::vector< devel::scheme::ObjectivePtr > & objectives;
//     // devel::scheme::EigenXform & scaffold_perturb;
//     shared_ptr<ScaffoldProvider> scaffold_provider;
// };


void compile_and_filter_results( 
        std::vector< SearchPointWithRots > & packed_results,
        std::vector< RifDockResult > & selected_results, 
        std::vector< RifDockResult > & allresults,
        RifDockData & rdd,
        omp_lock_t & dump_lock ) {


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

    typedef _RifDockResult<DirectorBase> RifDockResult;
    typedef typename ScaffoldProvider::ScaffoldIndex ScaffoldIndex;
    typedef tmplXRtriple<EigenXform, ScaffoldIndex> XRtriple;



    int64_t Nout = packed_results.size(); //samples.back().size(); //std::min(samples.back().size(),(size_t)1000000);
    Nout = std::min( (int64_t)rdd.opt.n_result_limit, Nout );

    std::vector< std::vector< RifDockResult > > allresults_pt( omp_max_threads() );

    //////////////////// brian
    // old
    // std::vector< XRtriple > selected_xforms;
    // selected_xforms.reserve(65536); // init big to reduce liklihood of resizes
    // float redundancy_filter_mag = rdd.opt.redundancy_filter_mag;
    // int nclose = 0;
    // new
    std::unordered_map< ScaffoldIndex, std::vector< XRtriple > > selected_xforms_map;  // default value here needs to be .reserve(65536)
    std::unordered_map< ScaffoldIndex, int > nclose_map; // default value here needs to be 0

    for ( uint64_t isamp = 0; isamp < Nout; isamp++ ) {
        ScaffoldIndex si = packed_results[isamp].index.scaffold_index;
        if ( selected_xforms_map.count(si) == 0 ) {
            selected_xforms_map[ si ].reserve(65536); // init big to reduce liklihood of resizes
            nclose_map[ si ] = 0;
        }
    }

    ////////////////////


    int nclosemax      = rdd.opt.force_output_if_close_to_input_num;
    float nclosethresh = rdd.opt.force_output_if_close_to_input;
    int n_pdb_out = rdd.opt.n_pdb_out;
    float redundancy_filter_mag = rdd.opt.redundancy_filter_mag;
    std::cout << "redundancy_filter_mag " << redundancy_filter_mag << "A \"rmsd\"" << std::endl;
    int64_t Nout_singlethread = std::min( (int64_t)10000, Nout );

    std::cout << "going throuth 10K results (1 thread): ";
    int64_t out_interval = 10000/81;
    for( int64_t isamp = 0; isamp < Nout_singlethread; ++isamp ){
        if( isamp%out_interval==0 ){ cout << '*'; cout.flush(); }

        ScaffoldIndex si = packed_results[isamp].index.scaffold_index;
        std::vector< XRtriple > & selected_xforms = selected_xforms_map.at( si );
        int & nclose = nclose_map.at( si );
        ScaffoldDataCacheOP sdc = rdd.scaffold_provider->get_data_cache_slow(si);
        float redundancy_filter_rg = sdc->get_redundancy_filter_rg( rdd.target_redundancy_filter_rg );
        EigenXform scaffold_perturb = sdc->scaffold_perturb;
        Eigen::Vector3f scaffold_center = sdc->scaffold_center;
                                        
        awful_compile_output_helper< EigenXform, ScenePtr, ObjectivePtr >(
            isamp, rdd.RESLS.size()-1, packed_results, rdd.scene_pt, rdd.director,
            redundancy_filter_rg, redundancy_filter_mag, scaffold_center,
            allresults_pt, selected_results, selected_xforms, n_pdb_out,
            #ifdef USE_OPENMP
                dump_lock,
            #endif
            rdd.objectives.back(), nclose, nclosemax, nclosethresh,
            scaffold_perturb
        );
    }
    std::cout << std::endl;

    std::cout << "going throuth all results (threaded): ";
    out_interval = Nout / 82;
    std::exception_ptr exception = nullptr;
    #ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic,8)
    #endif
    for( int64_t isamp = Nout_singlethread; isamp < Nout; ++isamp ){
        if( exception ) continue;
        try{
            if( isamp%out_interval==0 ){ cout << '*'; cout.flush(); }

            ScaffoldIndex si = packed_results[isamp].index.scaffold_index;
            std::vector< XRtriple > & selected_xforms = selected_xforms_map.at( si );
            int & nclose = nclose_map.at( si );
            ScaffoldDataCacheOP sdc = rdd.scaffold_provider->get_data_cache_slow(si);
            float redundancy_filter_rg = sdc->get_redundancy_filter_rg( rdd.target_redundancy_filter_rg );
            EigenXform scaffold_perturb = sdc->scaffold_perturb;
            Eigen::Vector3f scaffold_center = sdc->scaffold_center;

            awful_compile_output_helper< EigenXform, ScenePtr, ObjectivePtr >(
                isamp, rdd.RESLS.size()-1, packed_results, rdd.scene_pt, rdd.director,
                redundancy_filter_rg, redundancy_filter_mag, scaffold_center,
                allresults_pt, selected_results, selected_xforms, n_pdb_out,
                #ifdef USE_OPENMP
                    dump_lock,
                #endif
                rdd.objectives.back(), nclose, nclosemax, nclosethresh,
                scaffold_perturb
            );
        } catch(...) {
            #pragma omp critical
            exception = std::current_exception();
        }
    }
    if( exception ) std::rethrow_exception(exception);
    std::cout << std::endl;

    std::cout << "sort compiled results" << std::endl;
    BOOST_FOREACH( std::vector<RifDockResult> const & rs, allresults_pt ){
        BOOST_FOREACH( RifDockResult const & r, rs ){
            allresults.push_back( r );
        }
    }
    __gnu_parallel::sort( allresults.begin(), allresults.end() );

}

}}

#endif