// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:


#ifndef INCLUDED_riflib_rifdock_subroutines_util_hh
#define INCLUDED_riflib_rifdock_subroutines_util_hh


#include <riflib/types.hh>

#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>

#include <numeric/xyzTransform.hh>
#include <numeric/xyzVector.hh>
#include <utility/vector1.hh>
#include <scheme/kinematics/Director.hh>
#include <scheme/search/HackPack.hh>
#include <riflib/util.hh>

#include <rifdock.hh>
#include <riflib/rotamer_energy_tables.hh>
#include <riflib/RifFactory.hh>

// this is for the compile and filter results.
#include <unordered_map>



using ::scheme::make_shared;
using ::scheme::shared_ptr;

typedef int32_t intRot;

namespace devel {
namespace scheme {

inline
Eigen::Vector3f
pose_center(
    core::pose::Pose const & pose,
    utility::vector1<core::Size> const & useres = utility::vector1<core::Size>()
){
    typedef numeric::xyzVector<core::Real> Vec;
    Vec cen(0,0,0);
    int count = 0;
    for( int ir = 1; ir <= pose.size(); ++ir ) {
        if( useres.size()==0 || std::find(useres.begin(),useres.end(),ir)!=useres.end() ){
            for( int ia = 1; ia <= pose.residue(ir).nheavyatoms(); ++ia ){
                cen += pose.xyz(core::id::AtomID(ia,ir));
                ++count;
            }
        // } else {
            // std::cout << "pose_center skip " << ir << std::endl;
        }
    }
    cen /= double(count);
    // ::scheme::util::SimpleArray<3,float> center;
    Eigen::Vector3f center;
    center[0] = cen[0];
    center[1] = cen[1];
    center[2] = cen[2];
    return center;
}

inline
void
get_rg_radius(
    core::pose::Pose const & pose,
    float & rg,
    float & radius,
    utility::vector1<core::Size> const & useres = utility::vector1<core::Size>(),
    bool allatom = false
){
    Eigen::Vector3f centmp = pose_center( pose, useres );
    numeric::xyzVector<double> cen;
    float maxdis = -9e9, avgdis2 = 0.0;
    for( int i = 0; i < 3; ++i ) cen[i] = centmp[i];
    for( int iri = 1; iri <= useres.size(); ++iri ){
        int ir = useres[iri];
        if( allatom ){
            for( int ia = 1; ia <= pose.residue(ir).nheavyatoms(); ++ia ){
                numeric::xyzVector<double> coord = pose.residue(ir).xyz(ia);
                avgdis2 += cen.distance_squared( coord );
                maxdis = std::max( maxdis, (float)cen.distance( coord ) );
            }
        } else {
            numeric::xyzVector<double> coord;
            if(      pose.residue(ir).has("CB") ) coord = pose.residue(ir).xyz("CB");
            else if( pose.residue(ir).has("CA") ) coord = pose.residue(ir).xyz("CA");
            else                                  coord = pose.residue(ir).nbr_atom_xyz();
            avgdis2 += cen.distance_squared( coord );
            maxdis = std::max( maxdis, (float)cen.distance( coord ) );
        }
    }
    avgdis2 /= useres.size();
    rg = sqrt( avgdis2 );
    radius = maxdis;
}


inline
void xform_pose( core::pose::Pose & pose, numeric::xyzTransform<float> s, core::Size sres=1, core::Size eres=0 ) {
  if(eres==0) eres = pose.size();
  for(core::Size ir = sres; ir <= eres; ++ir) {
    for(core::Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia) {
      core::id::AtomID const aid(core::id::AtomID(ia,ir));
      pose.set_xyz( aid, s*pose.xyz(aid) );
    }
  }
}




template<class _DirectorBigIndex>
struct tmplRifDockResult {
    typedef _DirectorBigIndex DirectorBigIndex;
    typedef tmplRifDockResult<DirectorBigIndex> This;
    float dist0, packscore, nopackscore, rifscore, stericscore;
    uint64_t isamp;
    DirectorBigIndex scene_index;
    uint32_t prepack_rank;
    float cluster_score;
    bool operator< ( This const & o ) const { return packscore < o.packscore; }
    shared_ptr< std::vector< std::pair<intRot,intRot> > > rotamers_;
    core::pose::PoseOP pose_ = nullptr;
    size_t numrots() const { if(rotamers_==nullptr) return 0; return rotamers_->size(); }
    std::vector< std::pair<intRot,intRot> > const & rotamers() const { assert(rotamers_!=nullptr); return *rotamers_; }
};




#pragma pack (push, 4) // allows size to be 12 rather than 16
template<class _DirectorBigIndex>
struct tmplSearchPoint {
    typedef _DirectorBigIndex DirectorBigIndex;
    typedef tmplSearchPoint<DirectorBigIndex> This;
    float score;
    DirectorBigIndex index;
    tmplSearchPoint() : score(9e9) {}
    tmplSearchPoint(DirectorBigIndex i) : score(9e9), index(i) {}
    bool operator < (This const & o) const {
        return score < o.score;
    }
};
#pragma pack (pop)



template<class _DirectorBigIndex>
struct tmplSearchPointWithRots {
    typedef _DirectorBigIndex DirectorBigIndex;
    typedef tmplSearchPointWithRots<DirectorBigIndex> This;
    float score;
    uint32_t prepack_rank;
    DirectorBigIndex index;
    shared_ptr< std::vector< std::pair<intRot,intRot> > > rotamers_;
    core::pose::PoseOP pose_ = nullptr;
    tmplSearchPointWithRots() : score(9e9), prepack_rank(0), rotamers_(nullptr) {}
    tmplSearchPointWithRots(DirectorBigIndex i, uint32_t orank) : score(9e9), prepack_rank(orank), index(i), rotamers_(nullptr) {}
    // ~SearchPointWithRots() { delete rotamers_; }
    void checkinit() { if( rotamers_==nullptr ) rotamers_ = make_shared< std::vector< std::pair<intRot,intRot> > > ();  }
    std::vector< std::pair<intRot,intRot> > & rotamers() { checkinit(); return *rotamers_; }
    std::vector< std::pair<intRot,intRot> > const & rotamers() const { runtime_assert(rotamers_!=nullptr); return *rotamers_; }
    size_t numrots() const { if(rotamers_==nullptr) return 0; return rotamers_->size(); }
    bool operator < (This const & o) const {
        return score < o.score;
    }
};


// Convenience templates for the above templated containers

template <class __Director>
using _SearchPointWithRots = tmplSearchPointWithRots<_DirectorBigIndex<__Director>>;
typedef _SearchPointWithRots<DirectorBase> SearchPointWithRots;

template <class __Director>
using _RifDockResult = tmplRifDockResult<_DirectorBigIndex<__Director>>;
typedef _RifDockResult<DirectorBase> RifDockResult;

template <class __Director>
using _SearchPoint = tmplSearchPoint<_DirectorBigIndex<__Director>>;
typedef _SearchPoint<DirectorBase> SearchPoint;

    
    
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
                                std::vector< std::pair<EigenXform, uint64_t> > & selected_xforms,
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
            typedef std::pair< EigenXform, int64_t > XRpair;
            BOOST_FOREACH( XRpair const & xrp, selected_xforms ){
                EigenXform const & xsel = xrp.first;
                EigenXform const xdiff = xposition1inv * xsel;
                float diff = devel::scheme::xform_magnitude( xdiff, redundancy_filter_rg );
                if( diff < mindiff_candidate ){
                    mindiff_candidate = diff;
                    i_closest_result = xrp.second;
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
                    BOOST_FOREACH( XRpair const & xrp, selected_xforms ){
                        EigenXform const & xsel = xrp.first;
                        EigenXform const xdiff = xposition1inv * xsel;
                        mindiff_actual = std::min( mindiff_actual, devel::scheme::xform_magnitude( xdiff, redundancy_filter_rg ) );
                    }
                    if( mindiff_actual > redundancy_filter_mag || force_selected ){
                        if( redundancy_filter_mag > 0.0001 ) {
                            selected_xforms.push_back( std::make_pair(xposition1, (int64_t)selected_results.size()) );
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
    


}}



#endif
