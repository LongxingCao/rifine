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

#include <rif_dock_test.hh>
#include <riflib/rotamer_energy_tables.hh>
#include <riflib/RifFactory.hh>



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






struct RifDockData {
    RifDockOpt & opt;
    std::vector<float> & RESLS;
    DirectorBase & director;
    std::vector< ScenePtr > & scene_pt;
    ScenePtr & scene_minimal;
    std::vector<SimpleAtom> & target_simple_atoms;
    std::vector< VoxelArrayPtr > & target_field_by_atype;
    std::vector< std::vector< VoxelArrayPtr > > const * target_bounding_by_atype;
    std::vector< ::scheme::chemical::HBondRay > * target_donors;
    std::vector< ::scheme::chemical::HBondRay > * target_acceptors;
    float & target_redundancy_filter_rg;
    core::pose::Pose & target;
    shared_ptr< RotamerIndex > & rot_index_p;
    RotamerRFTablesManager & rotrf_table_manager;
    std::vector< ObjectivePtr > & objectives;
    ObjectivePtr & packing_objective;
    ::scheme::search::HackPackOpts & packopts;
    std::vector<shared_ptr<RifBase> > & rif_ptrs;
    RifSceneObjectiveConfig & rso_config;
    shared_ptr<RifFactory> & rif_factory;

    ScaffoldProviderOP scaffold_provider;
};


inline
void
sanity_check_rots(
    RifDockData & rdd, 
    RifDockIndex i,
    shared_ptr< std::vector< std::pair<intRot,intRot> > > rotamers,
    ScenePtr scene,
    bool original ) {

    devel::scheme::ScoreRotamerVsTarget<
        VoxelArrayPtr, ::scheme::chemical::HBondRay, ::devel::scheme::RotamerIndex
    > rot_tgt_scorer;
    rot_tgt_scorer.rot_index_p_ = rdd.rot_index_p;
    rot_tgt_scorer.target_field_by_atype_ = rdd.target_field_by_atype;
    rot_tgt_scorer.target_donors_ = *rdd.target_donors;
    rot_tgt_scorer.target_acceptors_ = *rdd.target_acceptors;
    rot_tgt_scorer.hbond_weight_ = rdd.packopts.hbond_weight;
    rot_tgt_scorer.upweight_iface_ = rdd.packopts.upweight_iface;
    rot_tgt_scorer.upweight_multi_hbond_ = rdd.packopts.upweight_multi_hbond;


    bool only_bad = true;
    bool all_missing = true;
    bool all_ala = true;

    for( int ipr = 0; ipr < rotamers->size(); ++ipr ){
        int irot = rotamers->at(ipr).second;

        BBActor bba = scene->template get_actor<BBActor>(1,rotamers->at(ipr).first);

        float rescore = rot_tgt_scorer.score_rotamer_v_target( irot, bba.position(), 10.0, 4 );
        if (rescore >= 0) {
        } else {
        }

        std::vector< std::pair< float, int > > rotscores;
        rdd.rif_ptrs.back()->get_rotamers_for_xform( bba.position(), rotscores );

        bool exists = false;
        for ( std::pair<float,int> const & p : rotscores ) {
            if (p.second == irot) {
                exists = true;
                break;
            } else {
            }
        }
        if (irot == 0) {
        } else {
            if (exists) {
                all_missing=false;
            } else {
                if (original) {
                    std::cout << "ZMISSING ROT";
                } else {
                    std::cout << "RECALCMISSING ROT";
                }
            }
            all_ala = false;
        }


    }

}



inline
void 
sanity_check_hackpack(
    RifDockData & rdd, 
    RifDockIndex i,
    shared_ptr< std::vector< std::pair<intRot,intRot> > > rotamers,
    ScenePtr scene ) {

    bool success = rdd.director->set_scene( i, rdd.RESLS.size()-1, *scene );
    if ( ! success ) {
        std::cout << "Bad index" << std::endl;
        return;
    }
    sanity_check_rots(rdd, i, rotamers, scene, true);

    rdd.director->set_scene( i, rdd.RESLS.size()-1, *scene );
    SearchPointWithRots temp;

    rdd.packing_objective->score_with_rotamers( *scene, temp.rotamers() );

    sanity_check_rots(rdd, i, temp.rotamers_, scene, false);


}



}}



#endif