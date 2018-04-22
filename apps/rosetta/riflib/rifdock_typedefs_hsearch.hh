// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.



#ifndef INCLUDED_riflib_rifdock_typedefs_hh
#define INCLUDED_riflib_rifdock_typedefs_hh

#include <scheme/types.hh>

#include <scheme/actor/BackboneActor.hh>
#include <scheme/actor/VoxelActor.hh>
#include <scheme/actor/Atom.hh>
#include <scheme/kinematics/Scene.hh>
#include <scheme/scaffold/ScaffoldProviderBase.hh>
#include <scheme/nest/NEST.hh>
#include <scheme/nest/pmap/OriTransMap.hh>
#include <scheme/kinematics/Director.hh>


#include <boost/mpl/vector.hpp>

namespace devel {
namespace scheme {


struct RIFAnchor {
    RIFAnchor() {}
};

struct ScaffoldDataCache;


// Typedefs that should never change

typedef ::scheme::actor::BackboneActor<EigenXform> BBActor;

typedef ::scheme::actor::VoxelActor<EigenXform,float> VoxelActor;

typedef ::scheme::actor::SimpleAtom< Eigen::Vector3f > SimpleAtom;


// Typedefs related to the Hierarchical Search Scene

typedef ::scheme::actor::Score_Voxel_vs_Atom<
        VoxelActor,
        SimpleAtom,
        false
    > MyClashScore;

typedef boost::mpl::vector<
            BBActor,
            SimpleAtom,
            VoxelActor,
            RIFAnchor
        > ParametricSceneContainers;

typedef ::scheme::kinematics::impl::Conformation<
    ParametricSceneContainers,
    ScaffoldDataCache > ParametricSceneConformation; 

typedef shared_ptr<ParametricSceneConformation> ParametricSceneConformationOP;
typedef shared_ptr<ParametricSceneConformation const > ParametricSceneConformationCOP;


typedef ::scheme::kinematics::Scene<
        ParametricSceneConformation,
        EigenXform
    > ParametricScene;
    

// Typedefs related to the Hierarchical Search ScaffoldProvider

typedef ::scheme::scaffold::TreeScaffoldProvider<ParametricSceneConformation> ScaffoldProvider;
typedef shared_ptr<ScaffoldProvider> ScaffoldProviderOP;

typedef typename ScaffoldProvider::ScaffoldIndex ScaffoldIndex;




// If you add something to the index, you must follow these rules
// 1. Keep your item lightweight
// 2. Add your item to the hash function

struct RifDockIndex {
    uint64_t nest_index;
    uint64_t seeding_index;
    ScaffoldIndex scaffold_index;

    RifDockIndex() :
      nest_index(::scheme::kinematics::director_index_default_value()) 
      {
          seeding_index = 0;
          scaffold_index = ::scheme::scaffold::scaffold_index_default_value(scaffold_index);
      }

    
    RifDockIndex(
                 uint64_t nest_index_in,
                 ScaffoldIndex scaffold_index_in
                 ) :
    nest_index( nest_index_in ),
    seeding_index ( 0 ),
    scaffold_index( scaffold_index_in ) {}
    
    RifDockIndex(
                 uint64_t nest_index_in,
                 uint64_t seeding_index_in
                 ) :
    nest_index( nest_index_in ),
    seeding_index ( seeding_index_in )
    {
        scaffold_index = ::scheme::scaffold::scaffold_index_default_value(scaffold_index);
    }
    
    
    RifDockIndex( 
        uint64_t nest_index_in,
        uint64_t seeding_index_in,
        ScaffoldIndex scaffold_index_in
        ) : 
        nest_index( nest_index_in ),
        seeding_index ( seeding_index_in ),
        scaffold_index( scaffold_index_in ) {}

    bool operator==(RifDockIndex const & o) {
      return (
        nest_index == o.nest_index &&
        seeding_index == o.seeding_index &&
        scaffold_index == o.scaffold_index
        );
    }

};




// Typedefs related to the Hierarchical Search Director

typedef ::scheme::nest::NEST< 6,
              devel::scheme::EigenXform,
              ::scheme::nest::pmap::OriTransMap,
              ::scheme::util::StoreNothing, // do not store a transform in the Nest
              uint64_t,
              float,
              false // do not inherit from NestBase
             > NestOriTrans6D;


// typedef ::scheme::kinematics::NestDirector< NestOriTrans6D > DirectorOriTrans6D;
typedef ::scheme::kinematics::NestDirector< NestOriTrans6D, uint64_t> RifDockNestDirector;

typedef shared_ptr<::scheme::kinematics::Director<EigenXform, uint64_t> > DirectorBase;





}
}

namespace std {

    template <>
    struct hash<devel::scheme::RifDockIndex>
    {
        std::size_t operator()(const devel::scheme::RifDockIndex& rdi) const {
            using std::size_t;
            using boost::hash;
            using boost::hash_combine;

            std::size_t seed = 0;

            boost::hash<int> hasher;
            hash_combine(seed, hasher(rdi.nest_index));
            std::hash<int>  seeding_index_hasher;
            hash_combine(seed, seeding_index_hasher(rdi.seeding_index));
            std::hash<devel::scheme::ScaffoldIndex> scaffold_index_hasher;
            hash_combine(seed, scaffold_index_hasher(rdi.scaffold_index));

            return seed;
        }
    };

}


#endif
