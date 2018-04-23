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


#endif
