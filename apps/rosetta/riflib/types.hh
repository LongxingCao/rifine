// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.



#ifndef INCLUDED_riflib_types_hh
#define INCLUDED_riflib_types_hh

#include <scheme/types.hh>
#include <stdint.h>
#include <Eigen/Geometry>
#include <scheme/objective/voxel/VoxelArray.hh>

namespace devel {
namespace scheme {

	typedef ::Eigen::Transform<float,3,Eigen::AffineCompact> EigenXform;
	typedef ::scheme::objective::voxel::VoxelArray<3,float> VoxelArray;
	typedef VoxelArray* VoxelArrayPtr;
	using ::scheme::shared_ptr;
	using ::scheme::make_shared;
	using ::scheme::enable_shared_from_this;


}
}

#endif
