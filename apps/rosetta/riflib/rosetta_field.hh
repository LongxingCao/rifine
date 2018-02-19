
// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols
#ifndef INCLUDED_riflib_rosetta_field_hh
#define INCLUDED_riflib_rosetta_field_hh

#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>
#include <scheme/actor/Atom.hh>
#include "scheme/objective/voxel/VoxelArray.hh"
#include <vector>
#include <core/types.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>

#include <Eigen/Dense>


namespace devel {
namespace scheme {

struct RosettaFieldOptions {
	bool block_hbond_sites = false;
	float field_resl = 0.25;
	float max_bounding_ratio = 4.0;
	std::string data_dir = "";
	int oversample = 2;
	bool repulsive_only_boundary = false; // atoms not in target_res are repulsive-only
	bool fail_if_no_cached_data = false;
	bool cache_mismatch_tolerance = 0.01;
	bool generate_only = false;
    int one_atype_only = 0;
	std::vector<core::id::AtomID> repulsive_atoms;
};

void
parse_atomids(
	core::pose::Pose const & refpose,
	utility::vector1<int> const & atomreslist,
	std::vector<core::id::AtomID> & atomids_out,
	std::string printme = ""
);

std::string
get_rosetta_fields(
	std::string const & target_fname,
	core::pose::Pose const & target,
	utility::vector1<core::Size> const & target_res,
	RosettaFieldOptions const & opts,
	std::vector<
		// shared_ptr<
			::scheme::objective::voxel::VoxelArray<3,float> *
		// >
	> & field_by_atype,
	bool verbose = true
);

std::string
get_rosetta_fields_specified_cache_prefix(
	std::string const & cache_prefix,
	std::string const & target_fname,
	core::pose::Pose const & target,
	utility::vector1<core::Size> const & target_res,
	RosettaFieldOptions const & opts,
	std::vector<
		// shared_ptr<
			::scheme::objective::voxel::VoxelArray<3,float> *
		// >
	> & field_by_atype,
	bool verbose = true
);


std::string
get_rosetta_bounding_fields(
	std::vector<float> const & RESLS,
	std::string const & target_fname,
	core::pose::Pose const & target,
	utility::vector1<core::Size> const & target_res,
	RosettaFieldOptions const & opts,
	std::vector<              ::scheme::objective::voxel::VoxelArray<3,float> *   > & field_by_atype,
	std::vector< std::vector< ::scheme::objective::voxel::VoxelArray<3,float> *	> > & bounding_by_atype,
	bool verbose = true
);

void
get_rosetta_bounding_fields_from_fba(
	std::vector<float> const & RESLS,
	std::string const & target_fname,
	core::pose::Pose const & target,
	utility::vector1<core::Size> const & target_res,
	RosettaFieldOptions const & opts,
	std::vector<              ::scheme::objective::voxel::VoxelArray<3,float> *   > & field_by_atype,
	std::vector< std::vector< ::scheme::objective::voxel::VoxelArray<3,float> *	> > & bounding_by_atype,
	bool verbose,
	std::string cache_prefix
);

void
get_scheme_atoms(
	core::pose::Pose const & target,
	utility::vector1<core::Size> const & target_res,
	std::vector< ::scheme::actor::Atom< Eigen::Vector3f > > & out,
	bool bbonly = false
);

void
get_scheme_atoms_cbonly(
	core::pose::Pose const & target,
	utility::vector1<core::Size> const & target_res,
	std::vector< ::scheme::actor::Atom< Eigen::Vector3f > > & out
);


void
get_scheme_atoms(
	core::pose::Pose const & target,
	std::vector< ::scheme::actor::Atom< Eigen::Vector3f > > & out,
	bool bbonly = false
);


}}

#endif
