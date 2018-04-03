// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.



#ifndef INCLUDED_riflib_rif_requirements_util_hh
#define INCLUDED_riflib_rif_requirements_util__hh

#include <string>
#include <vector>


namespace devel {
namespace scheme {
namespace rif {

struct HBondDefinition {
    std::string atom_name = "";
    int			    res_num = -1;
    std::vector< std::string > allowed_rot_names;
};

struct BidentateDefinition {
    std::string atom1_name = "";
    int res1_num = -1;
    std::string atom2_name = "";
    int res2_num = -1;
};

struct RequirementDefinition {
    int req_num;
    std::string require;
    std::vector< std::string > definition;
};

std::vector< HBondDefinition > get_hbond_definitions( std::string tuning_file );
std::vector< BidentateDefinition > get_bidentate_definitions( std::string tuning_file );
std::vector< RequirementDefinition > get_requirement_definitions( std::string tuning_file );

}
}
}

#endif