// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.



#ifndef INCLUDED_riflib_rif_RifGeneratorUserHotspots_hh
#define INCLUDED_riflib_rif_RifGeneratorUserHotspots_hh


#include <riflib/rif/RifGenerator.hh>
#include <scheme/chemical/RotamerIndex.hh>	


namespace devel {
namespace scheme {
namespace rif {

struct RifGeneratorUserHotspotsOpts {
	float hotspot_sample_cart_bound = 1.0;
	float hotspot_sample_angle_bound = 30.0;
    float hotspot_score_thresh = -0.5;
    int   hotspot_nsamples = 10000;
	float hbond_weight = 2.0;
	float upweight_multi_hbond = 0.0;
    int  dump_hotspot_samples = 0;
    bool single_file_hotspots_insertion = false;
	Eigen::Vector3f target_center;
	std::vector<std::string> hotspot_files;
};

struct RifGeneratorUserHotspots : public RifGenerator {

	RifGeneratorUserHotspotsOpts opts;

	RifGeneratorUserHotspots( RifGeneratorUserHotspotsOpts const & _opts) : opts(_opts) {}

	void generate_rif(
		RifAccumulatorP accumulator,
		RifGenParamsP params
	) override;

	void 
	modify_rotamer_spec(
		::scheme::chemical::RotamerIndexSpec& rot_spec
	) override;

};


}
}
}

#endif
