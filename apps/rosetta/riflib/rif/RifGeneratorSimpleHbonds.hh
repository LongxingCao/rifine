// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.



#ifndef INCLUDED_riflib_rif_RifGeneratorSimpleHbonds_hh
#define INCLUDED_riflib_rif_RifGeneratorSimpleHbonds_hh


#include <riflib/rif/RifGenerator.hh>

namespace devel {
namespace scheme {
namespace rif {

struct RifGeneratorSimpleHbondsOpts {
	float tip_tol_deg = 30.0;
	float rot_samp_resl = 6.0;
	float rot_samp_range = 360.0;
	float hbond_cart_sample_hack_range = 0.25;
	float hbond_cart_sample_hack_resl = 0.25;
	float score_threshold = 0.0;
	float dump_fraction = 0.0;
	float hbond_weight = 2.0;
	float upweight_multi_hbond = 0;
	bool debug = false;
};

struct RifGeneratorSimpleHbonds : public RifGenerator {

	utility::vector1<std::string> donresn_user;
	utility::vector1<std::string> accresn_user;
	RifGeneratorSimpleHbondsOpts opts;


	RifGeneratorSimpleHbonds(
		  utility::vector1<std::string> donresn_user
		, utility::vector1<std::string> accresn_user
		, RifGeneratorSimpleHbondsOpts opts
	)
	: donresn_user( donresn_user )
	, accresn_user( accresn_user )

	, opts( opts )
	{}

	void generate_rif(
		RifAccumulatorP accumulator,
		RifGenParamsP params
	) override;

};

}
}
}

#endif
