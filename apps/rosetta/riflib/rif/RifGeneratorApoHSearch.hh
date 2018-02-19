// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.



#ifndef INCLUDED_riflib_rif_RifGeneratorApoHSearch_hh
#define INCLUDED_riflib_rif_RifGeneratorApoHSearch_hh


#include <riflib/rif/RifGenerator.hh>

namespace devel {
namespace scheme {
namespace rif {

struct RifGeneratorApoHSearchOpts {
	float score_cut_adjust = 1.0;
	float abs_score_cut = 0.0;
	bool downweight_hydrophobics = false;
	float beam_size_M = 10000.0;
	float dump_fraction = 0.0;
};

struct RifGeneratorApoHSearch : public RifGenerator {

	RifGeneratorApoHSearchOpts opts;
	utility::vector1<std::string> apores;
	std::vector< std::vector< VoxelArray* > > const & bounding_by_atype;
	std::vector<float> RESLS;

	RifGeneratorApoHSearch(
		  utility::vector1<std::string> apores
		, std::vector< std::vector< VoxelArray* > > const & bounding_by_atype
		, std::vector< float > RESLS
		, RifGeneratorApoHSearchOpts opts
	)
	: opts( opts )
	, bounding_by_atype( bounding_by_atype )
	, RESLS( RESLS )
	, apores(apores)
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
