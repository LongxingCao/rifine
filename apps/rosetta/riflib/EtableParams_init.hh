#ifndef INCLUDED_rosetta_objective_EtableParams_init_hh
#define INCLUDED_rosetta_objective_EtableParams_init_hh

#include "scheme/rosetta/score/EtableParams.hh"
#include <vector>

namespace devel { namespace scheme {

struct EtableParamsInit {
///@brief horrible function to fill horrible rosetta datastructure of LJ/LK params
static void init_EtableParams(
	std::vector< ::scheme::rosetta::score::EtableParamsOnePair<float> > & analytic_parameters
);

};

}}


#endif
