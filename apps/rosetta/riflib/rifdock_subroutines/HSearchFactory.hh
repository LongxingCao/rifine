// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols

#ifndef INCLUDED_riflib_rifdock_subroutines_HSearchFactory_hh
#define INCLUDED_riflib_rifdock_subroutines_HSearchFactory_hh


#include <scheme/types.hh>
#include <riflib/rifdock_typedefs.hh>
#include <riflib/rifdock_subroutines/hsearch_original.hh>
#include <riflib/rifdock_subroutines/hsearch_morph_dive_pop.hh>


namespace devel {
namespace scheme {


bool
call_hsearch_protocol(
    RifDockData & rdd,
    HSearchData & d,
    shared_ptr<std::vector< SearchPointWithRots > > & hsearch_results_p
    ) {


    if (rdd.opt.scaff_search_mode == "default" || rdd.opt.scaff_search_mode == "nineA_baseline") {

        return hsearch_original( rdd, d, hsearch_results_p );

    } else if (rdd.opt.scaff_search_mode == "morph_dive_pop") {

        return hsearch_morph_dive_pop( rdd, d, hsearch_results_p );

    } else {
        utility_exit_with_message("Error!!!! -rif_dock:scaff_search_mode " + rdd.opt.scaff_search_mode + " does not name an HSearch protocol!");
    }

}

}}

#endif