// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols

// keep this file super lightweight. It gets included in rif_dock_test.hh

#ifndef INCLUDED_riflib_scaffold_nineA_util_hh
#define INCLUDED_riflib_scaffold_nineA_util_hh

#include <riflib/types.hh>



#include <string>
#include <vector>



namespace devel {
namespace scheme {




inline
std::vector<uint64_t>
parse_nineA_baseline_range(std::string const & nineA_baseline_range) {
    std::vector<uint64_t> cdindex_then_clusts;
    try {
        utility::vector1<std::string> cdindex_range = utility::string_split( nineA_baseline_range, ':' );
        cdindex_then_clusts.push_back(stoi(cdindex_range[1]));
        if ( cdindex_range.size() == 2 ) {
            utility::vector1<std::string> clust_elements = utility::string_split( cdindex_range[2], ',' );
            for ( std::string const & element : clust_elements ) {
                utility::vector1<std::string> maybe_range = utility::string_split( element, '-' );

                if ( maybe_range.size() == 1 ) {
                    cdindex_then_clusts.push_back(stoi(maybe_range[1]));
                } else if ( maybe_range.size() == 2 ) {
                    uint64_t lo = stoi(maybe_range[1]);
                    uint64_t hi = stoi(maybe_range[2]);
                    runtime_assert( hi >= lo );
                    for ( uint64_t i = lo; i < hi; i++ ) {
                        cdindex_then_clusts.push_back( i );
                    }
                } else {
                    utility_exit_with_message("nineA_baseline_range parse fail on element: " + element);
                }

            }

        }
    } catch ( std::exception const & ex ) {
    }
    if (cdindex_then_clusts.size() <= 1) {
        utility_exit_with_message("nineA_baseline_range parse fail!" );
    }
    return cdindex_then_clusts;
}




}}

#endif
