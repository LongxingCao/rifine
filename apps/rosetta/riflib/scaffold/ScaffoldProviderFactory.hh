// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols

#ifndef INCLUDED_riflib_scaffold_ScaffoldProviderFactory_hh
#define INCLUDED_riflib_scaffold_ScaffoldProviderFactory_hh


#include <scheme/types.hh>
#include <riflib/scaffold/Baseline9AScaffoldProvider.hh>
#include <riflib/scaffold/MorphingScaffoldProvider.hh>
#include <riflib/scaffold/SingleFileScaffoldProvider.hh>
#include <riflib/rifdock_typedefs.hh>


namespace devel {
namespace scheme {


ScaffoldProviderOP
get_scaffold_provider( 
        uint64_t iscaff,
        shared_ptr< RotamerIndex > rot_index_p, 
        RifDockOpt const & opt,
        MakeTwobodyOpts const & make2bopts,
        ::devel::scheme::RotamerRFTablesManager & rotrf_table_manager ) {


    if (opt.scaff_search_mode == "default") {

        return make_shared<SingleFileScaffoldProvider>(
                iscaff,
                rot_index_p,
                opt,
                make2bopts,
                rotrf_table_manager);

    } else if (opt.scaff_search_mode == "morph_dive_pop") {

        return make_shared<MorphingScaffoldProvider>(
                iscaff,
                rot_index_p,
                opt,
                make2bopts,
                rotrf_table_manager);

    } else if (opt.scaff_search_mode == "nineA_baseline") {

        return make_shared<Baseline9AScaffoldProvider>(
                iscaff,
                rot_index_p,
                opt,
                make2bopts,
                rotrf_table_manager);
    } else {
        utility_exit_with_message("Error!!!! -rif_dock:scaff_search_mode " + opt.scaff_search_mode + " does not name a ScaffoldProvider!");
    }

}

}}

#endif