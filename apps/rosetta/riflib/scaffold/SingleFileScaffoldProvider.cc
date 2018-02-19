// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols

#include <riflib/scaffold/SingleFileScaffoldProvider.hh>
#include <riflib/scaffold/util.hh>

#include <riflib/types.hh>
#include <scheme/numeric/rand_xform.hh>
#include <core/import_pose/import_pose.hh>
#include <utility/file/file_sys_util.hh>
#include <riflib/HSearchConstraints.hh>

#include <string>
#include <vector>
#include <boost/any.hpp>



namespace devel {
namespace scheme {


// SingleFileScaffoldProvider::SingleFileScaffoldProvider() {}

SingleFileScaffoldProvider::SingleFileScaffoldProvider( 
    uint64_t iscaff,
    shared_ptr< RotamerIndex > rot_index_p_in, 
    RifDockOpt const & opt_in,
    MakeTwobodyOpts const & make2bopts_in,
    ::devel::scheme::RotamerRFTablesManager & rotrf_table_manager_in ) :

    rot_index_p( rot_index_p_in), 
    opt(opt_in),
    make2bopts(make2bopts_in),
    rotrf_table_manager(rotrf_table_manager_in) {


    std::string scafftag;
    core::pose::Pose scaffold;
    utility::vector1<core::Size> scaffold_res;
    EigenXform scaffold_perturb;
    std::vector<CstBaseOP> csts;
    MorphRules morph_rules;

    get_info_for_iscaff( iscaff, opt, scafftag, scaffold, scaffold_res, scaffold_perturb, csts, morph_rules);

    ScaffoldDataCacheOP temp_data_cache = make_shared<ScaffoldDataCache>(
        scaffold,
        scaffold_res,
        scafftag,
        scaffold_perturb,
        rot_index_p,
        opt,
        csts);

    conformation_ = make_conformation_from_data_cache(temp_data_cache, false);

}


ScaffoldDataCacheOP 
SingleFileScaffoldProvider::get_data_cache_slow(::scheme::scaffold::TreeIndex i) {

    return get_scaffold(i)->cache_data_;

}


ParametricSceneConformationCOP 
SingleFileScaffoldProvider::get_scaffold(::scheme::scaffold::TreeIndex i) {
    if ( ! conformation_ ) {
        utility_exit_with_message("Conformation not intialized yet!!");
    }
    return conformation_;
}


::scheme::scaffold::TreeLimits 
SingleFileScaffoldProvider::get_scaffold_index_limits() const {
    return ::scheme::scaffold::TreeLimits(1, 1);
}


void 
SingleFileScaffoldProvider::set_fa_mode( bool fa ) {
    ScaffoldDataCacheOP cache = get_data_cache_slow( scaffold_index_default_value( ScaffoldIndex()) );
    if ( cache->conformation_is_fa != fa ) {
        conformation_ = make_conformation_from_data_cache(cache, fa);
    }
}


::scheme::scaffold::TreeIndex 
SingleFileScaffoldProvider::get_representative_scaffold_index() {
    return scaffold_index_default_value( ScaffoldIndex());
}


void 
SingleFileScaffoldProvider::setup_twobody_tables( ::scheme::scaffold::TreeIndex i ) {
    get_data_cache_slow( i )->setup_twobody_tables( rot_index_p, opt, make2bopts, rotrf_table_manager);
}



}}

