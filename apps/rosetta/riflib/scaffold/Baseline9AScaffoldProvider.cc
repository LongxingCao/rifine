// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols


#include <riflib/types.hh>
#include <riflib/scaffold/Baseline9AScaffoldProvider.hh>
#include <riflib/scaffold/util.hh>
#include <ObjexxFCL/format.hh>
#include <riflib/scaffold/nineA_util.hh>
#include <riflib/scaffold/NineAManager.hh>

#include <string>
#include <vector>
#include <boost/any.hpp>
#include <boost/format.hpp>

using ::scheme::scaffold::BOGUS_INDEX;
using ::scheme::scaffold::TreeIndex;
using ::scheme::scaffold::TreeLimits;


namespace devel {
namespace scheme {


Baseline9AScaffoldProvider::Baseline9AScaffoldProvider( 
    uint64_t iscaff,
    shared_ptr< RotamerIndex > rot_index_p_in, 
    RifDockOpt const & opt_in,
    MakeTwobodyOpts const & make2bopts_in,
    ::devel::scheme::RotamerRFTablesManager & rotrf_table_manager_in ) :

    rot_index_p( rot_index_p_in), 
    opt(opt_in),
    make2bopts(make2bopts_in),
    rotrf_table_manager(rotrf_table_manager_in) {

    nineA_manager_ = NineAManager::get_instance( rot_index_p, opt );


    std::vector<uint64_t> cdindex_then_clusts = parse_nineA_baseline_range(opt.nineA_baseline_range);
    
    uint64_t cdindex = cdindex_then_clusts[0];
    uint64_t this_clust = cdindex_then_clusts[ iscaff + 1];

    nmember_ = nineA_manager_->get_nineA_member( cdindex, this_clust );
}




ScaffoldDataCacheOP 
Baseline9AScaffoldProvider::get_data_cache_slow(::scheme::scaffold::TreeIndex i) {

    return get_scaffold(i)->cache_data_;
}


ParametricSceneConformationCOP 
Baseline9AScaffoldProvider::get_scaffold(::scheme::scaffold::TreeIndex i) {

    if ( ! nmember_.conformation ) {
        utility_exit_with_message("Conformation not intialized yet!!");
    }

    return nmember_.conformation;
}


::scheme::scaffold::TreeLimits 
Baseline9AScaffoldProvider::get_scaffold_index_limits() const {
    return TreeLimits(1, 1);
}


void 
Baseline9AScaffoldProvider::set_fa_mode( bool fa ) {
    ScaffoldDataCacheOP cache = get_data_cache_slow( scaffold_index_default_value( ScaffoldIndex()) );
    if ( cache->conformation_is_fa != fa ) {
        nmember_.conformation = make_conformation_from_data_cache(cache, fa);
    }
}


::scheme::scaffold::TreeIndex 
Baseline9AScaffoldProvider::get_representative_scaffold_index() {
    scaffold_index_default_value( ScaffoldIndex());
}



void 
Baseline9AScaffoldProvider::setup_twobody_tables( ::scheme::scaffold::TreeIndex i ) {
    get_data_cache_slow( i )->setup_twobody_tables( rot_index_p, opt, make2bopts, rotrf_table_manager);
}








}}

