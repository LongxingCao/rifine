// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols

#ifndef INCLUDED_riflib_scaffold_MorphingScaffoldProvider_hh
#define INCLUDED_riflib_scaffold_MorphingScaffoldProvider_hh

#include <riflib/types.hh>
#include <riflib/rifdock_typedefs.hh>
#include <riflib/scaffold/morph_util.hh>
#include <scheme/scaffold/ScaffoldProviderBase.hh>
#include <riflib/scaffold/ScaffoldDataCache.hh>

// temporary, delete this
#include <riflib/scaffold/NineAManager.hh>

#include <string>
#include <vector>
#include <boost/any.hpp>

#include <scheme/kinematics/Scene.hh>

#include <core/pose/Pose.hh>



namespace devel {
namespace scheme {




struct MorphingScaffoldProvider :
    public ::scheme::scaffold::TreeScaffoldProvider<ParametricSceneConformation> {

    // MorphingScaffoldProvider(); 

    MorphingScaffoldProvider( 
        uint64_t iscaff,
        shared_ptr< RotamerIndex > rot_index_p_in, 
        RifDockOpt const & opt_in,
        MakeTwobodyOpts const & make2bopts_in,
        ::devel::scheme::RotamerRFTablesManager & rotrf_table_manager_in );


    ParametricSceneConformationCOP get_scaffold(::scheme::scaffold::TreeIndex i) override;

    ::scheme::scaffold::TreeLimits get_scaffold_index_limits() const override;

    ScaffoldDataCacheOP get_data_cache_slow(::scheme::scaffold::TreeIndex i) override;

    void set_fa_mode( bool fa ) override;

    ::scheme::scaffold::TreeIndex get_representative_scaffold_index() override;

    void setup_twobody_tables( ::scheme::scaffold::TreeIndex i ) override;

    void modify_pose_for_output( ::scheme::scaffold::TreeIndex i, core::pose::Pose & pose ) override;

    void test_make_children(::scheme::scaffold::TreeIndex ti);

private:


    MorphMember & 
    get_morph_member(::scheme::scaffold::TreeIndex i);

    ::scheme::scaffold::TreeIndex
    add_morph_member( MorphMember mmember );

    // By definition, all members at depth 0 are provided by
    // the user
    std::vector< std::vector< MorphMember > > map_;


    shared_ptr< RotamerIndex > rot_index_p;
    RifDockOpt const & opt;
    MakeTwobodyOpts const & make2bopts;
    ::devel::scheme::RotamerRFTablesManager & rotrf_table_manager ;

};


}}

#endif
