// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols

#ifndef INCLUDED_riflib_scaffold_util_hh
#define INCLUDED_riflib_scaffold_util_hh

#include <riflib/types.hh>
#include <riflib/rifdock_typedefs.hh>
#include <rif_dock_test.hh>
#include <riflib/HSearchConstraints.hh>
#include <riflib/scaffold/morph_util.hh>

#include <protocols/minimization_packing/TaskAwareMinMover.hh>
#include <protocols/simple_moves/MutateResidue.hh>
#include <core/scoring/ScoreFunction.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>

#include <string>
#include <vector>


#ifdef USEHDF5
#include <protocols/indexed_structure_store/movers/DirectSegmentLookupMover.hh>
#else
#endif



namespace devel {
namespace scheme {

struct ScaffoldDataCache;

void
get_info_for_iscaff(
    uint64_t iscaff,
    RifDockOpt const & opt, 
    std::string & scafftag,
    core::pose::Pose & scaffold,
    utility::vector1<core::Size> & scaffold_res,
    EigenXform & scaffold_perturb,
    std::vector<CstBaseOP> & csts,
    MorphRules & morph_rules
    );

// historically, non_fa was used during HSearch and fa was used during hack pack
ParametricSceneConformationCOP
make_conformation_from_data_cache(shared_ptr<ScaffoldDataCache> cache, bool fa = false) ;


struct ApplyDSLMScratch {
    std::vector<core::scoring::ScoreFunctionOP> scorefxn_pt;
    std::vector<protocols::simple_moves::MutateResidueOP> to_ala_pt;
    std::vector<protocols::minimization_packing::TaskAwareMinMoverOP> hardmin_bb_pt;
    std::vector<core::scoring::ScoreFunctionOP> rep_scorefxn_pt;
};

#ifdef USEHDF5
std::vector<core::pose::PoseOP>
apply_direct_segment_lookup_mover( 
    protocols::indexed_structure_store::movers::DirectSegmentLookupMover & dsl_mover,
    core::pose::Pose const & pose,
    uint64_t low_cut_site,
    uint64_t high_cut_site,
    uint64_t minimum_loop_length,
    uint64_t max_structures,
    uint64_t max_rmsd,
    ApplyDSLMScratch & scratch   );
#endif

void
get_default_scaffold_res( core::pose::Pose const & pose,
    utility::vector1<core::Size> & scaffold_res);

void
add_pdbinfo_if_missing( core::pose::Pose & pose );


bool
internal_comparative_clash_check( core::scoring::Energies const & original_energies,
    core::pose::Pose const & to_check,
    core::scoring::ScoreFunctionOP scorefxn,
    core::Real max_fa_rep_delta,
    core::Size low_position,
    core::Size high_position );










}}

#endif
