// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:

// INC



#include <basic/options/option_macros.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <riflib/scaffold/nineA_util.hh>
#include <vector>


#ifndef INCLUDED_rif_dock_test_hh_1
#define INCLUDED_rif_dock_test_hh_1



OPT_1GRP_KEY(     StringVector , rif_dock, scaffolds )
OPT_1GRP_KEY(  StringVector, rif_dock, scaffold_res )
OPT_1GRP_KEY(  StringVector, rif_dock, scaffold_res_fixed )
OPT_1GRP_KEY(  Boolean     , rif_dock, scaffold_res_use_best_guess )
OPT_1GRP_KEY(  Boolean     , rif_dock, scaffold_to_ala )
OPT_1GRP_KEY(  Boolean     , rif_dock, scaffold_to_ala_selonly )
OPT_1GRP_KEY(  Boolean     , rif_dock, replace_orig_scaffold_res )
OPT_1GRP_KEY(  Boolean     , rif_dock, replace_all_with_ala_1bre )
OPT_1GRP_KEY(  Boolean     , rif_dock, random_perturb_scaffold )
OPT_1GRP_KEY(  Boolean     , rif_dock, dont_center_scaffold )

OPT_1GRP_KEY(  StringVector, rif_dock, target_bounding_xmaps )
OPT_1GRP_KEY(  String      , rif_dock, target_pdb )
OPT_1GRP_KEY(  String      , rif_dock, target_res )
OPT_1GRP_KEY(  String      , rif_dock, target_rif )
OPT_1GRP_KEY(  Real        , rif_dock, target_rf_resl )
OPT_1GRP_KEY(  Integer     , rif_dock, target_rf_oversample )
OPT_1GRP_KEY(  String      , rif_dock, target_rf_cache )

OPT_1GRP_KEY(  StringVector, rif_dock, data_cache_dir )

OPT_1GRP_KEY(  Real        , rif_dock, beam_size_M )
OPT_1GRP_KEY(  Real        , rif_dock, search_diameter )
OPT_1GRP_KEY(  Real        , rif_dock, hsearch_scale_factor )

OPT_1GRP_KEY(  Real        , rif_dock, max_rf_bounding_ratio )
OPT_1GRP_KEY(  Boolean     , rif_dock, make_bounding_plot_data )
OPT_1GRP_KEY(  Boolean     , rif_dock, align_output_to_scaffold )
OPT_1GRP_KEY(  Boolean     , rif_dock, output_scaffold_only )
OPT_1GRP_KEY(  Boolean     , rif_dock, output_full_scaffold_only )
OPT_1GRP_KEY(  Boolean     , rif_dock, output_full_scaffold )
OPT_1GRP_KEY(  Integer     , rif_dock, n_pdb_out )

OPT_1GRP_KEY(  Real        , rif_dock, rf_resl )
OPT_1GRP_KEY(  Integer     , rif_dock, rf_oversample )
OPT_1GRP_KEY(  Boolean     , rif_dock, downscale_atr_by_hierarchy )
OPT_1GRP_KEY(  Real        , rif_dock, favorable_1body_multiplier )
OPT_1GRP_KEY(  Real        , rif_dock, favorable_1body_multiplier_cutoff )
OPT_1GRP_KEY(  Real        , rif_dock, favorable_2body_multiplier )

OPT_1GRP_KEY(  Integer     , rif_dock, rotrf_oversample )
OPT_1GRP_KEY(  Real        , rif_dock, rotrf_resl )
OPT_1GRP_KEY(  Real        , rif_dock, rotrf_spread )
OPT_1GRP_KEY(  Real        , rif_dock, rotrf_scale_atr )
OPT_1GRP_KEY(  String      , rif_dock, rotrf_cache_dir )

OPT_1GRP_KEY(  Boolean     , rif_dock, hack_pack )
OPT_1GRP_KEY(  Real        , rif_dock, hack_pack_frac )
OPT_1GRP_KEY(  Real        , rif_dock, pack_iter_mult )
OPT_1GRP_KEY(  Integer     , rif_dock, pack_n_iters )
OPT_1GRP_KEY(  Real        , rif_dock, hbond_weight )
OPT_1GRP_KEY(  Real        , rif_dock, upweight_multi_hbond )
OPT_1GRP_KEY(  Real        , rif_dock, global_score_cut )

OPT_1GRP_KEY(  Integer     , rif_dock, n_result_limit )
OPT_1GRP_KEY(  Real        , rif_dock, redundancy_filter_mag )

OPT_1GRP_KEY(  Real        , rif_dock, force_output_if_close_to_input )
OPT_1GRP_KEY(  Integer     , rif_dock, force_output_if_close_to_input_num )

OPT_1GRP_KEY(  Real        , rif_dock, upweight_iface )

OPT_1GRP_KEY(  Boolean     , rif_dock, use_scaffold_bounding_grids )

OPT_1GRP_KEY(  Boolean     , rif_dock, restrict_to_native_scaffold_res )
OPT_1GRP_KEY(  Real        , rif_dock, bonus_to_native_scaffold_res )
OPT_1GRP_KEY(  Boolean     , rif_dock, add_native_scaffold_rots_when_packing )

OPT_1GRP_KEY(  Boolean     , rif_dock, dump_all_rif_rots )
OPT_1GRP_KEY(  Boolean     , rif_dock, dump_all_rif_rots_into_output )
OPT_1GRP_KEY(  Boolean     , rif_dock, rif_rots_as_chains )

OPT_1GRP_KEY(  String      , rif_dock, dump_rifgen_near_pdb )
OPT_1GRP_KEY(  Real        , rif_dock, dump_rifgen_near_pdb_dist )
OPT_1GRP_KEY(  Real        , rif_dock, dump_rifgen_near_pdb_frac )

OPT_1GRP_KEY(  String     , rif_dock, dokfile )
OPT_1GRP_KEY(  String     , rif_dock, outdir )
OPT_1GRP_KEY(  String     , rif_dock, output_tag )

OPT_1GRP_KEY(  Boolean    , rif_dock, dont_use_scaffold_loops )

OPT_1GRP_KEY(  Boolean    , rif_dock, dump_resfile )
OPT_1GRP_KEY(  Boolean    , rif_dock, pdb_info_pikaa )

OPT_1GRP_KEY(  Boolean    , rif_dock, cache_scaffold_data )

OPT_1GRP_KEY(  Real        , rif_dock, tether_to_input_position )

OPT_1GRP_KEY(  Boolean     , rif_dock, lowres_sterics_cbonly )

OPT_1GRP_KEY(  Integer     , rif_dock, require_satisfaction )
OPT_1GRP_KEY(  Integer     , rif_dock, require_n_rifres )

OPT_1GRP_KEY(  Real        , rif_dock, rosetta_score_fraction )
OPT_1GRP_KEY(  Real        , rif_dock, rosetta_score_then_min_below_thresh )
OPT_1GRP_KEY(  Integer     , rif_dock, rosetta_score_at_least )
OPT_1GRP_KEY(  Integer     , rif_dock, rosetta_score_at_most )
OPT_1GRP_KEY(  Integer     , rif_dock, rosetta_min_at_least )
OPT_1GRP_KEY(  Integer     , rif_dock, rosetta_min_at_most )
OPT_1GRP_KEY(  Real        , rif_dock, rosetta_min_fraction )
OPT_1GRP_KEY(  Boolean     , rif_dock, rosetta_min_fix_target )
OPT_1GRP_KEY(  Boolean     , rif_dock, rosetta_min_targetbb )
OPT_1GRP_KEY(  Boolean     , rif_dock, rosetta_min_scaffoldbb )
OPT_1GRP_KEY(  Boolean     , rif_dock, rosetta_min_allbb )
OPT_1GRP_KEY(  Real        , rif_dock, rosetta_score_cut )
OPT_1GRP_KEY(  Boolean     , rif_dock, rosetta_hard_min )
OPT_1GRP_KEY(  Boolean     , rif_dock, rosetta_score_total )
OPT_1GRP_KEY(  Boolean     , rif_dock, rosetta_score_ddg_only )
OPT_1GRP_KEY(  Real        , rif_dock, rosetta_score_rifres_rifres_weight )
OPT_1GRP_KEY(  Real        , rif_dock, rosetta_score_rifres_scaffold_weight )
OPT_1GRP_KEY(  String      , rif_dock, rosetta_soft_score )
OPT_1GRP_KEY(  String      , rif_dock, rosetta_hard_score )

OPT_1GRP_KEY(  Boolean     , rif_dock, extra_rotamers )
OPT_1GRP_KEY(  Boolean     , rif_dock, extra_rif_rotamers )
OPT_1GRP_KEY(  Integer     , rif_dock, always_available_rotamers_level )
OPT_1GRP_KEY(  Boolean     , rif_dock, packing_use_rif_rotamers )

OPT_1GRP_KEY(  Integer     , rif_dock, nfold_symmetry )
OPT_1GRP_KEY(  RealVector  , rif_dock, symmetry_axis )

OPT_1GRP_KEY(  Real        , rif_dock, user_rotamer_bonus_constant )
OPT_1GRP_KEY(  Real        , rif_dock, user_rotamer_bonus_per_chi )


OPT_1GRP_KEY(  Real        , rif_dock, resl0 )

OPT_1GRP_KEY(  Integer     , rif_dock, dump_x_frames_per_resl )
OPT_1GRP_KEY(  Boolean     , rif_dock, dump_only_best_frames )
OPT_1GRP_KEY(  Integer     , rif_dock, dump_only_best_stride )
OPT_1GRP_KEY(  String      , rif_dock, dump_prefix )

OPT_1GRP_KEY(  String      , rif_dock, scaff_search_mode )
OPT_1GRP_KEY(  String      , rif_dock, nineA_cluster_path )
OPT_1GRP_KEY(  String      , rif_dock, nineA_baseline_range )

OPT_1GRP_KEY(  Integer     , rif_dock, low_cut_site )
OPT_1GRP_KEY(  Integer     , rif_dock, high_cut_site )
OPT_1GRP_KEY(  Integer     , rif_dock, max_insertion )
OPT_1GRP_KEY(  Integer     , rif_dock, max_deletion )
OPT_1GRP_KEY(  Real        , rif_dock, fragment_cluster_tolerance )
OPT_1GRP_KEY(  Real        , rif_dock, fragment_max_rmsd )
OPT_1GRP_KEY(  Integer     , rif_dock, max_fragments )
OPT_1GRP_KEY(  StringVector, rif_dock, morph_rules_files )

OPT_1GRP_KEY(  Boolean     , rif_dock, include_parent )
OPT_1GRP_KEY(  Boolean     , rif_dock, use_parent_body_energies )

OPT_1GRP_KEY(  Integer     , rif_dock, dive_resl )
OPT_1GRP_KEY(  Integer     , rif_dock, pop_resl )
OPT_1GRP_KEY(  String      , rif_dock, match_this_pdb )
OPT_1GRP_KEY(  Real        , rif_dock, match_this_rmsd )
OPT_1GRP_KEY(  Real        , rif_dock, max_beam_multiplier )

OPT_1GRP_KEY(  String      , rif_dock, rot_spec_fname )
// constrain file
OPT_1GRP_KEY(  StringVector, rif_dock, cst_files )
OPT_1GRP_KEY(  StringVector, rif_dock, seeding_pos )
OPT_1GRP_KEY(  Boolean     , rif_dock, seeding_by_patchdock )
OPT_1GRP_KEY(  String      , rif_dock, xform_pos )
OPT_1GRP_KEY(  Boolean     , rif_dock, test_longxing )
OPT_1GRP_KEY(  Integer     , rif_dock, rosetta_score_each_seeding_at_least )
OPT_1GRP_KEY(  Real        , rif_dock, cluster_score_cut )
OPT_1GRP_KEY(  Real        , rif_dock, keep_top_clusters_frac )
OPT_1GRP_KEY(  Boolean     , rif_dock, only_dump_scaffold )
OPT_1GRP_KEY(  IntegerVector, rif_dock, requirements )

OPT_1GRP_KEY(  StringVector, rif_dock, dump_residue_names   )
OPT_1GRP_KEY(  Real        , rif_dock, dump_fraction        )
OPT_1GRP_KEY(  Real        , rif_dock, dump_score_cutoff        )
OPT_1GRP_KEY(  Integer     , rif_dock, dump_requirement_num      )
OPT_1GRP_KEY(  String      , rif_dock, dump_fname           )







void register_options() {
    using namespace basic::options;
    using namespace basic::options::OptionKeys;
    
    NEW_OPT(  rif_dock::scaffolds, "" , utility::vector1<std::string>() );
    NEW_OPT(  rif_dock::scaffold_res, "" , utility::vector1<std::string>() );
    NEW_OPT(  rif_dock::scaffold_res_fixed, "" , utility::vector1<std::string>() );
    NEW_OPT(  rif_dock::scaffold_res_use_best_guess, "" , false );
    NEW_OPT(  rif_dock::scaffold_to_ala, "" , false );
    NEW_OPT(  rif_dock::scaffold_to_ala_selonly, "" , true );
    NEW_OPT(  rif_dock::replace_orig_scaffold_res, "", true );
    NEW_OPT(  rif_dock::replace_all_with_ala_1bre, "" , false );
    NEW_OPT(  rif_dock::random_perturb_scaffold, "" , false );
    NEW_OPT(  rif_dock::dont_center_scaffold, "don't use this", false );
    
    NEW_OPT(  rif_dock::target_bounding_xmaps, "" , utility::vector1<std::string>() );
    NEW_OPT(  rif_dock::target_pdb, "" , "" );
    NEW_OPT(  rif_dock::target_res, "" , "" );
    NEW_OPT(  rif_dock::target_rif, "" , "" );
    NEW_OPT(  rif_dock::target_rf_resl, ""       , 0.25 );
    NEW_OPT(  rif_dock::target_rf_oversample, "" , 2 );
    NEW_OPT(  rif_dock::downscale_atr_by_hierarchy, "" , true );
    NEW_OPT(  rif_dock::favorable_1body_multiplier, "Anything with a one-body energy less than xx gets multiplied by this", 1 );
    NEW_OPT(  rif_dock::favorable_1body_multiplier_cutoff, "Anything with a one-body energy less than cutoff gets multiplied by this", 4 );
    NEW_OPT(  rif_dock::favorable_2body_multiplier, "Anything with a two-body energy less than 0 gets multiplied by this", 1 );
    
    NEW_OPT(  rif_dock::target_rf_cache, "" , "NO_CACHE_SPECIFIED_ON_COMMAND_LINE" );
    
    NEW_OPT(  rif_dock::data_cache_dir, "" , utility::vector1<std::string>(1,"./") );
    NEW_OPT(  rif_dock::beam_size_M, "" , 10.000000 );
    NEW_OPT(  rif_dock::max_rf_bounding_ratio, "" , 4 );
    NEW_OPT(  rif_dock::make_bounding_plot_data, "" , false );
    NEW_OPT(  rif_dock::align_output_to_scaffold, "" , false );
    NEW_OPT(  rif_dock::output_scaffold_only, "" , false );
    NEW_OPT(  rif_dock::output_full_scaffold_only, "" , false );
    NEW_OPT(  rif_dock::output_full_scaffold, "", false );
    NEW_OPT(  rif_dock::n_pdb_out, "" , 10 );
    
    NEW_OPT(  rif_dock::rf_resl, ""       , 0.25 );
    NEW_OPT(  rif_dock::rf_oversample, "" , 2 );
    
    NEW_OPT(  rif_dock::rotrf_oversample, "" , 2 );
    NEW_OPT(  rif_dock::rotrf_resl, "" , 0.3 );
    NEW_OPT(  rif_dock::rotrf_spread, "" , 0.0 );
    NEW_OPT(  rif_dock::rotrf_scale_atr, "" , 1.0 );
    NEW_OPT(  rif_dock::rotrf_cache_dir, "" , "./" );
    
    NEW_OPT(  rif_dock::hack_pack, "" , true );
    NEW_OPT(  rif_dock::hack_pack_frac, "" , 0.2 );
    NEW_OPT(  rif_dock::pack_iter_mult, "" , 2.0 );
    NEW_OPT(  rif_dock::pack_n_iters, "" , 1 );
    NEW_OPT(  rif_dock::hbond_weight, "" , 2.0 );
    NEW_OPT(  rif_dock::upweight_multi_hbond, "" , 0.0 );
    NEW_OPT(  rif_dock::global_score_cut, "" , 0.0 );
    
    NEW_OPT(  rif_dock::n_result_limit, "" , 2000000000 );
    
    NEW_OPT(  rif_dock::redundancy_filter_mag, "" , 1.0 );
    
    NEW_OPT(  rif_dock::force_output_if_close_to_input, "" , 1.0 );
    NEW_OPT(  rif_dock::force_output_if_close_to_input_num, "" , 0 );
    
    NEW_OPT(  rif_dock::upweight_iface, "", 1.2 );
    
    NEW_OPT(  rif_dock::use_scaffold_bounding_grids, "", false );
    
    NEW_OPT(  rif_dock::search_diameter, "", 150.0 );
    NEW_OPT(  rif_dock::hsearch_scale_factor, "global scaling of rotation/translation search grid", 1.0 );
    
    NEW_OPT(  rif_dock::restrict_to_native_scaffold_res, "aka structure prediction CHEAT", false );
    NEW_OPT(  rif_dock::bonus_to_native_scaffold_res, "aka favor native CHEAT", -0.3 );
    NEW_OPT(  rif_dock::add_native_scaffold_rots_when_packing, "CHEAT", false );
    
    NEW_OPT(  rif_dock::dump_all_rif_rots, "", false );
    NEW_OPT(  rif_dock::dump_all_rif_rots_into_output, "dump all rif rots into output", false);
    NEW_OPT(  rif_dock::rif_rots_as_chains, "dump rif rots as chains instead of models, loses resnum if true", false );
    
    NEW_OPT(  rif_dock::dump_rifgen_near_pdb, "dump rifgen rotamers with same AA type near this single residue", "");
    NEW_OPT(  rif_dock::dump_rifgen_near_pdb_dist, "", 1 );
    NEW_OPT(  rif_dock::dump_rifgen_near_pdb_frac, "", 1 );
    
    NEW_OPT(  rif_dock::dokfile, "", "default.dok" );
    NEW_OPT(  rif_dock::outdir, "", "./" );
    NEW_OPT(  rif_dock::output_tag, "", "" );
    
    NEW_OPT(  rif_dock::dont_use_scaffold_loops, "", false );
    
    NEW_OPT(  rif_dock::dump_resfile, "", false );
    NEW_OPT(  rif_dock::pdb_info_pikaa, "", false );
    
    NEW_OPT(  rif_dock::cache_scaffold_data, "", false );
    
    NEW_OPT(  rif_dock::tether_to_input_position, "", -1.0 );
    
    NEW_OPT(  rif_dock::lowres_sterics_cbonly, "", true );
    
    NEW_OPT(  rif_dock::require_satisfaction, "", 0 );
    NEW_OPT(  rif_dock::require_n_rifres, "This doesn't work during HackPack", 0 );
    
    NEW_OPT(  rif_dock::rosetta_score_fraction  , "",  0.00 );
    NEW_OPT(  rif_dock::rosetta_score_then_min_below_thresh, "", -9e9 );
    NEW_OPT(  rif_dock::rosetta_score_at_least, "", -1 );
    NEW_OPT(  rif_dock::rosetta_score_at_most, "", 999999999 );
    NEW_OPT(  rif_dock::rosetta_min_at_least, "", -1 );
    NEW_OPT(  rif_dock::rosetta_min_at_most, "", 999999999 );
    NEW_OPT(  rif_dock::rosetta_min_fraction  , "",  0.1 );
    NEW_OPT(  rif_dock::rosetta_min_targetbb  , "",  false );
    NEW_OPT(  rif_dock::rosetta_min_scaffoldbb  , "",  false );
    NEW_OPT(  rif_dock::rosetta_min_allbb  , "",  false );
    NEW_OPT(  rif_dock::rosetta_min_fix_target, "",  false );
    NEW_OPT(  rif_dock::rosetta_score_cut  , "", -10.0 );
    NEW_OPT(  rif_dock::rosetta_hard_min  , "", false );
    NEW_OPT(  rif_dock::rosetta_score_total  , "", false );
    NEW_OPT(  rif_dock::rosetta_score_ddg_only  , "", false );
    NEW_OPT(  rif_dock::rosetta_score_rifres_rifres_weight, "", 0.75 );
    NEW_OPT(  rif_dock::rosetta_score_rifres_scaffold_weight, "", 0.5 );
    NEW_OPT(  rif_dock::rosetta_soft_score, "", "beta_soft" );
    NEW_OPT(  rif_dock::rosetta_hard_score, "", "beta" );
    
    NEW_OPT(  rif_dock::extra_rotamers, "", true );
    NEW_OPT(  rif_dock::extra_rif_rotamers, "", true );
    NEW_OPT(  rif_dock::always_available_rotamers_level, "", 0 );
    NEW_OPT(  rif_dock::packing_use_rif_rotamers, "", true );
    
    NEW_OPT(  rif_dock::nfold_symmetry, "", 1 );
    NEW_OPT(  rif_dock::symmetry_axis, "", utility::vector1<double>() );
    
    NEW_OPT(  rif_dock::user_rotamer_bonus_constant, "", -2 );
    NEW_OPT(  rif_dock::user_rotamer_bonus_per_chi, "", -2 );
    
    NEW_OPT(  rif_dock::resl0, "", 16 );
    NEW_OPT(  rif_dock::dump_x_frames_per_resl, "Use this to make a movie", 0 );
    NEW_OPT(  rif_dock::dump_only_best_frames, "Only dump the best frames for the movie", false );
    NEW_OPT(  rif_dock::dump_only_best_stride, "When doing dump_only_best_frames, dump every Xth element of the best", 1 );
    NEW_OPT(  rif_dock::dump_prefix, "Convince Brian to make this autocreate the folder", "hsearch" );
    
    NEW_OPT(  rif_dock::scaff_search_mode, "Which scaffold mode and HSearch do you want? Options: default, morph_dive_pop, nineA_baseline", "default");
    NEW_OPT(  rif_dock::nineA_cluster_path, "Path to cluster database for nineA_baseline.", "" );
    NEW_OPT(  rif_dock::nineA_baseline_range, "format cdindex:low-high (python range style)", "");
    
    NEW_OPT(  rif_dock::low_cut_site, "The low cut point for fragment insertion, this res and the previous get minimized.", 0 );
    NEW_OPT(  rif_dock::high_cut_site, "The high cut point for fragment insertion, this res and the next get minimized.", 0 );
    NEW_OPT(  rif_dock::max_insertion, "Maximum number of residues to lengthen protein by.", 0 );
    NEW_OPT(  rif_dock::max_deletion, "Maximum number of residues to shorten protein by.", 0 );
    NEW_OPT(  rif_dock::fragment_cluster_tolerance, "RMSD cluster tolerance for fragments.", 0.5 );
    NEW_OPT(  rif_dock::fragment_max_rmsd , "Max RMSD to starting fragment.", 10000 );
    NEW_OPT(  rif_dock::max_fragments, "Maximum number of fragments to find.", 10000000 );
    NEW_OPT(  rif_dock::morph_rules_files, "List of files for each scaffold to specify morph regions", utility::vector1<std::string>() );
    
    NEW_OPT(  rif_dock::include_parent, "Include parent fragment in diversified scaffolds.", false );
    NEW_OPT(  rif_dock::use_parent_body_energies, "Don't recalculate 1-/2-body energies for fragment insertions", false );
    
    NEW_OPT(  rif_dock::dive_resl , "Dive to this depth before diversifying", 5 );
    NEW_OPT(  rif_dock::pop_resl , "Return to this depth after diversifying", 4 );
    NEW_OPT(  rif_dock::match_this_pdb, "Like tether to input position but applied at diversification time.", "" );
    NEW_OPT(  rif_dock::match_this_rmsd, "RMSD for match_this_pdb", 7 );
    NEW_OPT(  rif_dock::max_beam_multiplier, "Maximum beam multiplier after diversification. Otherwise defaults to number of fragments found.", 1 );
    
    NEW_OPT(  rif_dock::rot_spec_fname,"rot_spec_fname","NOT SPECIFIED");
    // constrain file names
    NEW_OPT(  rif_dock::cst_files, "" , utility::vector1<std::string>() );
    NEW_OPT(  rif_dock::seeding_pos, "" , utility::vector1<std::string>() );
    NEW_OPT(  rif_dock::seeding_by_patchdock, "The format of seeding file can be either Rosetta Xform or raw patchdock outputs", true );
    NEW_OPT(  rif_dock::xform_pos, "" , "" );
    NEW_OPT(  rif_dock::test_longxing, "" , false );
    NEW_OPT(  rif_dock::rosetta_score_each_seeding_at_least, "", -1 );
    NEW_OPT(  rif_dock::cluster_score_cut, "", 0);
    NEW_OPT(  rif_dock::keep_top_clusters_frac, "", 0.5);
    NEW_OPT(  rif_dock::only_dump_scaffold, "" , true );
    NEW_OPT(  rif_dock::requirements,        "which rif residue should be in the final output", utility::vector1< int >());
    
    NEW_OPT(  rif_dock::dump_residue_names,  "which residue should be dumped", utility::vector1<std::string>() );
    NEW_OPT(  rif_dock::dump_fraction,  "the fraction to be dumped", 0.0001 );
    NEW_OPT(  rif_dock::dump_score_cutoff,  "the cutoff value of the score to be dumped", 99999 );
    NEW_OPT(  rif_dock::dump_requirement_num,  "the requirement index", -1 );
    NEW_OPT(  rif_dock::dump_fname,  "the output file name", "random_dump.pdb" );
    
    
}

#endif


#ifndef INCLUDED_rif_dock_test_hh_3
#define INCLUDED_rif_dock_test_hh_3

struct RifDockOpt
{
    std::vector<std::string> scaffold_fnames;
    std::vector<std::string> scaffold_res_fnames;
    std::vector<std::string> data_cache_path;
    std::vector<std::string> rif_files;
    
    bool        VERBOSE                              ;
    double      resl0                                ;
    int64_t     DIM                                  ;
    int64_t     DIMPOW2                              ;
    int64_t     beam_size                            ;
    bool        replace_all_with_ala_1bre            ;
    bool        lowres_sterics_cbonly                ;
    float       tether_to_input_position_cut         ;
    bool        tether_to_input_position             ;
    float       global_score_cut                     ;
    std::string target_pdb                           ;
    std::string outdir                               ;
    std::string output_tag                           ;
    std::string dokfile_fname                        ;
    bool        dump_all_rif_rots                    ;
    bool        dump_all_rif_rots_into_output        ;
    bool        rif_rots_as_chains                   ;
    std::string dump_rifgen_near_pdb                 ;
    float       dump_rifgen_near_pdb_dist            ;
    float       dump_rifgen_near_pdb_frac            ;
    bool        add_native_scaffold_rots_when_packing;
    bool        restrict_to_native_scaffold_res      ;
    float       bonus_to_native_scaffold_res         ;
    float       hack_pack_frac                       ;
    float       hsearch_scale_factor                 ;
    float       search_diameter                      ;
    bool        use_scaffold_bounding_grids          ;
    bool        scaffold_res_use_best_guess          ;
    bool        scaff2ala                            ;
    bool        scaff2alaselonly                     ;
    bool        replace_orig_scaffold_res            ;
    int         require_satisfaction                 ;
    int         require_n_rifres                     ;
    float       target_rf_resl                       ;
    bool        align_to_scaffold                    ;
    bool        output_scaffold_only                 ;
    bool        output_full_scaffold_only            ;
    bool        output_full_scaffold                 ;
    bool        pdb_info_pikaa                       ;
    bool        dump_resfile                         ;
    std::string target_res_fname                     ;
    int         target_rf_oversample                 ;
    float       max_rf_bounding_ratio                ;
    std::string target_rf_cache                      ;
    bool        downscale_atr_by_hierarchy           ;
    float       favorable_1body_multiplier           ;
    float favorable_1body_multiplier_cutoff          ;
    float       favorable_2body_multiplier           ;
    bool        random_perturb_scaffold              ;
    bool        dont_center_scaffold				 ;
    bool        dont_use_scaffold_loops              ;
    bool        cache_scaffold_data                  ;
    float       rf_resl                              ;
    bool        hack_pack                            ;
    int         rf_oversample                        ;
    
    int         rotrf_oversample                     ;
    float       rotrf_resl                           ;
    float       rotrf_spread                         ;
    std::string rotrf_cache_dir                      ;
    float       rotrf_scale_atr                      ;
    
    float       pack_iter_mult                       ;
    int         pack_n_iters                         ;
    float       hbond_weight                         ;
    float       upweight_iface                       ;
    float       upweight_multi_hbond                 ;
    int         n_result_limit                       ;
    float       redundancy_filter_mag                ;
    int         force_output_if_close_to_input_num   ;
    float       force_output_if_close_to_input       ;
    int         n_pdb_out                            ;
    bool        extra_rotamers                       ;
    bool        extra_rif_rotamers                   ;
    int         always_available_rotamers_level      ;
    int         packing_use_rif_rotamers             ;
    
    float       rosetta_score_fraction               ;
    float       rosetta_score_then_min_below_thresh  ;
    float       rosetta_score_at_least               ;
    float       rosetta_score_at_most                ;
    float       rosetta_min_fraction                 ;
    bool        rosetta_min_fix_target               ;
    bool        rosetta_min_targetbb                 ;
    bool        rosetta_min_scaffoldbb               ;
    bool        rosetta_min_allbb                    ;
    float       rosetta_score_cut                    ;
    float       rosetta_hard_min                     ;
    bool        rosetta_score_total                  ;
    bool        rosetta_score_ddg_only               ;
    float       rosetta_score_rifres_rifres_weight   ;
    float       rosetta_score_rifres_scaffold_weight ;
    
    bool        rosetta_beta                         ;
    std::string rosetta_soft_score                   ;
    std::string rosetta_hard_score                   ;
    
    int         nfold_symmetry                       ;
    std::vector<float> symmetry_axis                 ;
    
    float       user_rotamer_bonus_constant		     ;
    float       user_rotamer_bonus_per_chi		     ;
    
    int         dump_x_frames_per_resl				 ;
    bool        dump_only_best_frames                ;
    int         dump_only_best_stride                ;
    std::string dump_prefix                          ;
    
    std::string scaff_search_mode					 ;
    std::string nineA_cluster_path					 ;
    std::string nineA_baseline_range				 ;
    
    int         low_cut_site                         ;
    int         high_cut_site                        ;
    int         max_insertion                        ;
    int         max_deletion                         ;
    float       fragment_cluster_tolerance           ;
    float       fragment_max_rmsd                    ;
    int         max_fragments                        ;
    std::vector<std::string> morph_rules_fnames      ;
    
    bool        include_parent                       ;
    bool        use_parent_body_energies             ;
    
    int         dive_resl                            ;
    int         pop_resl                             ;
    std::string match_this_pdb                       ;
    float       match_this_rmsd                      ;
    float       max_beam_multiplier                  ;
    
    std::string rot_spec_fname                       ;
    // constrain file names
    std::vector<std::string> cst_fnames              ;
    
    std::vector<std::string> seeding_fnames          ;
    std::string              xform_fname;
    bool                    test_longxing;
    float       rosetta_score_each_seeding_at_least  ;
    float       cluster_score_cut                    ;
    float       keep_top_clusters_frac               ;
    bool        only_dump_scaffold                   ;
    bool        seeding_by_patchdock                 ;
    
    std::vector<int> requirements;
    
    
    void init_from_cli();
    
    std::vector< std::string>    dump_residue_names   ;
    float                        dump_fraction        ;
    float                        dump_score_cutoff    ;
    int                          dump_requirement_num      ;
    std::string                  dump_fname           ;
    
    
    
    
};

#endif


#ifndef INCLUDED_rif_dock_test_hh_4
#define INCLUDED_rif_dock_test_hh_4

void RifDockOpt::init_from_cli()
{
    using basic::options::option;
    using namespace basic::options::OptionKeys;
    
    runtime_assert( option[rif_dock::target_rif].user() );
    
    VERBOSE                                = false;
    resl0                                  = option[rif_dock::resl0                              ]();
    DIM                                    = 6;
    DIMPOW2                                = 1<<DIM;
    beam_size                              = int64_t( option[rif_dock::beam_size_M]() * 1000000.0 / DIMPOW2 ) * DIMPOW2;
    replace_all_with_ala_1bre              = option[rif_dock::replace_all_with_ala_1bre          ]();
    
    target_pdb                             = option[rif_dock::target_pdb                         ]();
    lowres_sterics_cbonly                  = option[rif_dock::lowres_sterics_cbonly              ]();
    tether_to_input_position_cut           = option[rif_dock::tether_to_input_position           ]();
    tether_to_input_position               = tether_to_input_position_cut > 0.0;
    global_score_cut                       = option[rif_dock::global_score_cut                   ]();
    outdir                                 = option[rif_dock::outdir                             ]();
    output_tag                             = option[rif_dock::output_tag                         ]();
    dokfile_fname                          = outdir + "/" + option[rif_dock::dokfile             ]();
    dump_all_rif_rots                      = option[rif_dock::dump_all_rif_rots                  ]();
    dump_all_rif_rots_into_output		   = option[rif_dock::dump_all_rif_rots_into_output      ]();
    rif_rots_as_chains                     = option[rif_dock::rif_rots_as_chains                 ]();
    dump_rifgen_near_pdb                   = option[rif_dock::dump_rifgen_near_pdb               ]();
    dump_rifgen_near_pdb_dist              = option[rif_dock::dump_rifgen_near_pdb_dist          ]();
    dump_rifgen_near_pdb_frac              = option[rif_dock::dump_rifgen_near_pdb_frac          ]();
    add_native_scaffold_rots_when_packing  = option[rif_dock::add_native_scaffold_rots_when_packing ]();
    restrict_to_native_scaffold_res        = option[rif_dock::restrict_to_native_scaffold_res       ]();
    bonus_to_native_scaffold_res           = option[rif_dock::bonus_to_native_scaffold_res          ]();
    hack_pack_frac                         = option[rif_dock::hack_pack_frac                        ]();
    hsearch_scale_factor                   = option[rif_dock::hsearch_scale_factor                  ]();
    search_diameter                        = option[rif_dock::search_diameter                       ]();
    use_scaffold_bounding_grids            = option[rif_dock::use_scaffold_bounding_grids           ]();
    scaffold_res_use_best_guess            = option[rif_dock::scaffold_res_use_best_guess           ]();
    scaff2ala                              = option[rif_dock::scaffold_to_ala                       ]();
    scaff2alaselonly                       = option[rif_dock::scaffold_to_ala_selonly               ]();
    replace_orig_scaffold_res              = option[rif_dock::replace_orig_scaffold_res             ]();
    require_satisfaction                   = option[rif_dock::require_satisfaction                  ]();
    require_n_rifres                       = option[rif_dock::require_n_rifres                      ]();
    target_rf_resl                         = option[rif_dock::target_rf_resl                        ]();
    align_to_scaffold                      = option[rif_dock::align_output_to_scaffold              ]();
    output_scaffold_only                   = option[rif_dock::output_scaffold_only                  ]();
    output_full_scaffold_only              = option[rif_dock::output_full_scaffold_only             ]();
    output_full_scaffold                   = option[rif_dock::output_full_scaffold                  ]();
    pdb_info_pikaa                         = option[rif_dock::pdb_info_pikaa                        ]();
    dump_resfile                           = option[rif_dock::dump_resfile                          ]();
    target_res_fname                       = option[rif_dock::target_res                            ]();
    target_rf_oversample                   = option[rif_dock::target_rf_oversample                  ]();
    max_rf_bounding_ratio                  = option[rif_dock::max_rf_bounding_ratio                 ]();
    target_rf_cache                        = option[rif_dock::target_rf_cache                       ]();
    downscale_atr_by_hierarchy             = option[rif_dock::downscale_atr_by_hierarchy            ]();
    favorable_1body_multiplier             = option[rif_dock::favorable_1body_multiplier            ]();
    favorable_1body_multiplier_cutoff      = option[rif_dock::favorable_1body_multiplier_cutoff     ]();
    favorable_2body_multiplier             = option[rif_dock::favorable_2body_multiplier            ]();
    random_perturb_scaffold                = option[rif_dock::random_perturb_scaffold               ]();
    dont_center_scaffold				   = option[rif_dock::dont_center_scaffold					]();
    dont_use_scaffold_loops                = option[rif_dock::dont_use_scaffold_loops               ]();
    cache_scaffold_data                    = option[rif_dock::cache_scaffold_data                   ]();
    rf_resl                                = option[rif_dock::rf_resl                               ]();
    hack_pack                              = option[rif_dock::hack_pack                             ]();
    rf_oversample                          = option[rif_dock::rf_oversample                         ]();
    redundancy_filter_mag                  = option[rif_dock::redundancy_filter_mag                 ]();
    rotrf_oversample                       = option[rif_dock::rotrf_oversample                      ]();
    rotrf_resl                             = option[rif_dock::rotrf_resl                            ]();
    rotrf_spread                           = option[rif_dock::rotrf_spread                          ]();
    rotrf_cache_dir                        = option[rif_dock::rotrf_cache_dir                       ]();
    rotrf_scale_atr                        = option[rif_dock::rotrf_scale_atr                       ]();
    pack_iter_mult                         = option[rif_dock::pack_iter_mult                        ]();
    pack_n_iters                           = option[rif_dock::pack_n_iters                          ]();
    hbond_weight                           = option[rif_dock::hbond_weight                          ]();
    upweight_iface                         = option[rif_dock::upweight_iface                        ]();
    upweight_multi_hbond                   = option[rif_dock::upweight_multi_hbond                  ]();
    n_result_limit                         = option[rif_dock::n_result_limit                        ]();
    redundancy_filter_mag                  = option[rif_dock::redundancy_filter_mag                 ]();
    force_output_if_close_to_input_num     = option[rif_dock::force_output_if_close_to_input_num    ]();
    force_output_if_close_to_input         = option[rif_dock::force_output_if_close_to_input        ]();
    n_pdb_out                              = option[rif_dock::n_pdb_out                             ]();
    extra_rotamers                         = option[rif_dock::extra_rotamers                        ]();
    extra_rif_rotamers                     = option[rif_dock::extra_rif_rotamers                    ]();
    always_available_rotamers_level        = option[rif_dock::always_available_rotamers_level       ]();
    packing_use_rif_rotamers               = option[rif_dock::packing_use_rif_rotamers              ]();
    
  		rosetta_score_fraction                 = option[rif_dock::rosetta_score_fraction                ]();
  		rosetta_score_then_min_below_thresh    = option[rif_dock::rosetta_score_then_min_below_thresh   ]();
  		rosetta_score_at_least                 = option[rif_dock::rosetta_score_at_least                ]();
  		rosetta_score_at_most                  = option[rif_dock::rosetta_score_at_most                 ]();
  		rosetta_min_fraction                   = option[rif_dock::rosetta_min_fraction                  ]();
  		rosetta_min_fix_target                 = option[rif_dock::rosetta_min_fix_target                ]();
  		rosetta_min_targetbb                   = option[rif_dock::rosetta_min_targetbb                  ]();
  		rosetta_min_scaffoldbb                 = option[rif_dock::rosetta_min_scaffoldbb                ]();
  		rosetta_min_allbb                      = option[rif_dock::rosetta_min_allbb                     ]();
  		rosetta_score_cut                      = option[rif_dock::rosetta_score_cut                     ]();
  		rosetta_hard_min                       = option[rif_dock::rosetta_hard_min                      ]();
  		rosetta_score_total                    = option[rif_dock::rosetta_score_total                   ]();
  		rosetta_score_ddg_only                 = option[rif_dock::rosetta_score_ddg_only                ]();
  		rosetta_score_rifres_rifres_weight     = option[rif_dock::rosetta_score_rifres_rifres_weight    ]();
    rosetta_score_rifres_scaffold_weight   = option[rif_dock::rosetta_score_rifres_scaffold_weight  ]();
    rosetta_soft_score                     = option[rif_dock::rosetta_soft_score  					]();
    rosetta_hard_score                     = option[rif_dock::rosetta_hard_score 				    ]();
    rosetta_beta                           = option[corrections::beta 								]();
    user_rotamer_bonus_constant 		   = option[rif_dock::user_rotamer_bonus_constant 			]();
    user_rotamer_bonus_per_chi 			   = option[rif_dock::user_rotamer_bonus_per_chi 			]();
    
    dump_x_frames_per_resl				   = option[rif_dock::dump_x_frames_per_resl                ]();
    dump_only_best_frames				   = option[rif_dock::dump_only_best_frames                 ]();
    dump_only_best_stride                  = option[rif_dock::dump_only_best_stride                 ]();
    dump_prefix                            = option[rif_dock::dump_prefix                           ]();
    
    scaff_search_mode					   = option[rif_dock::scaff_search_mode   				    ]();
    nineA_cluster_path					   = option[rif_dock::nineA_cluster_path                    ]();
    nineA_baseline_range				   = option[rif_dock::nineA_baseline_range                  ]();
    
    low_cut_site                           = option[rif_dock::low_cut_site                          ]();
    high_cut_site                          = option[rif_dock::high_cut_site                         ]();
    max_insertion                          = option[rif_dock::max_insertion                         ]();
    max_deletion                           = option[rif_dock::max_deletion                          ]();
    fragment_cluster_tolerance             = option[rif_dock::fragment_cluster_tolerance            ]();
    fragment_max_rmsd                      = option[rif_dock::fragment_max_rmsd                     ]();
    max_fragments                          = option[rif_dock::max_fragments                         ]();
    
    include_parent                         = option[rif_dock::include_parent                        ]();
    use_parent_body_energies               = option[rif_dock::use_parent_body_energies              ]();
    
    dive_resl                              = option[rif_dock::dive_resl                             ]();
    pop_resl                               = option[rif_dock::pop_resl                              ]();
    match_this_pdb                         = option[rif_dock::match_this_pdb                        ]();
    match_this_rmsd                        = option[rif_dock::match_this_rmsd                       ]();
    max_beam_multiplier                    = option[rif_dock::max_beam_multiplier                   ]();
    
    rot_spec_fname						   = option[rif_dock::rot_spec_fname                        ]();
    
    
    
    
    for( std::string s : option[rif_dock::scaffolds     ]() )     scaffold_fnames.push_back(s);
    for( std::string s : option[rif_dock::scaffold_res  ]() ) scaffold_res_fnames.push_back(s);
    for( std::string s : option[rif_dock::data_cache_dir]() )     data_cache_path.push_back(s);
    
    for( std::string fn : option[rif_dock::target_bounding_xmaps]() ) rif_files.push_back(fn);
    rif_files.push_back( option[rif_dock::target_rif]() );
    
    if( scaff2ala && scaff2alaselonly &&  option[rif_dock::scaffold_to_ala_selonly].user() ){
        std::cout << "WARNING: -scaffold_to_ala overrides -scaffold_to_ala_selonly!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
    }
    
    if( rosetta_score_total && rosetta_score_ddg_only ){
        std::cout << "WARNING: rosetta_score_total overrives rosetta_score_ddg_only" << std::endl;
        rosetta_score_ddg_only = false;
    }
    
    
    
    nfold_symmetry = option[rif_dock::nfold_symmetry]();
    symmetry_axis.clear();
    if( option[rif_dock::symmetry_axis]().size() == 3 ){
        symmetry_axis.push_back( option[rif_dock::symmetry_axis]()[1] );
        symmetry_axis.push_back( option[rif_dock::symmetry_axis]()[2] );
        symmetry_axis.push_back( option[rif_dock::symmetry_axis]()[3] );
    } else if( option[rif_dock::symmetry_axis]().size() == 0 ){
        symmetry_axis.push_back(0);
        symmetry_axis.push_back(0);
        symmetry_axis.push_back(1);
    } else {
        std::cout << "bad rif_dock::symmetry_axis option" << std::endl;
        std::exit(-1);
    }
    
    
    // Brian
    
    
    if (option[rif_dock::use_scaffold_bounding_grids]()) {
        std::cout << "ERROR: use_scaffold_bounding_grids no longer supported. Email bcov@uw.edu" << std::endl;
        std::exit(-1);
    }
    
    
    if (option[rif_dock::nfold_symmetry]() > 1) {
        std::cout << "ERROR: nfold_symmetry not currently supported. Email bcov@uw.edu" << std::endl;
        std::exit(-1);
    }
    
    
    if ( scaff_search_mode == "nineA_baseline" ) {
        if ( scaffold_fnames.size() > 0 ) {
            std::cout << "ERROR: can't use -scaffolds with nineA_baseline." << std::endl;
            std::exit(-1);
        }
        
        std::vector<uint64_t> cdindex_then_clusts = devel::scheme::parse_nineA_baseline_range( nineA_baseline_range );
        uint64_t num_scaffolds = cdindex_then_clusts.size() - 1;
        runtime_assert( num_scaffolds > 0 );
        scaffold_fnames.resize(num_scaffolds);
        
        dont_center_scaffold = true;
    }
    
    for( std::string s : option[rif_dock::morph_rules_files ]() ) morph_rules_fnames.push_back(s);
    
    // constrain file names
    for( std::string s : option[rif_dock::cst_files  ]() ) cst_fnames.push_back(s);
    
    for( std::string s : option[rif_dock::seeding_pos                                              ]() ) seeding_fnames.push_back(s);
    seeding_by_patchdock                    = option[rif_dock::seeding_by_patchdock                ]();
    xform_fname                             = option[rif_dock::xform_pos                           ]();
    test_longxing                           = option[rif_dock::test_longxing                       ]();
    rosetta_score_each_seeding_at_least     = option[rif_dock::rosetta_score_each_seeding_at_least ]();
    cluster_score_cut                       = option[rif_dock::cluster_score_cut                   ]();
    keep_top_clusters_frac                  = option[rif_dock::keep_top_clusters_frac              ]();
    only_dump_scaffold                      = option[rif_dock::only_dump_scaffold                  ]();
    
    
    for( int req : option[rif_dock::requirements]() ) requirements.push_back(req);
    
    
    
    for( std::string res : option[rif_dock::dump_residue_names]() ) dump_residue_names.push_back(res);
    dump_fraction = option[rif_dock::dump_fraction]();
    dump_score_cutoff = option[rif_dock::dump_score_cutoff]();
    dump_requirement_num = option[rif_dock::dump_requirement_num]();
    dump_fname = option[rif_dock::dump_fname]();
    
    
    
    
    
}

#endif


	#include <numeric/random/random.hh>

	#include <ObjexxFCL/format.hh>

	#include <boost/foreach.hpp>
	#include <boost/lexical_cast.hpp>
	// #include <boost/random/mersenne_twister.hpp>

	#include <core/id/AtomID.hh>
	#include <core/import_pose/import_pose.hh>
        #include <core/chemical/ChemicalManager.hh>
        #include <core/chemical/ResidueTypeSet.hh>
	#include <core/pose/Pose.hh>
	#include <core/pose/PDBInfo.hh>
	#include <core/pose/util.hh>
	#include <core/scoring/EnergyGraph.hh>
	#include <core/scoring/ScoreFunction.hh>
	#include <core/scoring/ScoreFunctionFactory.hh>
	#include <core/scoring/hbonds/HBondOptions.hh>
	#include <core/scoring/methods/EnergyMethodOptions.hh>
	#include <core/conformation/ResidueFactory.hh>
	#include <core/kinematics/MoveMap.hh>
	#include <core/scoring/Energies.hh>
        #include <protocols/minimization_packing/MinMover.hh>

	#include <devel/init.hh>
	#include <riflib/RotamerGenerator.hh>
	#include <riflib/rosetta_field.hh>
	#include <riflib/util.hh>
	#include <riflib/rotamer_energy_tables.hh>

	// #include <numeric/alignment/QCP_Kernel.hh>
	#include <parallel/algorithm>
	#include <exception>
	#include <stdexcept>

	#include <scheme/actor/Atom.hh>
	#include <scheme/actor/BackboneActor.hh>
	#include <scheme/actor/VoxelActor.hh>
	#include <scheme/kinematics/Director.hh>
	#include <scheme/kinematics/SceneBase.hh>
	#include <scheme/nest/pmap/OriTransMap.hh>
	#include <scheme/numeric/rand_xform.hh>
	// #include <scheme/objective/ObjectiveFunction.hh>
	#include <scheme/objective/voxel/FieldCache.hh>
	// #include <scheme/objective/voxel/VoxelArray.hh>
	// #include <scheme/objective/hash/XformMap.hh>
	// #include <scheme/objective/storage/RotamerScores.hh>
	#include <scheme/util/StoragePolicy.hh>
	#include <scheme/search/HackPack.hh>
	#include <scheme/objective/integration/SceneObjective.hh>

	#include <riflib/RifFactory.hh>

	#include <utility/file/file_sys_util.hh>
	#include <utility/io/izstream.hh>
	#include <utility/io/ozstream.hh>

	#include <chrono>
	#include <random>


/// Brian
	#include <scheme/objective/hash/XformHash.hh>
	#include <riflib/scaffold/ScaffoldDataCache.hh>
	#include <riflib/scaffold/ScaffoldProviderFactory.hh>


// refactor
	#include <riflib/rifdock_subroutines/util.hh>
	
//#include <riflib/rifdock_subroutines/HSearchFactory.hh>	

	#include <riflib/rifdock_subroutines/hack_pack.hh>
//#include <riflib/rifdock_subroutines/rosetta_rescore.hh>
	#include <riflib/rifdock_subroutines/compile_and_filter_results.hh>
	#include <riflib/rifdock_subroutines/output_results.hh>

	#include <riflib/HSearchConstraints.hh>



using ::scheme::make_shared;
using ::scheme::shared_ptr;

typedef int32_t intRot;


int main(int argc, char *argv[]) {


	register_options();
	devel::init(argc,argv);


	devel::scheme::print_header( "setup global options" );
	RifDockOpt opt;
	opt.init_from_cli();
	utility::file::create_directory_recursive( opt.outdir );



	#ifdef USE_OPENMP
		omp_lock_t cout_lock, dump_lock;
		omp_init_lock( &cout_lock );
		omp_init_lock( &dump_lock );
	#endif


	using namespace core::scoring;
		using std::cout;
		using std::endl;
		using namespace devel::scheme;
		typedef numeric::xyzVector<core::Real> Vec;
		typedef numeric::xyzMatrix<core::Real> Mat;
		// typedef numeric::xyzTransform<core::Real> Xform;
		using ObjexxFCL::format::F;
		using ObjexxFCL::format::I;
		using devel::scheme::print_header;
		using ::devel::scheme::RotamerIndex;

	/////////////////////////////////////////////////////////////////////////////////
	/////////////////////// static shit
	////////////////////////////////////////////////////////////////////////////////
	typedef ::scheme::util::SimpleArray<3,float> F3;
	typedef ::scheme::util::SimpleArray<3,int> I3;


		std::string const rif_type = get_rif_type_from_file( opt.rif_files.back() );
		BOOST_FOREACH( std::string fn, opt.rif_files ){
			std::string rif_type2 = get_rif_type_from_file( fn );
			runtime_assert_msg( rif_type==rif_type2, "mismatched rif types, expect: " + rif_type + " got: " + rif_type2 + " for " + fn );
		}




		std::mt19937 rng( 0);//std::random_device{}() );




		devel::scheme::RifFactoryConfig rif_factory_config;
		rif_factory_config.rif_type = rif_type;
		shared_ptr<RifFactory> rif_factory = ::devel::scheme::create_rif_factory( rif_factory_config );






	print_header( "create rotamer index" );
		
		std::cout << "Loading " << opt.rot_spec_fname << "..." << std::endl;
		std::string rot_index_spec_file = opt.rot_spec_fname;
		shared_ptr< RotamerIndex > rot_index_p = ::devel::scheme::get_rotamer_index( rot_index_spec_file );
		RotamerIndex & rot_index( *rot_index_p );


			print_header( "read in RIFs" ); /////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::vector<shared_ptr<RifBase> > rif_ptrs;
	std::vector<bool> rif_using_rot;
	{
		std::vector<std::string> rif_descriptions( opt.rif_files.size() );
		rif_ptrs.resize( opt.rif_files.size() );
		std::exception_ptr exception = nullptr;
		#ifdef USE_OPENMP
		#pragma omp parallel for schedule(dynamic,1)
		#endif
		for( int i_readmap = opt.rif_files.size() - 1; i_readmap < opt.rif_files.size(); ++i_readmap ){
			if( exception ) continue;
			try {
				std::string const & rif_file = opt.rif_files[i_readmap];
				std::string & rif_dscr = rif_descriptions[i_readmap];
				shared_ptr<RifBase> & rif_ptr = rif_ptrs[i_readmap];
				rif_ptr = rif_factory->create_rif_from_file( rif_file, rif_dscr );
				runtime_assert_msg( rif_ptrs[i_readmap] , "rif creation from file failed! " + rif_file );
				if( opt.VERBOSE ){
					#ifdef USE_OPENMP
					#pragma omp critical
					#endif
					std::cout << "================= read " << rif_file << "=================" << std::endl
					          << "description:" << std::endl << rif_dscr << std::endl
					          << "load factor: " << rif_ptr->load_factor() << std::endl;
				}
				
				if (i_readmap == opt.rif_files.size() -1 ) {
					// #pragma omp criticial
					// rif_ptr->super_print( std::cout, rot_index_p );
					// std::ofstream out_file;
					// out_file.open("rif.txt");
					// rif_ptr->super_print( out_file, rot_index_p );
					// out_file.close();

				}
			} catch( std::exception const & ex ) {
				#ifdef USE_OPENMP
				#pragma omp critical
				#endif
				exception = std::current_exception();
			}
		}
		if( exception ) std::rethrow_exception(exception);

		std::cout << "RIF description:" << std::endl << rif_descriptions.back() << std::endl;
		std::cout << "load factor: " << rif_ptrs.back()->load_factor() << std::endl;
		std::cout << "size of value-type: " << rif_ptrs.back()->sizeof_value_type() << std::endl;
		std::cout << "mem_use: " << ::devel::scheme::KMGT( rif_ptrs.back()->mem_use() ) << std::endl;
		std::cout << "===================================================================================" << std::endl;

		rif_using_rot.resize( rot_index_p->size(), false );
		rif_using_rot[ rot_index.ala_rot() ] = true; // always include ala
		rif_ptrs.back()->get_rotamer_ids_in_use( rif_using_rot );
		int Nusingrot = 0;
		for( int i = 0; i < rif_using_rot.size(); ++i ){
			Nusingrot += rif_using_rot[i] ? 1 : 0;
		}
		std::cout << "rif uses: " << Nusingrot << " rotamers " << std::endl;


		if (opt.dump_rifgen_near_pdb.length() > 0) {
			float dump_dist = opt.dump_rifgen_near_pdb_dist;
			float dump_frac = opt.dump_rifgen_near_pdb_frac;
			core::pose::Pose pose = *(core::import_pose::pose_from_file(opt.dump_rifgen_near_pdb));
			core::conformation::Residue const & res = pose.residue(1);
			numeric::xyzVector<core::Real> n_xyz = res.xyz("N");

			std::stringstream fname;
			fname << "rifgen_dump_" << opt.dump_rifgen_near_pdb << "_" << boost::str(boost::format("%.2f") % dump_dist ) << ".pdb.gz";

			rif_ptrs.back()->dump_rotamers_near_res( res, fname.str(), dump_dist, dump_frac, rot_index_p );

		}

       
		std::vector< std::string > to_dump;
        if (opt.dump_residue_names.size() == 0) {
            to_dump.push_back("*");
        } else {
            to_dump = opt.dump_residue_names;
        }
		
		rif_ptrs.back()->random_dump_rotamers( to_dump, opt.dump_fname, rot_index_p, opt.dump_fraction, opt.dump_score_cutoff, opt.dump_requirement_num );
	}

  
	


			

	



	#ifdef USE_OPENMP
		omp_destroy_lock( &cout_lock );
		omp_destroy_lock( &dump_lock );
	#endif

	std::cout << "rif_dock_test_DONE" << std::endl;

	return 0;
 }
