// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols

#include <riflib/scaffold/util.hh>
#include <riflib/types.hh>
#include <riflib/rifdock_typedefs.hh>
#include <scheme/scaffold/ScaffoldProviderBase.hh>
#include <riflib/scaffold/ScaffoldDataCache.hh>
#include <scheme/numeric/rand_xform.hh>
#include <core/pose/PDBInfo.hh>

// includes to emulate minimize_segment.xml
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/OperateOnResidueSubset.hh>
#include <core/select/residue_selector/NotResidueSelector.hh>
#include <protocols/residue_selectors/StoredResidueSubsetSelector.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <basic/options/option.hh>
#include <protocols/moves/MoverStatus.hh>

#ifdef USEHDF5
#include <basic/options/keys/indexed_structure_store.OptionKeys.gen.hh>
#endif

#include <string>
#include <vector>
#include <boost/any.hpp>

#include <scheme/kinematics/Scene.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <utility/file/file_sys_util.hh>



namespace devel {
namespace scheme {



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
    ) {

    std::string scaff_fname = opt.scaffold_fnames.at(iscaff);
    scafftag = utility::file_basename( utility::file::file_basename( scaff_fname ) );

    std::cout << "!!!!!!!!!!!!!!!name:: " << scaff_fname << std::endl;

    core::import_pose::pose_from_file( scaffold, scaff_fname );

    scaffold_perturb = EigenXform::Identity();
    if( opt.random_perturb_scaffold ){
        runtime_assert_msg( !opt.use_scaffold_bounding_grids,
            "opt.use_scaffold_bounding_grids incompatible with random_perturb_scaffold" );
        std::mt19937 rng( 0);// std::random_device{}() );
        ::scheme::numeric::rand_xform(rng,scaffold_perturb);
        xform_pose( scaffold, eigen2xyz(scaffold_perturb) );
    }


    scaffold_res.clear();
    std::string scaff_res_fname = "";
    if( opt.scaffold_res_fnames.size() ){
        if( opt.scaffold_res_fnames.size() == opt.scaffold_fnames.size() ){
            scaff_res_fname = opt.scaffold_res_fnames.at(iscaff);
        } else if( opt.scaffold_res_fnames.size() == 1 ){
            scaff_res_fname = opt.scaffold_res_fnames.front();
        } else {
            utility_exit_with_message( "-scaffold_res list not same length as -scaffolds list" );
        }
        if( opt.scaffold_res_use_best_guess ){
            utility_exit_with_message("should only use -scaffold_res_use_best_guess true iff not specifying scaffold_res");
        }
        scaffold_res = devel::scheme::get_res( scaff_res_fname , scaffold );
    } else if (opt.scaffold_res_use_best_guess ){
        scaffold_res = devel::scheme::get_designable_positions_best_guess( scaffold, opt.dont_use_scaffold_loops );
        std::cout << "using scaffold residues: ";
        for(auto ir:scaffold_res) std::cout << " " << ir << scaffold.residue(ir).name3();
        std::cout << std::endl;
    } else {
        get_default_scaffold_res( scaffold, scaffold_res );
    }


    if( opt.cst_fnames.size() ){
        std::string cst_fname = "";
        if( opt.cst_fnames.size() == opt.scaffold_fnames.size() ){
            cst_fname = opt.cst_fnames.at(iscaff);
        } else if( opt.cst_fnames.size() == 1 ){
            cst_fname = opt.cst_fnames.front();
        } else {
            utility_exit_with_message( "-cst_files list not same length as -scaffolds list" );
        }
        runtime_assert_msg(parse_constrains_file(cst_fname, csts), "Faild to parse the constrain file!!!");
    }

    if( opt.morph_rules_fnames.size() ){
        std::string morph_rules_fname = "";
        if( opt.morph_rules_fnames.size() == opt.scaffold_fnames.size() ){
            morph_rules_fname = opt.morph_rules_fnames.at(iscaff);
        } else if( opt.morph_rules_fnames.size() == 1 ){
            morph_rules_fname = opt.morph_rules_fnames.front();
        } else {
            utility_exit_with_message( "-morph_rules_files list not same length as -scaffolds list" );
        }
        runtime_assert_msg(parse_morph_rules_file(morph_rules_fname, morph_rules, opt), "Faild to parse the morph_rules file!!!");
    } else {
        std::cout << "Morph rules file not specified, using command-line options" << std::endl;
        morph_rules.push_back(morph_rule_from_options(opt));
    }

}


void
get_default_scaffold_res( core::pose::Pose const & pose,
    utility::vector1<core::Size> & scaffold_res ) {

    for( int ir = 1; ir <= pose.size(); ++ir){
        if( !pose.residue(ir).is_protein() ) continue;
        scaffold_res.push_back(ir);
    }
}

// historically, non_fa was used during HSearch and fa was used during hack pack
ParametricSceneConformationCOP
make_conformation_from_data_cache(ScaffoldDataCacheOP cache, bool fa /*= false*/) {
    typedef numeric::xyzVector<core::Real> Vec;
    ParametricScene scene(1);

    core::pose::Pose const & scaffold_centered = *(cache->scaffold_centered_p);
    std::vector<int> const & scaffres_g2l = *(cache->scaffres_g2l_p);
    utility::vector1<core::Size> const & scaffold_res = *(cache->scaffold_res_p);
    std::vector< SimpleAtom > const & scaffold_simple_atoms = *(cache->scaffold_simple_atoms_p);
    std::vector< SimpleAtom > const & scaffold_simple_atoms_all = *(cache->scaffold_simple_atoms_all_p);


    for( int ir = 1; ir <= scaffold_centered.size(); ++ir ){
        Vec N  = scaffold_centered.residue(ir).xyz("N" );
        Vec CA = scaffold_centered.residue(ir).xyz("CA");
        Vec C  = scaffold_centered.residue(ir).xyz("C" );

        // todo map res indices, must also edit onebody_energies
        BBActor bbactor( N, CA, C, '-', '-', scaffres_g2l[ir-1] );
        runtime_assert( bbactor.index_ == scaffres_g2l[ir-1] );


        if( std::find(scaffold_res.begin(),scaffold_res.end(),ir)!=scaffold_res.end() ){
            scene.add_actor(0,bbactor);
        }
    }

    if (fa) {
        BOOST_FOREACH( SimpleAtom const & sa, scaffold_simple_atoms_all ) scene.add_actor( 0, sa );
        runtime_assert( scene.template num_actors<SimpleAtom>(0) == scaffold_simple_atoms_all.size() );
    } else {
        BOOST_FOREACH( SimpleAtom const & sa, scaffold_simple_atoms ) scene.add_actor( 0, sa );
        runtime_assert( scene.template num_actors<SimpleAtom>(0) == scaffold_simple_atoms.size() );
    }

    ParametricSceneConformationCOP conformation = scene.conformation_ptr(0);

    ParametricSceneConformationOP conformation_mutable = std::const_pointer_cast<ParametricSceneConformation>( conformation );

    uint64_t sanity = cache->debug_sanity + 1;
    cache->debug_sanity = sanity;

    conformation_mutable->cache_data_ = cache;
    conformation_mutable->cache_data_->conformation_is_fa = fa;

    runtime_assert( conformation->cache_data_->debug_sanity == sanity );

    std::cout << "FA status is: " << (fa ? "True" : "False") << std::endl;

    return conformation;


}

#ifdef USEHDF5
std::vector<core::pose::PoseOP>
apply_direct_segment_lookup_mover( 
    protocols::indexed_structure_store::movers::DirectSegmentLookupMover & dsl_mover,
    core::pose::Pose const & in_pose,
    uint64_t low_cut_site,
    uint64_t high_cut_site,
    uint64_t minimum_loop_length,
    uint64_t max_structures,
    uint64_t max_rmsd,
    ApplyDSLMScratch & scratch ) {

    using namespace core::pack::task::operation;
    using namespace core::select::residue_selector;
    using namespace protocols::residue_selectors;


    core::pose::Pose pose;
    utility::vector1<core::Size> rmsd_range;
    core::select::residue_selector::ResidueSubset rmsd_region(in_pose.size(), false);

    for (uint64_t i = 1; i <= in_pose.size(); i++) {
        if (i <= low_cut_site || i >= high_cut_site + 1) {
            pose.append_residue_by_bond( in_pose.residue(i) );
        } else if ( i == high_cut_site ) {
            pose.append_residue_by_jump( in_pose.residue(i), 1, "N", "N", true );
        } else {
            rmsd_region[i] = true;
        }
    }




    const std::string stored_subset_name = "inserted_lookup_segment";
    dsl_mover.stored_subset_name( stored_subset_name );
    dsl_mover.structure_store_path( basic::options::option[basic::options::OptionKeys::indexed_structure_store::fragment_store]() );


    if (scratch.scorefxn_pt.size() == 0) {
        scratch.scorefxn_pt.resize(omp_max_threads());
        scratch.to_ala_pt.resize(omp_max_threads());
        scratch.hardmin_bb_pt.resize(omp_max_threads());
        scratch.rep_scorefxn_pt.resize(omp_max_threads());

        std::cout << "Generating per-thread scaffold relaxers" << std::endl;

        for ( int i = 0; i < omp_max_threads(); i++ ) {

            //<SCOREFXNS>

            core::scoring::ScoreFunctionOP scorefxn = make_shared<core::scoring::ScoreFunction>();
            // core::scoring::ScoreFunctionOP scorefxn = core::scoring::ScoreFunctionFactory::create_score_function("beta_soft");
            scorefxn->set_weight(core::scoring::coordinate_constraint, 2);
            scorefxn->set_weight(core::scoring::cart_bonded, 2);
            scorefxn->set_weight(core::scoring::pro_close, 0);

            //<TASKOPERATIONS>

            // RestrictAbsentCanonicalAASOP ala_only( new RestrictAbsentCanonicalAAS() );
            // ala_only->include_residue( 0 );
            // ala_only->keep_aas( "A" );

            ResidueSelectorCOP stored_residue_subset( new StoredResidueSubsetSelector( stored_subset_name ) );

            ResLvlTaskOperationOP prevent_repacking_rlt( new PreventRepackingRLT() );
            OperateOnResidueSubsetOP only_lookup_segment( new OperateOnResidueSubset( prevent_repacking_rlt, stored_residue_subset, true ) );


            protocols::simple_moves::MutateResidueOP to_ala( new protocols::simple_moves::MutateResidue() );
            to_ala->set_res_name(core::chemical::AA::aa_ala);
            to_ala->set_selector(stored_residue_subset);

            // a = protocols.simple_moves.MutateResidue()

            // a.set_res_name(core.chemical.AA.aa_ala)


            // <MOVERS>

            // core::pack::task::TaskFactoryOP to_ala_tf( new core::pack::task::TaskFactory() );
            // to_ala_tf->push_back( only_lookup_segment );
            // to_ala_tf->push_back( ala_only );

            // protocols::minimization_packing::PackRotamersMoverOP to_ala ( new protocols::minimization_packing::PackRotamersMover() );
            // to_ala->score_function(scorefxn);
            // to_ala->task_factory( to_ala_tf );


            core::pack::task::TaskFactoryOP hardmin_bb_tf( new core::pack::task::TaskFactory() );
            hardmin_bb_tf->push_back( only_lookup_segment );

            protocols::minimization_packing::MinMoverOP min_mover( new protocols::minimization_packing::MinMover() );
            min_mover->tolerance( 0.0001 );
            min_mover->min_type( "lbfgs_armijo_nonmonotone" );
            min_mover->cartesian( true );
            min_mover->score_function( scorefxn );

            protocols::minimization_packing::TaskAwareMinMoverOP hardmin_bb( new protocols::minimization_packing::TaskAwareMinMover( min_mover, hardmin_bb_tf ));
            hardmin_bb->chi( true );
            hardmin_bb->bb( true );


            //<PROTOCOLS>


            core::scoring::ScoreFunctionOP rep_scorefxn = make_shared<core::scoring::ScoreFunction>();
            rep_scorefxn->set_weight(core::scoring::fa_rep, 1);

            scratch.scorefxn_pt[i] = scorefxn;
            scratch.to_ala_pt[i] = to_ala;
            scratch.hardmin_bb_pt[i] = hardmin_bb;
            scratch.rep_scorefxn_pt[i] = rep_scorefxn;

        }
    }

    core::pose::Pose clash_check_reference = pose;
    ::devel::scheme::pose_to_ala( clash_check_reference );

    (*(scratch.rep_scorefxn_pt.front()))(clash_check_reference);
    core::scoring::Energies const & clash_check_reference_energies = clash_check_reference.energies();

    std::vector<core::scoring::EnergiesOP> clash_check_reference_energies_pt(omp_max_threads());
    std::vector<core::pose::PoseOP> in_pose_pt(omp_max_threads());
    for ( int i = 0; i < omp_max_threads(); i++ ) {
        clash_check_reference_energies_pt[i] = clash_check_reference_energies.clone();
        in_pose_pt[i] = in_pose.clone();
    }

    core::Size pose_size = pose.size();

    core::pose::PoseOP result = pose.clone();
    dsl_mover.apply( *result );

    std::vector<core::pose::PoseOP> pre_results;
    pre_results.push_back(result);

    int fetch_size = std::min((uint64_t)(omp_max_threads() * 2), max_structures * 3 / 2);
    int next_fetch = 0;
    bool out_of_results = dsl_mover.get_last_move_status() == protocols::moves::MS_FAIL_DO_NOT_RETRY;
    std::vector<core::pose::PoseOP> mid_results(fetch_size);
    uint64_t found_results = 0;
    uint64_t seen_pre_results = 0;

    while ( ! out_of_results && found_results < max_structures ) {

        std::cout << "Grabbing pre-results" << std::endl;

        next_fetch += fetch_size;
        while (pre_results.size() < next_fetch) {
            result = dsl_mover.get_additional_output();
            if (! result) {
                out_of_results = true;
                break;
            }
            core::pose::PoseOP to_insert = make_shared<core::pose::Pose>();
            to_insert->detached_copy( *result );
            pre_results.push_back(to_insert);
        } 

        int pre_results_size = pre_results.size();

        if (pre_results.size() > max_structures) {
            std::cout << "Found more fragments than max_fragments: " << max_structures << "+" << std::endl;
        } else {
            std::cout << "Only found " << pre_results_size << " fragment insertions" << std::endl;
        }


        mid_results.resize(pre_results.size()); // this should always be growing it
        std::exception_ptr exception = nullptr;

        #ifdef USE_OPENMP
        #pragma omp parallel for schedule(dynamic,1)
        #endif
        for( int ipre = seen_pre_results; ipre < pre_results_size; ++ipre )
        {
            if (found_results >= max_structures) continue;  // burn out the loop when we're done
            try {

                int const ithread = omp_get_thread_num();

                core::pose::PoseOP iresult = pre_results[ipre];

                core::scoring::ScoreFunctionOP iscorefxn = scratch.scorefxn_pt[ithread];
                protocols::simple_moves::MutateResidueOP ito_ala = scratch.to_ala_pt[ithread];
                protocols::minimization_packing::TaskAwareMinMoverOP ihardmin_bb = scratch.hardmin_bb_pt[ithread];
                core::scoring::ScoreFunctionOP irep_scorefxn = scratch.rep_scorefxn_pt[ithread];
                core::scoring::EnergiesOP iclash_check_reference_energies = clash_check_reference_energies_pt[ithread];
                core::pose::PoseOP iin_pose = in_pose_pt[ithread];

                runtime_assert(iresult);
                runtime_assert(ito_ala);
                runtime_assert(ihardmin_bb);
                runtime_assert(irep_scorefxn);
                runtime_assert(iclash_check_reference_energies);
                runtime_assert(iin_pose);
                // runtime_assert();
                // runtime_assert();

                if ( iresult->num_chains() > 1 ) {
                    std::cout << "Broken pose" << std::endl;
                    continue;
                }
                if ( iresult->size() - pose_size < minimum_loop_length ) {
                    std::cout << "Loop too short" << std::endl;
                    continue;
                }
                ito_ala->apply( *iresult );
                ihardmin_bb->apply( *iresult );

                if ( internal_comparative_clash_check( *iclash_check_reference_energies, *iresult, irep_scorefxn,
                        8, low_cut_site - 1, high_cut_site + 1 ) ) {
                    std::cout << "Internal Clash!!" << std::endl;
                    continue;
                }

        // this doesn't handle insertions or deletion!!!!!!!!
                core::Real rmsd = subset_CA_rmsd(*iresult, *iin_pose, rmsd_region, false );
                std::cout << "RMSD: " << rmsd << std::endl;
                // std::cout << rmsd_region << std::endl;
                if ( rmsd > max_rmsd ) {
                    std::cout << "RMSD too high" << std::endl;
                    continue;
                }


                mid_results[ipre] = iresult;

                #pragma omp critical
                found_results ++;

                } catch(...) {
                #pragma omp critical
                exception = std::current_exception();
            }

        } // end of OMP loop
        if( exception ) std::rethrow_exception(exception);

        seen_pre_results = pre_results.size();

    }

    // by doing it this way, we make sure we keep the correct ordering of the results
    std::vector<core::pose::PoseOP> results;
    for ( int i = 0; i < mid_results.size(); i++ ) {
        if ( mid_results[i] ) {
            results.push_back(mid_results[i]);
        }
        if (results.size() >= max_structures) break;
    }


//     do {
//         if ( result->num_chains() > 1 ) {
//             std::cout << "Broken pose" << std::endl;
//             continue;
//         }
//         if ( result->size() - pose.size() < minimum_loop_length ) {
//             std::cout << "Loop too short" << std::endl;
//             continue;
//         }
//         to_ala.apply( *result );
//         hardmin_bb.apply( *result );

//         if ( internal_comparative_clash_check( clash_check_reference_energies, *result, rep_scorefxn,
//                 8, low_cut_site - 1, high_cut_site + 1 ) ) {
//             std::cout << "Internal Clash!!" << std::endl;
//             continue;
//         }

// // this doesn't handle insertions or deletion!!!!!!!!
//         core::Real rmsd = subset_CA_rmsd(*result, in_pose, rmsd_region, false );
//         std::cout << "RMSD: " << rmsd << std::endl;
//         std::cout << rmsd_region << std::endl;
//         if ( rmsd > max_rmsd ) {
//             std::cout << "RMSD too high" << std::endl;
//             continue;
//         }

//         results.push_back( result );
//     } while ( results.size() < max_structures && ( result = dsl_mover.get_additional_output() ) );


    return results;
}
#endif


void
add_pdbinfo_if_missing( core::pose::Pose & pose ) {
    if ( ! pose.pdb_info() ) {
        core::pose::PDBInfoOP pdb_info = make_shared<core::pose::PDBInfo>( pose );
        pose.pdb_info(pdb_info);
    }
}


bool
internal_comparative_clash_check( core::scoring::Energies const & original_energies,
    core::pose::Pose const & to_check,
    core::scoring::ScoreFunctionOP scorefxn,
    core::Real max_fa_rep_delta,
    core::Size low_position,
    core::Size high_position ) {

    core::pose::Pose pose = to_check;
    ::devel::scheme::pose_to_ala(pose);

    (*scorefxn)(pose);

    core::scoring::Energies const & energies = pose.energies();

    for ( core::Size i = low_position; i <= high_position; i++ ) {
        core::Real original = original_energies.residue_total_energies(i)[core::scoring::fa_rep];
        core::Real new_score = energies.residue_total_energies(i)[core::scoring::fa_rep];

        if ( new_score - original > max_fa_rep_delta ) {
            return true;
        }
    }

    return false;

}



}}

