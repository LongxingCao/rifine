
#define GLOBAL_VARIABLES_ARE_BAD
#include <rifine.hh>
#undef GLOBAL_VARIABLES_ARE_BAD

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

// rifdock subroutine util.hh
#include <riflib/rifdock_typedefs_rifine.hh>
#include <riflib/rifdock_subroutines/util_rifine.hh>
#include <riflib/scaffold/ScaffoldDataCache.hh>

#include <utility/file/file_sys_util.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

#include <chrono>
#include <random>
#include <random>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>


int main( int argc, char *argv[] )
{
    #ifdef USE_OPENMP
        omp_lock_t cout_lock, dump_lock;
        omp_init_lock( &cout_lock );
        omp_init_lock( &dump_lock );
    #endif
    
    
    using namespace core::scoring;
    using namespace devel::scheme;
    typedef numeric::xyzVector<core::Real> Vec;
    typedef numeric::xyzMatrix<core::Real> Mat;
    using ObjexxFCL::format::F;
    using ObjexxFCL::format::I;
    using devel::scheme::print_header;
    using devel::scheme::RotamerIndex;
    
    /////////////////////////////////////////////////////////////////////////////////
    /////////////////////// static shit ??
    /////////////////////////////////////////////////////////////////////////////////
    typedef ::scheme::util::SimpleArray< 3, float > F3;
    typedef ::scheme::util::SimpleArray< 3, int > I3;
    
    
    
    register_options();
    devel::init(argc, argv);
    print_header("setup global options");
    RifDockOpt opt;
    opt.init_from_cli();
    utility::file::create_directory_recursive( opt.outdir );
    
    std::string const rif_type = get_rif_type_from_file( opt.rif_files.back() );
    for( std::string fn : opt.rif_files )
    {
        std::string rif_type2 = get_rif_type_from_file( fn );
        runtime_assert_msg( rif_type == rif_type2, "mismatched rif types, expect: " + rif_type + " got: " + rif_type2 + " for " + fn );
    }
    std::cout << "read RIF type: " << rif_type << std::endl;
    
    // It seems the resl0 is hard coded to 16.
    std::cout << "Search Resls: " << opt.resl0;
    std::vector< float > RESLS(1, opt.resl0);
    for (int i = 1; i < opt.rif_files.size(); ++i) {
        RESLS.push_back( RESLS.back()/2.0 );
        std::cout << " " << RESLS.back();
    }
    std::cout << std::endl;
    
    double time_rif = 0, time_pack = 0, time_ros = 0;
    // to generage some random numbers.
    std::mt19937 rng( time(0) );
    boost::uniform_real<> uniform;
    
    // create dokfile
    {
        std::string dokfile_name_orig = opt.dokfile_fname;
        int i = 2;
        while ( utility::file::file_exists(opt.dokfile_fname)) {
            opt.dokfile_fname = dokfile_name_orig + "." + str( i );
            ++i;
        }
        if( i != 2)
            std::cout << "WARNING!" << dokfile_name_orig << " already exists, using " << opt.dokfile_fname << " instead!" << std::endl;
        else
            std::cout << "output scores to " << opt.dokfile_fname << std::endl;
    }
    utility::io::ozstream dokout( opt.dokfile_fname );
    
    // create rif factory
    devel::scheme::RifFactoryConfig rif_factory_config;
    rif_factory_config.rif_type = rif_type;
    shared_ptr< RifFactory > rif_factory = ::devel::scheme::create_rif_factory( rif_factory_config );
    
    // create rotamer index
    print_header("create rotamer index");
    std::cout << "Loading " << opt.rot_spec_fname << "..." << std::endl;
    std::string rot_index_spec_file = opt.rot_spec_fname;
    shared_ptr< RotamerIndex > rot_index_p = ::devel::scheme::get_rotamer_index( rot_index_spec_file );
    RotamerIndex & rot_index ( *rot_index_p );
    
    // pack options
    ::scheme::search::HackPackOpts packopts;
    packopts.pack_n_iters         = opt.pack_n_iters;
    packopts.pack_iter_mult       = opt.pack_iter_mult;
    packopts.hbond_weight         = opt.hbond_weight;
    packopts.upweight_iface       = opt.upweight_iface;
    packopts.upweight_multi_hbond = opt.upweight_multi_hbond;
    packopts.use_extra_rotamers   = opt.extra_rotamers;
    packopts.always_available_rotamers_level = opt.always_available_rotamers_level;
    packopts.packing_use_rif_rotamers = opt.packing_use_rif_rotamers;
    packopts.add_native_scaffold_rots_when_packing = opt.add_native_scaffold_rots_when_packing;
    packopts.rotamer_inclusion_threshold = -0.5;//-0.5
    packopts.rotamer_onebody_inclusion_threshold = 5;//5
    packopts.init_with_best_1be_rots = true;
    packopts.user_rotamer_bonus_constant=opt.user_rotamer_bonus_constant;
    packopts.user_rotamer_bonus_per_chi=opt.user_rotamer_bonus_per_chi;
    
    // rotamer rosetta field operations.
    RotamerRFOpts rotrfopts;
    rotrfopts.oversample     = opt.rotrf_oversample;
    rotrfopts.field_resl     = opt.rotrf_resl;
    rotrfopts.field_spread   = opt.rotrf_spread;
    rotrfopts.data_dir       = opt.rotrf_cache_dir;
    rotrfopts.scale_atr      = opt.rotrf_scale_atr;
    ::devel::scheme::RotamerRFTablesManager rotrf_table_manager( rot_index_p, rotrfopts );
    // rotrf_table_manager.preinit_all();
    
    // make two body options.
    MakeTwobodyOpts make2bopts;
    // hacked by brian             VVVV
    // Why set the onebody_so high
    make2bopts.onebody_threshold = 4.0;
    make2bopts.distance_cut = 15.0;
    make2bopts.hbond_weight = packopts.hbond_weight;
    make2bopts.favorable_2body_multiplier = opt.favorable_2body_multiplier;
    
    // read and prepare target structure
    print_header( "read and prepare target structure ");
    core::pose::Pose target;
    std::vector< SimpleAtom > target_simple_atoms;
    utility::vector1< core::Size > target_res;
    std::vector< HBondRay > target_donors, target_acceptors;
    float rif_radius = 0.0, target_redundancy_filter_rg = 0.0;
    {
        core::import_pose::pose_from_file( target, opt.target_pdb );
        
        target_res = devel::scheme::get_res( opt.target_res_fname, target, /*nocgp*/ false );
        get_rg_radius( target, target_redundancy_filter_rg, rif_radius, target_res, true );
        rif_radius += 7; //the maximum distance between the atom and the center of the target residues, plus 7!
        
        ::devel::scheme::HBRayOpts hbopt;
        hbopt.satisfied_atoms = ::devel::scheme::get_satisfied_atoms( target );
        for ( core::Size ir : target_res )
        {
            ::devel::scheme::get_donor_rays     ( target, ir, hbopt, target_donors );
            ::devel::scheme::get_acceptor_rays  ( target, ir, hbopt, target_acceptors );
        }
        std::cout << "target_donors.size(): " << target_donors.size() << "    target_acceptors.size(): " << target_acceptors.size() << std::endl;
    }
    
    // rosetta field
    std::vector< VoxelArrayPtr > target_field_by_atype;
    std::vector< std::vector< VoxelArrayPtr > > target_bounding_by_atype;
    {
        // target field by atype
        target_bounding_by_atype.resize( RESLS.size() );
        devel::scheme::RosettaFieldOptions rfopts;
        rfopts.data_dir = "DUMMY_DATA_DIR_FIXME";
        rfopts.oversample = opt.target_rf_oversample;
        rfopts.block_hbond_sites = false;
        rfopts.max_bounding_ratio = opt.max_rf_bounding_ratio;
        rfopts.fail_if_no_cached_data = true;
        rfopts.repulsive_only_boundary = true;
        rfopts.cache_mismatch_tolerance = 0.01; // this is kinda loose...
        std::string cache_prefix = opt.target_rf_cache;
        devel::scheme::get_rosetta_fields_specified_cache_prefix(
            cache_prefix,
            opt.target_pdb,
            target,
            target_res,
            rfopts,
            target_field_by_atype,
            false
            );
        
        // std::cout << "using target bounding grids, generating (or loading) them" << std::endl;
        if ( true ) {
            devel::scheme::RosettaFieldOptions rfopts;
            rfopts.field_resl = opt.target_rf_resl;
            rfopts.data_dir = "DUMMY_DATA_DIR_FIXME";
            rfopts.oversample = opt.target_rf_oversample;
            rfopts.block_hbond_sites = false;
            rfopts.max_bounding_ratio = opt.max_rf_bounding_ratio;
            rfopts.fail_if_no_cached_data = true;
            rfopts.repulsive_only_boundary = true; // default
            devel::scheme::get_rosetta_bounding_fields_from_fba(
                RESLS,
                opt.target_pdb,
                target,
                target_res,
                rfopts,
                target_field_by_atype,
                target_bounding_by_atype,
                false,
                cache_prefix
            );
            runtime_assert( target_bounding_by_atype.size() == RESLS.size() );
            // now scale down the any positive component by 1/RESL if RESL > 1
            if ( opt.downscale_atr_by_hierarchy ) {
                std::cout << "downscale_atr_by_hierarchy on target bounding steric grids" << std::endl;
                for ( int iresl = 0; iresl < RESLS.size(); ++iresl ) {
                    float correction = 1.0 / RESLS[iresl];
                    if (correction > 1.0 ) break;
                    for ( VoxelArrayPtr vap: target_bounding_by_atype[iresl] ) {
                        if ( vap != nullptr ) {
                            std::exception_ptr exception = nullptr;
                            #ifdef USE_OPENMP
                            #pragma omp parallel for schedule(dynamic,64)
                            #endif
                            for ( int k = 0; k < vap->num_elements(); ++k ) {
                                if( exception ) continue;
                                try {
                                    float & dat = vap->data()[k];
                                    if( dat < 0 ){
                                        dat = dat * correction;
                                        // dat = 0; // for testing w/o attractive sterics
                                    }
                                } catch( std::exception const & ex ) {
                                    #pragma omp critical
                                    exception = std::current_exception();
                                } // try ... catch ...
                            } // if val != nullptr
                            if( exception ) std::rethrow_exception(exception);
                        } // for each voxel array
                    }
                } // for each resolution
            } // downscale_atr_by_hierarchy
            
        } // end of if(true)
    } // end of load bounding grid
    
    // prepare and read RIF
    print_header( "read in RIFs" );
    std::vector< shared_ptr< RifBase> > rif_ptrs;
    std::vector< bool > rif_using_rot;
    {
        std::vector< std::string > rif_descriptions( opt.rif_files.size() );
        rif_ptrs.resize( opt.rif_files.size() );
        std::exception_ptr exception = nullptr;
        #ifdef USE_OPENMP
        #pragma omp parallel for schedule( dynamic, 1)
        #endif
        for ( int i_readmap = 0; i_readmap < opt.rif_files.size(); ++i_readmap ) {
            if ( exception ) continue;
            try {
                std::string const & rif_file = opt.rif_files[i_readmap];
                std::string & rif_decr = rif_descriptions[i_readmap];
                shared_ptr<RifBase> & rif_ptr = rif_ptrs[i_readmap];
                rif_ptr = rif_factory->create_rif_from_file( rif_file, rif_decr );
                runtime_assert_msg( rif_ptrs[i_readmap] , "rif creation from file failed! " + rif_file );
                if( opt.VERBOSE ){
                    #ifdef USE_OPENMP
                    #pragma omp critical
                    #endif
                    std::cout << "================= read " << rif_file << "=================" << std::endl
                    << "description:" << std::endl << rif_decr << std::endl
                    << "load factor: " << rif_ptr->load_factor() << std::endl;
                }
                #ifdef USE_OPENMP
                #pragma omp critical
                #endif
                std::cout << "loaded RIF score for resl " << F(7,3,RESLS[i_readmap])
                          << " raw cart_resl: " << F(7,3,rif_ptr->cart_resl() )
                          << " raw ang_resl: " << F(7,3,rif_ptr->ang_resl() ) << std::endl;
            } catch ( std::exception const & ex ) {
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
        rif_using_rot[ rot_index.ala_rot() ] = true; // always include ala ??
        rif_ptrs.back()->get_rotamer_ids_in_use( rif_using_rot );
        int Nusingrot = 0;
        for ( int i = 0; i < rif_using_rot.size(); ++i ) {
            Nusingrot += rif_using_rot[i] ? 1 : 0;
        }
        std::cout << "rif uses: " << Nusingrot << " rotamers " << std::endl;
    }
    
    if ( 0 == opt.scaffold_fnames.size() ){
        std::cout << "WARNING: NO SCAFFOLDS!!!!!!!!!!"  << std::endl;
    }
    
    for ( int iscaff = 0; iscaff < opt.scaffold_fnames.size(); ++iscaff )
    {
        std::string scaff_fname = opt.scaffold_fnames.at( iscaff );
        std::vector< std::string > scaffold_sequence_glob0; // Scaffold sequence in name3 space
        utility::vector1< core::Size > scaffold_res; // seqpose of residue to design, default whold scaffold
        
        try {
            // the docking process happens here!!!
            runtime_assert( rot_index_p );
            std::string scafftag = utility::file_basename( utility::file::file_basename( scaff_fname ) );
            std::cout << "/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl;
            std::cout << "/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl;
            std::cout << "//////   begin scaffold " << scafftag << " " << iscaff << " of " << opt.scaffold_fnames.size() << std::endl;
            std::cout << "/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl;
            std::cout << "/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl;
            
            
            core::pose::Pose scaffold, scaffold_centered, scaffold_full_centered, both_pose, both_full_pose, scaffold_only_pose;
            float scaff_radius = 0.0;
            float redundancy_filter_rg = 0.0;
            
            std::vector< int > scaffres_g2l, scaffres_l2g;
            std::vector< bool > scaffuseres;
            Eigen::Vector3f scaffold_center;
            std::vector< Vec > scaffca;
            std::vector< std::vector< float > > scaffold_onebody_glob0;
            shared_ptr<std::vector<std::vector<float> > > local_onebody_p = make_shared<std::vector<std::vector<float> > >();
            std::vector< std::vector< float > > & local_onebody ( *local_onebody_p );
            shared_ptr<std::vector< std::pair<int,int> > > local_rotamers_p = make_shared<std::vector< std::pair<int,int>>>();;
            std::vector< std::pair< int, int > > & local_rotamers( *local_rotamers_p );
            typedef ::scheme::objective::storage::TwoBodyTable< float > TBT;
            shared_ptr<TBT> scaffold_twobody = make_shared< TBT > ( scaffold.size(), rot_index.size() );
            shared_ptr<TBT> local_twobody;
            EigenXform scaffold_perturb = EigenXform::Identity();
            
            // scaffold initialization, scaffold residues, scaffold_twobody, onebody /////////////////////////
            {
                core::import_pose::pose_from_file( scaffold, scaff_fname );
                if( opt.random_perturb_scaffold ){
                    if( opt.random_perturb_scaffold ){
                        runtime_assert_msg( !opt.use_scaffold_bounding_grids,
                            "opt.use_scaffold_bounding_grids incompatible with random_perturb_scaffold" );
                    }
                    ::scheme::numeric::rand_xform( rng, scaffold_perturb );
                    xform_pose( scaffold, eigen2xyz(scaffold_perturb ) );
                }
                
                scaffold_full_centered = scaffold;
                
                for ( int ir = 1; ir <= scaffold.size(); ++ir ) {
                    scaffold_sequence_glob0.push_back( scaffold.residue(ir).name3() );
                }
                
                // scaffold residues
                std::string scaff_res_fname = "";
                if ( opt.scaffold_res_fnames.size() ) {
                    if ( opt.scaffold_res_fnames.size() == opt.scaffold_fnames.size() ) {
                        scaff_res_fname = opt.scaffold_res_fnames.at( iscaff );
                    } else if ( opt.scaffold_res_fnames.size() == 1 ) {
                        scaff_res_fname = opt.scaffold_res_fnames.front();
                    } else {
                        utility_exit_with_message( "-scaffold_res list not the same length as -scaffolds list" );
                    }
                    if ( opt.scaffold_res_use_best_guess ) {
                        utility_exit_with_message( "should only use -scaffold_res_use_best_guess true iff not specifying scaffold_res" );
                    }
                    scaffold_res = devel::scheme::get_res( scaff_res_fname, scaffold );
                } else if ( opt.scaffold_res_use_best_guess ) {
                    scaffold_res = devel::scheme::get_designable_positions_best_guess( scaffold, opt.dont_use_scaffold_loops );
                    std::cout << "using scaffold residues: ";
                    for ( auto ir : scaffold_res ) std::cout << " " << ir << scaffold.residue(ir).name3();
                    std::cout << std::endl;
                } else {
                    for( int ir = 1; ir <= scaffold.size(); ++ir){
                        if( !scaffold.residue(ir).is_protein() ) continue;
                        if( scaffold.residue(ir).name3() == "PRO" ) continue;
                        if( scaffold.residue(ir).name3() == "GLY" ) continue;
                        if( scaffold.residue(ir).name3() == "CYS" ) continue;
                        scaffold_res.push_back(ir);
                    }
                } // scaffold residues selection
                
                if ( opt.scaff2ala )              ::devel::scheme::pose_to_ala( scaffold );
                else if ( opt.scaff2alaselonly )  ::devel::scheme::pose_to_ala( scaffold, scaffold_res );
                std::cout << "rifdock scaffold_res: " << scaffold_res << std::endl;
                
                
                float scaff_redundancy_filter_rg = 0;
                get_rg_radius( scaffold, scaff_redundancy_filter_rg, scaff_radius, scaffold_res, false); // not allatom for scaffold.
                redundancy_filter_rg = std::min( scaff_redundancy_filter_rg, target_redundancy_filter_rg );
                std::cout << "scaffold selected region rg: " << scaff_redundancy_filter_rg << ", radius: " << scaff_radius << std::endl;
                std::cout << "using redundancy_filter_rg: " << redundancy_filter_rg << std::endl;
                
                
                // scaffold global_to_local, local_to_global
                int count = 0;
                scaffres_g2l.resize( scaffold.size(), -1   );
                scaffuseres.resize ( scaffold.size(), false);
                for ( auto ir : scaffold_res ){
                    scaffres_g2l[ir-1] = count++;
                    scaffres_l2g.push_back(ir-1);
                    scaffuseres[ir-1] = true;
                }
                
                // test score before and after
                
                
                // scaffold center and move the scaffold to the center. ########## setup the pose..
                scaffold_center = pose_center( scaffold, scaffold_res );
                scaffold_centered = scaffold;
                scaffold_full_centered = scaffold;
                
                for (int ir = 1; ir <= scaffold.size(); ++ir ) {
                    Vec tmp( scaffold_center[0], scaffold_center[1], scaffold_center[2] );
                    for ( int ia = 1; ia <= scaffold.residue_type(ir).natoms(); ++ia) {
                        core::id::AtomID aid(ia, ir);
                        scaffold_centered.set_xyz( aid, scaffold.xyz(aid) - tmp );
                    }
                    for( int ia = 1; ia <= scaffold_full_centered.residue_type(ir).natoms(); ++ia ){
                        core::id::AtomID aid(ia,ir);
                        scaffold_full_centered.set_xyz( aid, scaffold_full_centered.xyz(aid) - tmp );
                    }
                }
                
                add_pdbinfo_if_missing( scaffold );
                add_pdbinfo_if_missing( scaffold_centered );
                add_pdbinfo_if_missing( scaffold_full_centered );
                
                both_pose      = scaffold_centered;
                both_full_pose = scaffold_full_centered;
                scaffold_only_pose = scaffold_centered;
                ::devel::scheme::append_pose_to_pose( both_pose, target );
                ::devel::scheme::append_pose_to_pose( both_full_pose, target );
                runtime_assert( both_pose.size() == scaffold.size() + target.size() );
                runtime_assert( both_pose.size() == both_full_pose.size() );

                // fill the scaffca vector using the CA coordinates.
                for (int ir = 1; ir <= scaffold.size(); ++ir ) {
                    scaffca.push_back( scaffold.residue(ir).xyz("CA") );
                }
                
                // setup onebody table: onebody_global and onebody_local.
                std::string scaff_tag = utility::file_basename( utility::file::file_basename( scaff_fname ) );
                std::string scaff_res_hashstr = ::devel::scheme::get_res_list_hash( scaffold_res );
                std::string cachefile_1be = "__1BE_"+scaff_tag+(opt.replace_all_with_ala_1bre?"_ALLALA":"")+"_reshash"+scaff_res_hashstr+".bin.gz";
                if ( ! opt.cache_scaffold_data ) cachefile_1be = "";
                std::cout << "rifdock: get_onebody_rotamer_energies" << std::endl;
                
                
                get_onebody_rotamer_energies(
                        scaffold_centered,
                        scaffold_res,
                        rot_index,
                        scaffold_onebody_glob0,
                        opt.data_cache_path,
                        cachefile_1be,
                        opt.replace_all_with_ala_1bre
                    );
                if ( opt.restrict_to_native_scaffold_res ) {
                    std::cout << "KILLING NON-NATIVE ROTAMERS ON SCAFFOLD!!!" << std::endl;
                    for ( int ir = 0; ir < scaffold_onebody_glob0.size(); ++ir ){
                        for ( int irot = 0; irot < rot_index.size(); ++irot ) {
                            if ( rot_index.resname(irot) != scaffold_sequence_glob0.at(ir) && rot_index.resname(irot) != "ALA" ) {
                                scaffold_onebody_glob0[ir][irot] = 9e9;
                            }
                        }
                    }
                }
                if( opt.bonus_to_native_scaffold_res != 0 ){
                    std::cout << "adding to native scaffold res 1BE " << opt.bonus_to_native_scaffold_res << std::endl;
                    for( int ir = 0; ir < scaffold_onebody_glob0.size(); ++ir ){
                        for( int irot = 0; irot < rot_index.size(); ++irot ){
                            if( rot_index.resname(irot) == scaffold_sequence_glob0.at(ir) ){
                                scaffold_onebody_glob0[ir][irot] += opt.bonus_to_native_scaffold_res;
                            }
                        }
                    }
                }
                
                
                // move the onebody multiplier here.
                if ( opt.favorable_1body_multiplier != 1 ) {
                    for( int ir = 0; ir < scaffold_onebody_glob0.size(); ++ir ){
                        for( int irot = 0; irot < rot_index.size(); ++irot ){
                            if ( scaffold_onebody_glob0[ir][irot] <= opt.favorable_1body_multiplier_cutoff ) {
                                scaffold_onebody_glob0[ir][irot] *= opt.favorable_1body_multiplier;
                            }
                        }
                    }
                }
                // onebody multiplyer
        
                
                
                for (int i =0; i < scaffres_l2g.size(); ++i) {
                    local_onebody.push_back( scaffold_onebody_glob0.at( scaffres_l2g.at(i) ) );
                }
                for( int i = 0; i < scaffres_g2l.size(); ++i ){
                    if( scaffres_g2l[i] < 0 ){
                        for( float & f : scaffold_onebody_glob0[i] ) f = 9e9;
                    }
                }
                
                // setup twobody table: twobody_global and twobody_local
                std::cout << "rifdock: get_twobody_tables" << std::endl;
                std::string cachefile2b = "__2BE_" + scaff_tag + "_reshash" + scaff_res_hashstr + ".bin.gz";
                if( ! opt.cache_scaffold_data || opt.extra_rotamers ) cachefile2b = "";
                std::string dscrtmp;
                get_twobody_tables(
                    opt.data_cache_path,
                    cachefile2b,
                    dscrtmp,
                    scaffold,
                    rot_index,
                    scaffold_onebody_glob0,
                    rotrf_table_manager,
                    make2bopts,
                    *scaffold_twobody
                );
                std::cout << "rifdock: twobody memuse: " << (float)scaffold_twobody->twobody_mem_use()/1000.0/1000.0 << "M" << std::endl;
                {
                    std::cout << "rifdock: onebody dimension: " << scaffold_onebody_glob0.size() << " " << scaffold_onebody_glob0.front().size() << std::endl;
                    int onebody_n_allowed = 0;
                    for( auto const & t : scaffold_onebody_glob0 ){
                        for( auto const & v : t ){
                            if( v < make2bopts.onebody_threshold ) onebody_n_allowed++;
                        }
                    }
                    std::cout << "rifdock: onebody Nallowed: " << onebody_n_allowed << std::endl;
                }
                local_twobody = scaffold_twobody->create_subtable( scaffuseres, scaffold_onebody_glob0, make2bopts.onebody_threshold );
                std::cout << "filt_2b memuse: " << (float)local_twobody->twobody_mem_use()/1000.0/1000.0 << "M" << std::endl;
                
                
                
                // store the local rotamers, the rotamer bounds.
                local_rotamers.clear();
                for( int i = 0; i < local_onebody.size(); ++i ) {
                    int iresglocal = scaffres_l2g.at(i);
                    std::string name3 = scaffold_sequence_glob0.at(iresglocal);
                    std::pair< int, int > ib = rot_index.index_bounds( name3 );
                    local_rotamers.push_back( ib );
                }
                
                
            } // scaffold residues, scaffold_twobody, onebody table initialization  //////////////////////////
            
            // scaffold simple atoms ??????? There are used to set the bb-actor?? Mainly for target-scaffold interaction.
            std::vector< SimpleAtom > scaffold_simple_atoms, scaffold_simple_atoms_all;
            {
                for (int ir = 1; ir <= scaffold_centered.size(); ++ir ) {
                    utility::vector1<core::Size> resids( 1, ir );
                    {
                        std::vector<SchemeAtom> scaff_res_atoms;
                        if (!opt.lowres_sterics_cbonly && std::find( scaffold_res.begin(), scaffold_res.end(), ir) != scaffold_res.end() ) {
                            // backbone only
                            devel::scheme::get_scheme_atoms( scaffold_centered, resids, scaff_res_atoms, true );
                        } else {
                            // CB only
                            devel::scheme::get_scheme_atoms_cbonly( scaffold_centered, resids, scaff_res_atoms );
                        }
                        int restype = rot_index.chem_index_.resname2num( scaffold_centered.residue(ir).name3() );
                        for( int ia = 0; ia < scaff_res_atoms.size(); ++ia){
                            SchemeAtom const & a( scaff_res_atoms[ia] );
                            runtime_assert( a.type() > 0 );
                            if( a.type() >= 21 ) continue;
                            SimpleAtom sa( a.position(), a.type(), restype, ia );
                            scaffold_simple_atoms.push_back(sa);
                        }
                    }
                    {
                        std::vector<SchemeAtom> all_scaff_res_atoms;
                        devel::scheme::get_scheme_atoms( scaffold_centered, resids, all_scaff_res_atoms, false );
                        int restype = rot_index.chem_index_.resname2num( scaffold_centered.residue(ir).name3() ); // for UNK will be -1
                        for( int ia = 0; ia < all_scaff_res_atoms.size(); ++ia){
                            SchemeAtom const & a( all_scaff_res_atoms[ia] );
                            runtime_assert( a.type() > 0 );
                            if( a.type() >= 21 ) continue;
                            SimpleAtom sa( a.position(), a.type(), restype, ia );
                            scaffold_simple_atoms_all.push_back(sa);
                        }
                    }
                }
                std::cout << "scaffold_simple_atoms " << scaffold_simple_atoms.size() << std::endl;
            }
            
            
            // RifSceneObjective config. To create scene and objective.
            RifSceneObjectiveConfig rso_config;
                rso_config.packopts = &packopts;
                rso_config.rif_ptrs = rif_ptrs;
                rso_config.target_bounding_by_atype = &target_bounding_by_atype;
                rso_config.target_field_by_atype = &target_field_by_atype;
                rso_config.rot_index_p = rot_index_p;
                rso_config.target_donors = &target_donors;
                rso_config.target_acceptors = &target_acceptors;
                rso_config.n_sat_groups = 1000;//target_donors.size() + target_acceptors.size();
                rso_config.require_satisfaction = opt.require_satisfaction;
                rso_config.require_n_rifres = opt.require_n_rifres;
            
            ScenePtr scene_prototype;
            std::vector< ObjectivePtr > objectives;
            ObjectivePtr packing_objective;
            runtime_assert( rif_factory->create_objectives( rso_config, objectives, packing_objective ) );
            scene_prototype = rif_factory->create_scene();
            runtime_assert_msg( objectives.front()->is_compatible( *scene_prototype ), "objective and scene types not compatible!" );
            
            
            // setup scene_minimal and scene_full. The scene_minimal is used during the Hsearch and scene_full is for Hpack.
            // it seems the scene_full is never used.
            ScenePtr scene_minimal( scene_prototype->clone_deep() );
            ScenePtr scene_full( scene_prototype->clone_deep() );
            {
                // setup the bb actor
                for (int ir = 1; ir <= scaffold.size(); ++ir) {
                    Vec N = scaffold_centered.residue(ir).xyz("N");
                    Vec CA = scaffold_centered.residue(ir).xyz("CA");
                    Vec C  = scaffold_centered.residue(ir).xyz("C" );
                    
                    // todo map res indices, must also edit onebody_energies
                    BBActor bbactor( N, CA, C, '-', '-', scaffres_g2l[ir-1] );
                    runtime_assert( bbactor.index_ == scaffres_g2l[ir-1] );
                    
                    if( std::find(scaffold_res.begin(),scaffold_res.end(),ir)!=scaffold_res.end() ){
                        scene_full->add_actor(1,bbactor);
                        scene_minimal->add_actor(1,bbactor);
                    }
                }
                
                // setup the backbones.
                {
                    for( SimpleAtom const & sa : scaffold_simple_atoms ) scene_minimal->add_actor( 1, sa );
                    for( SimpleAtom const & sa : scaffold_simple_atoms_all ) scene_full->add_actor( 1, sa );
                    runtime_assert( scene_minimal->template num_actors<SimpleAtom>(1) == scaffold_simple_atoms.size() );
                    runtime_assert( scene_full->template num_actors<SimpleAtom>(1) == scaffold_simple_atoms_all.size() );
                    scene_minimal->add_actor( 0, VoxelActor(target_bounding_by_atype) );
                    scene_full->add_actor( 0, VoxelActor(target_bounding_by_atype) );
                }
            }
            
            
            // To make is compatible with Brian's #include <riflib/scaffold/ScaffoldDataCache.hh>
            ParametricSceneConformationCOP conformation_minimal = (dynamic_cast<ParametricScene*>(scene_minimal.get()))->conformation_ptr(1);
            ParametricSceneConformationOP conformation_minimal_mutable = std::const_pointer_cast<ParametricSceneConformation>( conformation_minimal );
            ScaffoldDataCacheOP cache_minimal = make_shared< ScaffoldDataCache >();
            cache_minimal->local_onebody_p = local_onebody_p;
            cache_minimal->local_twobody_p = local_twobody;
            cache_minimal->local_rotamers_p = local_rotamers_p;
            cache_minimal->conformation_is_fa = false;
            conformation_minimal_mutable->cache_data_ = cache_minimal;
            
            
            ParametricSceneConformationCOP conformation_full = (dynamic_cast<ParametricScene*>(scene_full.get()))->conformation_ptr(1);
            ParametricSceneConformationOP conformation_full_mutable = std::const_pointer_cast<ParametricSceneConformation>( conformation_full );
            ScaffoldDataCacheOP cache_full = make_shared< ScaffoldDataCache >();
            cache_full->local_onebody_p = local_onebody_p;
            cache_full->local_twobody_p = local_twobody;
            cache_full->local_rotamers_p = local_rotamers_p;
            cache_full->conformation_is_fa = true;
            conformation_full_mutable->cache_data_ = cache_full;
            
            
            // parse the seeding files.
            std::string seeding_fname = "";
            shared_ptr<std::vector<EigenXform> > seeding_positions_p = make_shared<std::vector<EigenXform> >();
            std::vector<EigenXform> & seeding_positions(*seeding_positions_p);
            bool use_seeding = false;
            {
                if( opt.seeding_fnames.size() ){
                    if( opt.seeding_fnames.size() == opt.scaffold_fnames.size() ){
                        seeding_fname = opt.seeding_fnames.at(iscaff);
                    }  else {
                        utility_exit_with_message( "-seeding_files list not same length as -scaffolds list" );
                    }
                    runtime_assert_msg(parse_seeding_file(seeding_fname, seeding_positions), "Faild to parse the seeding file!!!");
                    use_seeding = true;
                }else {
                    use_seeding = false;
                }
            }
            
            
            // scores for the initial positions.
            if ( opt.test_longxing )
            {
                std::cout << "scores for the initial seeding positions" << std::endl;
                if ( use_seeding ) {
                    for ( int i = 0; i < seeding_positions.size(); ++i ) {
                        scene_minimal->set_position( 1, seeding_positions[i] );
                        // output the score of each position:      Seeding position #i: sc sc sc sc sc
                        std::cout << "Seeding position " << I(5,i + 1) << ":";
                        for( int j = 0; j < RESLS.size(); ++j ){
                            double sc = objectives[j]->score(*scene_minimal);
                            std::cout << " " << F(10,3,sc);
                        }
                        std::cout << std::endl;
                    }
                }
            }
            
            
            
            // test case for the original position.
            if ( opt.test_longxing )
            {
                std::cout << "scores for scaffold in original position: " << std::endl;
                EigenXform x(EigenXform::Identity());
                x.translation() = scaffold_center;
                scene_minimal->set_position(1,x);
                for(int i = 0; i < RESLS.size(); ++i){
                    std::vector<float> sc = objectives[i]->scores(*scene_minimal);
                    std::cout << "input bounding score " << i << " " << F(7,3,RESLS[i]) << " "
                    << F( 7, 3, sc[0]+sc[1] ) << " "
                    << F( 7, 3, sc[0]       ) << " "
                    << F( 7, 3, sc[1]       ) << std::endl;
                }
                
            }
            
            // I would only do a refinement of the seeding positions.
            // This is fixed, never change this part so that it matchs the xform_positions.
            DirectorBase director;
            shared_ptr<RifDockNestDirector> n_director;
            {
                // start resolution from 1... must be compatible with the rif resolution (suggested by Brian).
                // I think this fit perfect with the rif resolution.
                //double resl0 = 1;
                // It should be 1 to match the hight resolution rif, and I set it to 2 to test.
                double resl0 = 1;
                
                float const body_radius = std::min( scaff_radius, rif_radius );
                double const cart_grid = resl0*opt.hsearch_scale_factor/sqrt(3); // 1.5 is a big hack here.... 2 would be more "correct"
                double const hackysin = std::min( 1.0, resl0*opt.hsearch_scale_factor/2.0/ body_radius );
                runtime_assert( hackysin > 0.0 );
                double const rot_resl_deg0 = asin( hackysin ) * 180.0 / M_PI;
                int nside = std::ceil( opt.search_diameter / cart_grid );
                std::cout << "search dia.    : " << opt.search_diameter << std::endl;
                std::cout << "nside          : " << nside        << std::endl;
                std::cout << "resl0:           " << resl0 << std::endl;
                std::cout << "body_radius:     " << body_radius << std::endl;
                std::cout << "rif_radius:      " << rif_radius << std::endl;
                std::cout << "scaffold_radius: " << scaff_radius << std::endl;
                std::cout << "cart_grid:       " << cart_grid  << std::endl;
                std::cout << "rot_resl_deg0:   " << rot_resl_deg0 << std::endl;
                I3 nc( nside, nside, nside );
                F3 lb = F3( 0.0, 0.0, 0.0) + F3( -cart_grid*nside/2.0, -cart_grid*nside/2.0, -cart_grid*nside/2.0 );
                F3 ub = F3( 0.0, 0.0, 0.0) + F3(  cart_grid*nside/2.0,  cart_grid*nside/2.0,  cart_grid*nside/2.0 );
                std::cout << "cart grid ub " << ub << std::endl;
                std::cout << "cart grid lb " << lb << std::endl;
                std::cout << "(ub-lb/nc) = " << ((ub-lb)/nc.template cast<float>()) << std::endl;
                std::cout << "cartcen to corner (cart. covering radius): " << sqrt(3.0)*cart_grid/2.0 << std::endl;
                shared_ptr<RifDockNestDirector> nest_director = make_shared<RifDockNestDirector>( rot_resl_deg0, lb, ub, nc, 1 );
                std::cout << "NestDirector:" << std::endl << *nest_director << std::endl;
                std::cout << "nest size:    " << float(nest_director->size(0, RifDockIndex()).nest_index) << std::endl;
                std::cout << "size of search space: ~" << float(nest_director->size(0, RifDockIndex()).nest_index) << " grid points" << std::endl;
                
                shared_ptr<RifDockSeedingDirector> seeding_director = make_shared<RifDockSeedingDirector>(seeding_positions_p, 1, -1 /*maximum allowed ang*/);
                
                std::vector<DirectorBase> director_list;
                director_list.push_back( nest_director );  // Nest director must come first!!!!
                director_list.push_back( seeding_director );
                
                director = make_shared<RifDockDirector>(director_list);
                
                n_director = nest_director;
                
            }
            
            // parse and read in the exhausitive searching positions.
            std::vector< std::pair< int64_t, EigenXform > > xform_positions;
            {
                runtime_assert_msg(parse_exhausitive_searching_file(opt.xform_fname, xform_positions, 10), "Faild to parse the xform file!!!");
            }
            
            // nowo I have ererything (seeding positions and exhausitive searching positions ), I can just do the calcusition without the director!!
            
            std::vector< std::vector< SearchPointWithRots > > packed_results;
            std::vector< ScenePtr > scene_pt( omp_max_threads_1() );
            std::vector< ScenePtr > scene_hpack_pt( omp_max_threads_1() );
            {
                std::vector< std::vector< SearchPoint > > samples;
                
                int64_t seeding_size = director->size(0, RifDockIndex() ).seeding_index;
                int64_t nest_size = xform_positions.size();
                
                std::cout << std::endl;
                std::cout << "========================================== Now the Hsearch stage ============================================" << std::endl << std::endl;
                std::cout << "Rifine stage: begin threaded exhausitive searching, " << KMGT( nest_size ) << " samples, " << KMGT( seeding_size ) << " seeding positions." << std::endl;
                
                
                int64_t const out_interval = seeding_size * nest_size / 109;
                std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
                start = std::chrono::high_resolution_clock::now();
                
                // prefill the samples vector
                samples.resize( seeding_size );
                for ( int64_t ii = 0; ii < seeding_size; ++ii) samples[ii].resize( nest_size );
                
                
                // set up scene for each thread
                BOOST_FOREACH( ScenePtr & s, scene_pt ) s = scene_minimal->clone_deep();
                
                // do a real exhausitive search of all the searching points.
                #ifdef USE_OPENMP
                #pragma omp parallel for schedule(dynamic,64)
                #endif
                for ( int64_t i_seed = 0; i_seed < seeding_size; ++i_seed ) {
                    for ( int64_t i_samp = 0; i_samp < nest_size; ++i_samp ) {
                        
                        if ( ( i_seed * nest_size + i_samp ) % out_interval == 0 ) { std::cout << '*'; std::cout.flush(); }
                        
                        ScenePtr tscene( scene_pt[omp_get_thread_num()] );
                        
                        EigenXform p(EigenXform::Identity());
                        p.rotate( xform_positions[i_samp].second.rotation() * seeding_positions[i_seed].rotation() );
                        p.translation() = xform_positions[i_samp].second.translation() + seeding_positions[i_seed].translation();
                        
                        
                        tscene->set_position(1, p);
                        
                        samples[i_seed][i_samp].index.nest_index = xform_positions[i_samp].first;
                        samples[i_seed][i_samp].index.seeding_index = i_seed;
                        samples[i_seed][i_samp].score = objectives.back()->score(*tscene);
                        
                    } // end loop of each xform position
                } // end loop of each seeding position
                
                // selecting the best for hack pack.
                std::cout << "" << std::endl;
                std::cout << "Done with refinement stage, now selecting the best results for hpack." << std::endl;
                int64_t const len = nest_size * seeding_size;
                
                
                #ifdef USE_OPENMP
                #pragma omp parallel for schedule(dynamic,64)
                #endif
                for ( int64_t i_seed = 0; i_seed < seeding_size; ++i_seed ) {
                    if ( len > opt.beam_size /* opt.beam_size / opt.DIMPOW2 */ ) {
                        int64_t len_each = int( std::ceil ( opt.beam_size / seeding_size ) );
                        std::nth_element( samples[i_seed].begin(), samples[i_seed].begin()+len_each, samples[i_seed].end() );
                        samples[i_seed].resize( len_each );
                    }
                    std::sort(samples[i_seed].begin(), samples[i_seed].end());
                    
                }
                
                
                end = std::chrono::high_resolution_clock::now();
                std::cout << "Rifine and selection stage done, and it took "  << (end - start).count() << " CPU ticks." << std::endl << std::endl;
                
                ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                /////////////////////////////////                 HACK PACK                  ///////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                // TODO: add code for the case without packing. But I think that should never happen.
                
                if ( opt.hack_pack /* opt.hackpack, as I think I should always do hack pack. There is now way to skip this step. */ ) {
                    std::cout << "========================================= Now the HackPack stage ============================================" << std::endl << std::endl;
                    
                    // change the scene_minimal
                    {
                        BOOST_FOREACH( ScenePtr & s, scene_hpack_pt ) s = scene_full->clone_deep();
                    }
                    
                    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
                    start = std::chrono::high_resolution_clock::now();
                    
                    
                    
                    // As Brian said, the hackpack is very fast, so Maybe I should pack all of them?????????
                    // So here I changed the logic of pack samples selection.
                    // It would be a little slow to hackpack all of them.
                    // But finally, I found the hack pack is much slower, so I should pick the top x percent for the hpack stage.
                    packed_results.resize(seeding_size);
                    
                    // it is a little bit tricky to setup the progress  bar for the hackpack step, as the number of the samples in each seeding position is different.
                    // also here setup the number of samples to hackpack.
                    std::vector<std::vector< bool > > progress_bar(seeding_size);
                    {
                        int64_t out_interval = 0;
                        int64_t total_npack = 0;
                        for (int64_t i_seed = 0; i_seed < seeding_size; ++i_seed) {
                            int cutoff_score_pos = -1;
                            int each_seed_size = samples[i_seed].size();
                            for ( int64_t i_samp = 0; i_samp < each_seed_size; ++i_samp) {
                                // here I change this to cluster_score_cut, it is better as .......
                                // cluster score here, good .......
                                if ( samples[i_seed][i_samp].score >= opt.cluster_score_cut ) {
                                    cutoff_score_pos = i_samp;
                                    break;
                                }
                            }
                            // define the number of samples in each seeding position.
                            int n_sample = int( std::ceil ( ((-1 == cutoff_score_pos) ? each_seed_size: cutoff_score_pos) * opt.hack_pack_frac ) );
                            
                            packed_results[i_seed].resize(n_sample);
                            progress_bar[i_seed].resize( n_sample );
                            total_npack += n_sample;
                        }
                        out_interval = total_npack / 109;
                        out_interval = out_interval > 0 ? out_interval : 1;
                        total_npack = 0;
                        for (int64_t i_seed = 0; i_seed < seeding_size; ++i_seed) {
                            for (int64_t i_samp = 0; i_samp < packed_results[i_seed].size(); ++i_samp) {
                                total_npack += 1;
                                if (total_npack % out_interval == 0) {
                                    progress_bar[i_seed][i_samp] = true;
                                } else {
                                    progress_bar[i_seed][i_samp] = false;
                                }
                            }
                        }
                    }
                    
                    
                    // the real hackpack part.
                    
                    #ifdef USE_OPENMP
                    #pragma omp parallel for schedule(dynamic,64)
                    #endif
                    for ( int64_t i_seed = 0; i_seed < seeding_size; ++i_seed ) {
                        for (int64_t i_sample = 0; i_sample < packed_results[i_seed].size(); ++i_sample) {
                            
                            if ( true == progress_bar[i_seed][i_sample] ) { std::cout << '*'; std::cout.flush(); }
                            
                            packed_results[i_seed][i_sample].index = samples[i_seed][i_sample].index;
                            packed_results[i_seed][i_sample].prepack_rank = i_sample;
                            
                            // choose the scene.
                            ScenePtr tscene( scene_hpack_pt[omp_get_thread_num()] );
                            bool director_success = director->set_scene( packed_results[i_seed][i_sample].index, 0, *tscene );
                            if ( ! director_success ) {
                                // packed_results[ ipack ].rotamers(); // this initializes it to blank
                                packed_results[ i_seed ][ i_sample ].score = 9e9;
                                continue;
                            }
                            //packed_results[ i_seed ][ i_sample ].score = samples[i_seed][i_sample].score;
                            packed_results[ i_seed ][ i_sample ].score = packing_objective->score_with_rotamers( *tscene, packed_results[i_seed][ i_sample ].rotamers() );
                        }
                    }
                    
                    // now sorting the results.
                    #ifdef USE_OPENMP
                    #pragma omp parallel for schedule(dynamic,64)
                    #endif
                    for (int64_t i_seed=0; i_seed < seeding_size; ++i_seed) {
                        std::sort(packed_results[i_seed].begin(), packed_results[i_seed].end());
                    }
                    
                    
                    end = std::chrono::high_resolution_clock::now();
                    std::cout << std::endl;
                    std::cout << "Hackpack stage done, and it took "  << (end - start).count() << " CPU ticks." << std::endl << std::endl;
                } // end if of hack pack
                
            } // end block of the position refinement and hpack.
            
            
            
            // condense and select the good results, I would also use the redundancy filter to remove the bad ones.
            std::vector< std::vector< SearchPointWithRots > > scoremin_results;
            {
                print_header( "find good seeding pos by cluster score and filter redundant hackpack results" );
                int64_t const seeding_size = director->size(0, RifDockIndex()).seeding_index;
                
                // cluster score in each seeding position.
                std::vector< int > cluster_score(seeding_size, -1);
                for (int64_t i_seed = 0; i_seed < seeding_size; ++i_seed) {
                    int c_score = 0;
                    for (int64_t i_samp = 0; i_samp < packed_results[i_seed].size(); ++i_samp) {
                        if ( packed_results[i_seed][i_samp].score > opt.cluster_score_cut ) {
                            break;
                        }
                        c_score += 1;
                    }
                    cluster_score[i_seed] = c_score;
                }
                
                
                // print the cluster score of each seeding position.
                std::cout << "Cluster score of each seeding pos: " << std::endl;
                for (int64_t i_seed = 0; i_seed < seeding_size; ++i_seed) {
                    std::cout << i_seed << ":" << cluster_score[i_seed] << " ";
                }
                std::cout << std::endl;
                
                std::vector<int> cluster_score_nth;
                for ( auto a : cluster_score ) if(a > 0) cluster_score_nth.push_back(a);
                int c_score_cut = int( std::floor ( (1 - opt.keep_top_clusters_frac) * cluster_score_nth.size() ) );
                std::nth_element( cluster_score_nth.begin(), cluster_score_nth.begin()+c_score_cut, cluster_score_nth.end() );
                int cutoff_of_num = cluster_score_nth[c_score_cut];
                
                std::cout << std::endl << "Cutoff value of clusters: " << cutoff_of_num << std::endl;
                
                for (int i_seed = 0; i_seed < seeding_size; ++i_seed) {
                    if ( cluster_score[i_seed] >= cutoff_of_num && packed_results[i_seed].size() > 0 && packed_results[i_seed][0].score <= opt.global_score_cut ) {
                        scoremin_results.resize( scoremin_results.size() + 1 );
                        scoremin_results.back().push_back( packed_results[i_seed][0] );
                        for (int i_samp = 1; i_samp < packed_results[i_seed].size() && packed_results[i_seed][i_samp].score < opt.global_score_cut; ++i_samp) {
                            bool is_redundant = false;
                            for ( int64_t i_selected = 0; i_selected < scoremin_results.back().size(); ++i_selected ) {
                                director->set_scene( packed_results[i_seed][i_samp].index, 0, *scene_minimal );
                                director->set_scene( scoremin_results.back()[i_selected].index, 0, *scene_full );
                                EigenXform p1 = scene_minimal->position(1);
                                EigenXform p2 = scene_full->position(1);
                                float mag = devel::scheme::xform_magnitude( p2 * p1.inverse(), redundancy_filter_rg );
                                if ( mag <= opt.redundancy_filter_mag ) {
                                    is_redundant = true;
                                    break;
                                }
                            }
                            if ( !is_redundant ) {
                                scoremin_results.back().push_back( packed_results[i_seed][i_samp] );
                            }
                        }
                    }
                }
                
                // report the total number of results left
                int64_t total_left = 0;
                std::cout << std::endl << "Num of samples left in selected seeding position: " << std::endl;
                for ( auto & a : scoremin_results ){
                    int num = a.size();
                    total_left += num;
                    std::cout << a[0].index.seeding_index << ":" << num << " ";
                }
                std::cout << std::endl << std::endl;
                std::cout << "Number of total searching points left: " << total_left << std::endl << std::endl;
            }
            

            
            
            
            
            
            
            
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////       Rosetta score and min     /////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // rossetta score and min, you should always do rosetta score and min, there is no way to skip this step
            // TODO: add code if skip this step, should I do this. Just for fun and a waste of time??
            // Maybe here I can use the rosetta_score_each_seeding_at_least flag, as it is useful for each seeding position.
            std::vector<SearchPointWithRots> rifine_results;
            bool const do_rosetta_score = opt.rosetta_score_fraction > 0 || opt.rosetta_score_then_min_below_thresh > -9e8 || opt.rosetta_score_each_seeding_at_least > 0;
            if ( do_rosetta_score && opt.hack_pack ) {
                
                std::cout << "=========================================  Rosetta score and min ============================================" << std::endl << std::endl;
                std::chrono::time_point<std::chrono::high_resolution_clock> start_rosetta = std::chrono::high_resolution_clock::now();
                
               
                
                
                int do_min = 2;
                if ( opt.rosetta_min_fraction == 0.0 ) do_min = 1;
                
                std::vector<bool> is_scaffold_fixed_res(scaffold.size()+1,true);
                for(int designable : scaffold_res){
                    is_scaffold_fixed_res[designable] = false;
                }
                
                
                // I don't know how to name the final result. Maybe just rifine.
                int total_score_min = 0;
                int out_interval = 1;
                {
                    for (int64_t i_seed = 0; i_seed < scoremin_results.size(); ++i_seed) {
                        int n_sample = int( std::round( scoremin_results[i_seed].size() * opt.rosetta_score_fraction ) );
                        n_sample = n_sample > opt.rosetta_score_each_seeding_at_least ? n_sample : opt.rosetta_score_each_seeding_at_least;
                        n_sample = scoremin_results[i_seed].size() > n_sample ? n_sample : scoremin_results[i_seed].size();
                        total_score_min += n_sample;
                        for ( int64_t i_samp = 0; i_samp < n_sample; ++i_samp) {
                            rifine_results.push_back( scoremin_results[i_seed][i_samp] );
                        }
                    }
                }
                
                
                for (int minimizing = 0; minimizing < do_min; ++minimizing) {

                    
                    std::vector<core::kinematics::MoveMapOP> movemap_pt( omp_max_threads() );
                    std::vector<protocols::minimization_packing::MinMoverOP> minmover_pt( omp_max_threads() );
                    std::vector<core::scoring::ScoreFunctionOP> scorefunc_pt( omp_max_threads() );
                    std::vector<core::pose::Pose> work_pose_pt        ( omp_max_threads() );
                    std::vector<core::pose::Pose> both_per_thread     (omp_max_threads());
                    for ( int i = 0; i < omp_max_threads(); ++i ) {
                        // both_full_per_thread[i] = both_full_pose;
                        if( opt.replace_orig_scaffold_res ){
                            both_per_thread[i] = both_full_pose;
                        } else {
                            both_per_thread[i] = both_pose;
                        }
                        // create score function
                        // scorefunc_pt[i] = core::scoring::ScoreFucntionFactory::create_score_function( opt.rosetta_soft_score );
                        if ( minimizing ) {
                            if ( opt.rosetta_hard_min ) {
                                scorefunc_pt[i] = core::scoring::ScoreFunctionFactory::create_score_function( opt.rosetta_hard_score );
                            } else {
                                 scorefunc_pt[i] = core::scoring::ScoreFunctionFactory::create_score_function( opt.rosetta_soft_score );
                            }
                        } else if ( do_min == 2 ) {
                            // not minimizing, but will do minimization.
                            scorefunc_pt[i] = core::scoring::ScoreFunctionFactory::create_score_function( opt.rosetta_soft_score );
                            if( !opt.rosetta_hard_min ){
                                scorefunc_pt[i]->set_weight( core::scoring::fa_rep, scorefunc_pt[i]->get_weight(core::scoring::fa_rep)*0.7 );
                                scorefunc_pt[i]->set_weight( core::scoring::fa_dun, scorefunc_pt[i]->get_weight(core::scoring::fa_dun)*0.7 );
                            }
                        } else {
                            // not minimizing at all, score pass only.
                            scorefunc_pt[i] = core::scoring::ScoreFunctionFactory::create_score_function( opt.rosetta_soft_score );
                            if( !opt.rosetta_hard_min ){
                                scorefunc_pt[i]->set_weight( core::scoring::fa_rep, scorefunc_pt[i]->get_weight(core::scoring::fa_rep)*1.0 );
                                scorefunc_pt[i]->set_weight( core::scoring::fa_dun, scorefunc_pt[i]->get_weight(core::scoring::fa_dun)*1.0 );
                            }
                        }
                        if( target.size() == 1 ){
                            // assume this is a ligand, so hbonding is important
                            scorefunc_pt[i]->set_weight( core::scoring::fa_elec    , scorefunc_pt[i]->get_weight(core::scoring::fa_elec    )*2.0 );
                            scorefunc_pt[i]->set_weight( core::scoring::hbond_sc   , scorefunc_pt[i]->get_weight(core::scoring::hbond_sc   )*2.0 );
                            scorefunc_pt[i]->set_weight( core::scoring::hbond_bb_sc, scorefunc_pt[i]->get_weight(core::scoring::hbond_bb_sc)*2.0 );
                        } else {
                            scorefunc_pt[i]->set_weight( core::scoring::fa_elec    , scorefunc_pt[i]->get_weight(core::scoring::fa_elec    )*1.0 );
                            scorefunc_pt[i]->set_weight( core::scoring::hbond_sc   , scorefunc_pt[i]->get_weight(core::scoring::hbond_sc   )*1.0 );
                            scorefunc_pt[i]->set_weight( core::scoring::hbond_bb_sc, scorefunc_pt[i]->get_weight(core::scoring::hbond_bb_sc)*1.0 );
                        }
                        
                        // set up the minmover
                        movemap_pt[i] = core::kinematics::MoveMapOP( new core::kinematics::MoveMap() );
                        movemap_pt[i]->set_chi(true);
                        movemap_pt[i]->set_jump(true);
                        for(int ir = 1; ir <= both_full_pose.size(); ++ir){
                            bool is_scaffold = ir <= scaffold.size();
                            if( is_scaffold ) movemap_pt[i]->set_bb(ir, opt.rosetta_min_allbb || opt.rosetta_min_scaffoldbb );
                            else              movemap_pt[i]->set_bb(ir, opt.rosetta_min_allbb || opt.rosetta_min_targetbb );
                            if( opt.rosetta_min_fix_target && !is_scaffold ){
                                movemap_pt[i]->set_chi(ir,false);
                            }
                        }
                        minmover_pt[i] = protocols::minimization_packing::MinMoverOP(new protocols::minimization_packing::MinMover( movemap_pt[i], scorefunc_pt[i], "dfpmin_armijo_nonmonotone", 0.001, true ) );
                    } // end loop of initialize of scorefunctions.
                    
                    
                    // decide the numbers to score and min. I think here I am using the global score cut to remove the "bad" designs.
                    {
                        if ( minimizing ) {
                            
                            int num_to_min = int ( std::ceil ( total_score_min * opt.rosetta_min_fraction ) );
                            if ( num_to_min < opt.rosetta_min_at_least ) {
                                total_score_min = opt.rosetta_min_at_least > total_score_min ? total_score_min : opt.rosetta_min_at_least ;
                            } else {
                                total_score_min = num_to_min;
                            }
                            
                            out_interval = total_score_min / 109;
                            out_interval = out_interval > 0 ? out_interval : 1;
                            std::cout << std::endl << "rosetta min on: "   << KMGT(total_score_min) << std::endl;
                        } else {
                            // I just think 109 is beautiful, there is no reason why.
                            out_interval = total_score_min / 109;
                            out_interval = out_interval > 0 ? out_interval : 1;
                            std::cout << "rosetta score on: " << KMGT(total_score_min) << std::endl;
                        }
                    } // end block of the code to decide the numbers to score and min.
                    
                    
                    // the real score and min code.
                    #ifdef USE_OPENMP
                    #pragma omp parallel for schedule(dynamic,1)
                    #endif
                    for (int64_t i_samp = 0; i_samp < total_score_min; ++i_samp) {
                            
                            if ( i_samp % out_interval == 0 ) { std::cout << '*'; std::cout.flush();  };
                            
                            // get the current thread num
                            int const ithread = omp_get_thread_num();
                            
                            // get the transform of the scaffold
                            director->set_scene(rifine_results[i_samp].index, 0, *scene_pt[ithread]);
                            // I don't think there is any need to align everything to the scaffold.
                            EigenXform xposition1 = scene_pt[ithread]->position(1);
                            
                            // transform the scaffold to the position
                            core::pose::Pose & pose_to_min( work_pose_pt[ithread] );
                            pose_to_min = both_per_thread[ithread];
                            xform_pose( pose_to_min, eigen2xyz(xposition1) , 1 , scaffold.size() );
                            
                            // place the rotamers.
                            core::chemical::ResidueTypeSetCAP rts = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");
                            std::vector<bool> is_rif_res(pose_to_min.size(),false);
                            for ( int ipr = 0; ipr < rifine_results[i_samp].numrots(); ++ipr ) {
                                int ires = scaffres_l2g.at( rifine_results[i_samp].rotamers().at(ipr).first );
                                int irot =                  rifine_results[i_samp].rotamers().at(ipr).second;
                                core::conformation::ResidueOP newrsd = core::conformation::ResidueFactory::create_residue( rts.lock()->name_map(rot_index.resname(irot)) );
                                pose_to_min.replace_residue( ires+1, *newrsd, true );
                                is_rif_res[ires] = true;
                                for ( int ichi = 0; ichi < rot_index.nchi(irot); ++ichi ) {
                                    pose_to_min.set_chi( ichi+1, ires+1, rot_index.chi( irot, ichi ) );
                                }
                            }
                            
                            
                            EigenXform Xtorifframe = EigenXform::Identity();
                            
                            // replace the clash residues
                            std::vector<int> replaced_scaffold_res, rifres;
                            auto alaop = core::conformation::ResidueFactory::create_residue( rts.lock()->name_map("ALA" ) );
                            
                            std::vector<int> rifatypemap = get_rif_atype_map();
                            for (int ir = 1; ir <= scaffold.size(); ++ir) {
                                auto const & ires = pose_to_min.residue(ir);
                                if ( !ires.is_protein() ) continue;
                                if( ires.aa()==core::chemical::aa_gly ||
                                   ires.aa()==core::chemical::aa_ala ||
                                   ires.aa()==core::chemical::aa_pro ) continue;
                                if(is_rif_res[ir-1]){
                                    rifres.push_back(ir);
                                    continue;
                                }
                                if(is_scaffold_fixed_res[ir]) continue;
                                
                                bool ir_clash = false;
                                float evtarget = 0.0;
                                for (int ia = 6; ia <= ires.nheavyatoms(); ++ia) {
                                    auto const & ixyz = ires.xyz(ia);
                                    
                                    Eigen::Vector3f satm;
                                    for(int k = 0; k < 3; ++k) satm[k] = ixyz[k];
                                    int const irifatype = rifatypemap[ires.atom_type_index(ia)];
                                    evtarget += target_field_by_atype.at(irifatype)->at(satm);
                                }
                                if( evtarget > 3.0f ) ir_clash = true;
                                
                                // check against other rif res
                                for( int ipr = 0; ipr < rifine_results[i_samp].numrots(); ++ipr ){
                                    int jr = 1+scaffres_l2g.at( rifine_results[i_samp].rotamers().at(ipr).first );
                                    auto const & jres = pose_to_min.residue(jr);
                                    // should do rsd nbr check... but speed not critical here ATM...
                                    for( int ia = 6; ia <= ires.nheavyatoms(); ++ia){
                                        auto const & ixyz = ires.xyz(ia);
                                        for( int ja = 6; ja <= jres.nheavyatoms(); ++ja){
                                            auto const & jxyz = jres.xyz(ja);
                                            if( ixyz.distance_squared(jxyz) < 9.0 ){
                                                ir_clash = true;
                                            }
                                        }
                                    }
                                }
                                if ( ir_clash ) {
                                    pose_to_min.replace_residue(ir, *alaop, true);
                                    replaced_scaffold_res.push_back(ir);
                                }
                            } // check clash and replace to AlA
                            
                            
                            // the minimizing block
                            if ( minimizing ) {
                                minmover_pt[ithread]->apply( pose_to_min );
                            } else {
                                scorefunc_pt[ithread]->score( pose_to_min );
                            }// the minimizing block
                            
                            // store the score
                            if ( opt.rosetta_score_total ) {
                                rifine_results[i_samp].score = pose_to_min.energies().total_energy();
                            } else {
                                double rosetta_score = 0.0;
                                auto const & weights = pose_to_min.energies().weights();
                                if( !opt.rosetta_score_ddg_only ){
                                    for( int ir = 1; ir <= scaffold.size(); ++ir ){
                                        if( is_rif_res[ir-1] ){
                                            rosetta_score += pose_to_min.energies().onebody_energies(ir).dot(weights);
                                        }
                                    }
                                }
                                if( !opt.rosetta_score_ddg_only && target.size()==1 ){
                                    // is ligand, add it's internal energy
                                    rosetta_score += pose_to_min.energies().onebody_energies(pose_to_min.size()).dot(weights);
                                }
                                auto const & egraph = pose_to_min.energies().energy_graph();
                                for(int ir = 1; ir <= egraph.num_nodes(); ++ir){
                                    for ( utility::graph::Graph::EdgeListConstIter
                                         iru  = egraph.get_node(ir)->const_upper_edge_list_begin(),
                                         irue = egraph.get_node(ir)->const_upper_edge_list_end();
                                         iru != irue; ++iru
                                         ){
                                        EnergyEdge const & edge( static_cast< EnergyEdge const & > (**iru) );
                                        int jr = edge.get_second_node_ind();
                                        
                                        // this is DDG
                                        if( ir <= scaffold.size() && jr > scaffold.size() ){
                                            // ir in scaff, jr in target
                                            rosetta_score += edge.dot(weights);
                                        }
                                        if( !opt.rosetta_score_ddg_only && jr <= scaffold.size() ){
                                            // ir & jr in scaffold
                                            if( is_rif_res[ir-1] || is_rif_res[jr-1] ){
                                                double const edgescore = edge.dot(weights);
                                                if( edgescore > 0.0 ){
                                                    // always assess full score for bad interactions
                                                    rosetta_score += edgescore;
                                                } else if( is_rif_res[ir-1] && is_rif_res[jr-1] ){
                                                    // both rif residues
                                                    rosetta_score += opt.rosetta_score_rifres_rifres_weight * edgescore;
                                                    // bonus for hbonds between rif residues
                                                    rosetta_score += edge[core::scoring::hbond_sc];
                                                } else {
                                                    // rest: one rif res, one other scaff res
                                                    rosetta_score += opt.rosetta_score_rifres_scaffold_weight * edgescore;
                                                }
                                            } else {
                                                // scaffold / scaffold ignored
                                            }
                                        }
                                    }
                                }
                                rifine_results[i_samp].score = rosetta_score;
                            }
                            
                            // store the results?
                            if( (minimizing+1 == do_min)     && rifine_results[i_samp].score < opt.rosetta_score_cut ){
                                rifine_results[i_samp].pose_ = core::pose::PoseOP( new core::pose::Pose(pose_to_min) );
                                for(int ir : rifres){
                                    rifine_results[i_samp].pose_->pdb_info()->add_reslabel(ir, "RIFRES" );
                                }
                                for(int ir : replaced_scaffold_res){
                                    rifine_results[i_samp].pose_->pdb_info()->add_reslabel(ir, "PRUNED" );
                                }
                            }
                            
                        } // end loop of each sample in each seeding position.
                
                
                    // sort the results
                    __gnu_parallel::sort( rifine_results.begin(), rifine_results.begin() + total_score_min );
                 // end loop of minimizing
                }
                std::chrono::duration<double> elapsed_seconds_rosetta = std::chrono::high_resolution_clock::now()-start_rosetta;
                std::cout << std::endl << "Rosetta score and min done! Total " << elapsed_seconds_rosetta.count() << " CPU ticks spend." << std::endl << std::endl;
            }
        
            for (int64_t i_samp = 0; i_samp < rifine_results.size(); ++i_samp) {
                if ( i_samp > 20 || rifine_results[i_samp].pose_ == nullptr ) {
                    break;
                }
                std::cout << "SeedP: " << rifine_results[i_samp].index.seeding_index << "_" << rifine_results[i_samp].index.nest_index << " " << rifine_results[i_samp].score << std::endl;
                    rifine_results[i_samp].pose_->dump_pdb( "seed_" + str(rifine_results[i_samp].index.seeding_index) + "_" + str(rifine_results[i_samp].index.nest_index) + ".pdb" );
            } // end if of the rosetta score and min
            
            
            
            
            /*
             
             
             {
             
             // Here I just should output some of the results to make sure everything is fine.
             
             director->set_scene( packed_results[0][0].index, 0, *scene_minimal );
             director->set_scene( packed_results[0][1].index, 0, *scene_full );
             EigenXform p1 = scene_minimal->position(1);
             EigenXform p2 = scene_full->position(1);
             float mag = xform_magnitude( p2 * p1.inverse(),  redundancy_filter_rg);
             
             
             // It seem it perfectly recap
             SearchPointWithRots selected_result;
             
             
             //for ( int64_t i_samp = 0; i_samp <; ++i_samp) {
             int current_rank;
             for (int64_t i_samp = 0; i_samp < packed_results[0].size(); ++i_samp) {
             if (packed_results[0][i_samp].index.nest_index == 2370540 ) {
             selected_result = packed_results[0][i_samp];
             current_rank = i_samp;
             break;
             }
             }
             std::cout << "The xform_magnitute of p1 and p2: " << mag << std::endl;
             std::cout << "Current rank of the right solution: " << current_rank << std::endl;
             std::cout << "Prerank of the right solution: " << selected_result.prepack_rank << std::endl;
             std::cout << "Score of the right solution: " << selected_result.score << std::endl;
             
             
             //selected_result = packed_results[2][200];
             
             for (int i = 0; i < 2; ++i) {
             selected_result = packed_results[0][i];
             core::pose::Pose dump_pose = both_full_pose;
             director->set_scene( selected_result.index, 0, *scene_full );
             EigenXform p = scene_full->position(1);
             xform_pose( dump_pose, eigen2xyz(p), 1, scaffold.size() );
             
             
             core::chemical::ResidueTypeSetCAP rts = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");
             for ( int ipr = 0; ipr < selected_result.numrots(); ++ipr ) {
             int ires = scaffres_l2g.at( selected_result.rotamers().at(ipr).first );
             int irot =                  selected_result.rotamers().at(ipr).second;
             core::conformation::ResidueOP newrsd = core::conformation::ResidueFactory::create_residue( rts.lock()->name_map(rot_index.resname(irot)) );
             dump_pose.replace_residue( ires+1, *newrsd, true );
             for ( int ichi = 0; ichi < rot_index.nchi(irot); ++ichi ) {
             dump_pose.set_chi( ichi+1, ires+1, rot_index.chi( irot, ichi ) );
             }
             }
             //std::cout << "i_samp: " << str(i_samp) << " score: " << selected_result.score  << std::endl;
             dump_pose.dump_pdb( "rank_" + str(i) +  ".pdb");
             
             
             //scaffold_onebody_glob0
             std::cout << "The one body energy of " << str(i) << std::endl;
             for ( int ipr = 0; ipr < selected_result.numrots(); ++ipr ) {
             int ires = scaffres_l2g.at( selected_result.rotamers().at(ipr).first );
             int irot =                  selected_result.rotamers().at(ipr).second;
             std::cout << "res: " << str(ires+1) << " " << rot_index.resname(irot) << " irot "  << irot << " " << scaffold_onebody_glob0[ires][irot] << std::endl;
             }
             
             // scaffold_twobody
             std::cout << "Print out the two body table and see the energy" << std::endl;
             for ( int ipr_1 = 0; ipr_1 < selected_result.numrots() - 1; ++ipr_1) {
             for (int ipr_2 = 1; ipr_2 < selected_result.numrots(); ++ipr_2) {
             int ires_1 = scaffres_l2g.at( selected_result.rotamers().at(ipr_1).first );
             int irot_1 =                  selected_result.rotamers().at(ipr_1).second;
             
             int ires_2 = scaffres_l2g.at( selected_result.rotamers().at(ipr_2).first );
             int irot_2 =                  selected_result.rotamers().at(ipr_2).second;
             
             
             std::cout << "res1 " << str(ires_1 + 1) << " " << rot_index.resname(irot_1) << "res2 " << str(ires_2 + 1) << " " << rot_index.resname(irot_2) << " " << scaffold_twobody->twobody(ires_1, ires_2, irot_1, irot_2) << std::endl;
             }
             
             }
             
             // }
             
             }
             
             for (int i = 0; i < 20; ++i) {
             std::cout << "i_samp: " << str(i) << " score: " << packed_results[0][i].score << " prepack_rank: " << packed_results[0][i].prepack_rank << std::endl;
             }
             
             }
             
             
             
             
             
             {
             // Here I just should output some of the results to make sure everything is fine.
             // It seem it perfectly recap
             SearchPointWithRots selected_result;
             
             for ( int64_t i_samp = 0; i_samp < 20; ++i_samp) {
             
             
             selected_result = packed_results[2][i_samp];
             
             core::pose::Pose dump_pose = both_full_pose;
             director->set_scene( selected_result.index, 0, *scene_full );
             EigenXform p = scene_full->position(1);
             xform_pose( dump_pose, eigen2xyz(p), 1, scaffold.size() );
             
             
             core::chemical::ResidueTypeSetCAP rts = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");
             for ( int ipr = 0; ipr < selected_result.numrots(); ++ipr ) {
             int ires = scaffres_l2g.at( selected_result.rotamers().at(ipr).first );
             int irot =                  selected_result.rotamers().at(ipr).second;
             core::conformation::ResidueOP newrsd = core::conformation::ResidueFactory::create_residue( rts.lock()->name_map(rot_index.resname(irot)) );
             dump_pose.replace_residue( ires+1, *newrsd, true );
             for ( int ichi = 0; ichi < rot_index.nchi(irot); ++ichi ) {
             dump_pose.set_chi( ichi+1, ires+1, rot_index.chi( irot, ichi ) );
             }
             }
             std::cout << "SeedP: " << 2 << " " << packed_results[2][i_samp].score << " " << packed_results[2][i_samp].prepack_rank << " " << i_samp << std::endl;
             //dump_pose.dump_pdb( "Seeding_" + str(i_samp) + "_best.pdb");
             }
             }
             
             // check the outputs of hpack
             {
             // Here I just should output some of the results to make sure everything is fine.
             // It seem it perfectly recap
             SearchPointWithRots selected_result;
             
             for ( int64_t i_samp = 0; i_samp < 20; ++i_samp) {
             
             
             selected_result = packed_results[2][i_samp];
             
             core::pose::Pose dump_pose = both_full_pose;
             director->set_scene( selected_result.index, 0, *scene_full );
             EigenXform p = scene_full->position(1);
             xform_pose( dump_pose, eigen2xyz(p), 1, scaffold.size() );
             
             
             core::chemical::ResidueTypeSetCAP rts = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");
             for ( int ipr = 0; ipr < selected_result.numrots(); ++ipr ) {
             int ires = scaffres_l2g.at( selected_result.rotamers().at(ipr).first );
             int irot =                  selected_result.rotamers().at(ipr).second;
             core::conformation::ResidueOP newrsd = core::conformation::ResidueFactory::create_residue( rts.lock()->name_map(rot_index.resname(irot)) );
             dump_pose.replace_residue( ires+1, *newrsd, true );
             for ( int ichi = 0; ichi < rot_index.nchi(irot); ++ichi ) {
             dump_pose.set_chi( ichi+1, ires+1, rot_index.chi( irot, ichi ) );
             }
             }
             std::cout << "i_samp: " << str(i_samp) << " score: " << selected_result.score  << std::endl;
             dump_pose.dump_pdb( "Seeding_" + str(i_samp) + "_best.pdb");
             }
             
             for ( int64_t i_samp1 = 0; i_samp1 < 20; ++i_samp1) {
             for ( int64_t i_samp2 = i_samp1 + 1; i_samp2 < 20; ++i_samp2) {
             director->set_scene( packed_results[2][i_samp1].index, 0, *scene_minimal );
             director->set_scene( packed_results[2][i_samp2].index, 0, *scene_full );
             EigenXform p1 = scene_minimal->position(1);
             EigenXform p2 = scene_full->position(1);
             
             float xmag =  xform_magnitude( p2 * p1.inverse(), redundancy_filter_rg );
             
             std::cout << xmag << "  ";
             }
             std::cout << std::endl;
             }
             
             }
             
             for (int64_t i_seed = 0; i_seed < seeding_size; ++i_seed) {
                if (packed_results[i_seed].size() > 0) {
                    std::cout << "SeedP: " << i_seed << " " << packed_results[i_seed][0].score << " " << packed_results[i_seed][0].prepack_rank << std::endl;
                }
             }
            
            // Here I just should output some of the results to make sure everything is fine.
            // It seem it perfectly recap
            {
                SearchPointWithRots selected_result;
                
                for (int64_t i_samp = 0; i_samp < packed_results[2].size(); ++i_samp) {
                    if (0 == packed_results[2][i_samp].prepack_rank) {
                        selected_result = packed_results[2][i_samp];
                    }
                }
                
                core::pose::Pose dump_pose = both_full_pose;
                director->set_scene( selected_result.index, 0, *scene_full );
                EigenXform p = scene_full->position(1);
                xform_pose( dump_pose, eigen2xyz(p), 1, scaffold.size() );
                
                core::chemical::ResidueTypeSetCAP rts = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");
                for ( int ipr = 0; ipr < selected_result.numrots(); ++ipr ) {
                    int ires = scaffres_l2g.at( selected_result.rotamers().at(ipr).first );
                    int irot =                  selected_result.rotamers().at(ipr).second;
                    core::conformation::ResidueOP newrsd = core::conformation::ResidueFactory::create_residue( rts.lock()->name_map(rot_index.resname(irot)) );
                    dump_pose.replace_residue( ires+1, *newrsd, true );
                    for ( int ichi = 0; ichi < rot_index.nchi(irot); ++ichi ) {
                        dump_pose.set_chi( ichi+1, ires+1, rot_index.chi( irot, ichi ) );
                    }
                }
                
                dump_pose.dump_pdb("best.pdb");
            }
                
            
            
            for (int64_t i_seed = 0; i_seed < seeding_size; ++i_seed) {
                std::cout << "SeedP: " << i_seed << " " << samples[i_seed][0].score << " " << samples[i_seed][10].score << " " << samples[i_seed][20].score << std::endl;
            }
             
            RifDockIndex rdi;
            for (int i = 0; i < nest_size; ++i) {
                rdi.seeding_index = 1;
                rdi.nest_index = i;
                
                bool director_success = n_director->set_scene( rdi, 0, *scene_minimal );
                
                if ( director_success ) {
                    EigenXform p = scene_minimal->position(1);
                    double ang = Eigen::AngleAxisf( p.rotation() ).angle();
                    if ( std::abs( ang ) <= 20 ) {
                        std::cout << "SP" << " " << i << " " << ang << " "
                                  << p.linear().row(0) << " " << p.linear().row(1)<< " " << p.linear().row(2) << " "
                                  << p.translation().x() << " " << p.translation().y() << " " << p.translation().z() << std::endl;
                    }
                }
            }
            
            
            
            // The code here seems augly, that's because I must be in compatible with the current current rifdock master.
            // the real searching, packing and scoring happen here.
            bool serach_failed = false;
            int64_t seeding_size = director->size(0, RifDockIndex() ).seeding_index;
            int64_t nest_size = director->size(0, RifDockIndex()).nest_index;
            std::vector< std::vector< SearchPoint > > samples( seeding_size );
            
            
            
            std::cout << "Seeding size: " << seeding_size << std::endl;
            std::cout << "Nest size" << nest_size << std::endl;
            
            {
                // do an exausitive searching.
                RifDockIndex rdi;
                for ( int64_t i_seed = 0; i_seed < seeding_size; ++i_seed ) {
                    rdi.seeding_index = i_seed;
                    // for each sampling position
                    for ( int64_t i_samp = 0; i_samp < nest_size; ++i_samp ) {
                        rdi.nest_index = i_samp;
                        
                        
                        bool director_success = director->set_scene( rdi, 0, *scene_minimal );
                        if ( ! director_success ) {
                            continue;
                        }
                        
                        
                        SearchPoint sp;
                        sp.index = rdi;
                        // sp.score = objectives.back()->score(*scene_minimal);
                        // samples[i_seed].push_back( sp );
                    }
                }
                for ( int64_t i_seed = 0; i_seed < seeding_size; ++i_seed ) {
                    std::cout << samples[i_seed].size() << std::endl;
                }
            }
            
            
            
            RifDockIndex rdi;
            for (int i = 0; i < nest_size; ++i) {
                rdi.seeding_index = 1;
                rdi.nest_index = i;
                std::cout << "I am here scoring the first time for " << i << std::endl;
                bool director_success = director->set_scene( rdi, 0, *scene_minimal );
                if ( director_success ) {
                    std::cout << "Score another three times: "<< i << std::endl;
                    std::cout << director->set_scene( rdi, 0, *scene_minimal ) << std::endl;
                    std::cout << director->set_scene( rdi, 0, *scene_minimal ) << std::endl;
                    std::cout << director->set_scene( rdi, 0, *scene_minimal ) << std::endl;
                }
            }
            
            
            
            if ( false )
            {
                core::pose::Pose final_pose = scaffold_centered;
                EigenXform final_pos;
                double final_score = 9e9;
                int64_t lowsc_index = 0;
                double current_sc;
                RifDockIndex rdi;
                for ( int64_t isamp = 0; isamp < nest_size0; ++isamp ) {
                    rdi.nest_index = isamp;
                    bool director_success = director->set_scene( rdi, 0, *scene_minimal );
                    if ( ! director_success ) {
                        continue;
                    }
                    
                    current_sc = objectives[3]->score(*scene_minimal);
                    if ( current_sc < final_score ) {
                        final_score = current_sc;
                        final_pos = scene_minimal->position(1);
                        lowsc_index = isamp;
                        
                    }
                    
                    if ( uniform(rng) < 0.01 ) {
                        std::string dump_name = str(isamp) + "_" + str(objectives[3]->score(*scene_minimal)) + ".pdb";
                        core::pose::Pose pose_dump( scaffold_centered );
                        xform_pose( pose_dump, eigen2xyz( scene_minimal->position(1) ) );
                        pose_dump.dump_pdb(dump_name);
                    }
                }
                xform_pose( final_pose, eigen2xyz( final_pos ) );
                final_pose.dump_pdb( "LowSc_" + str(lowsc_index) + "_" + str(final_score) + ".pdb" );
            }
            */
            
            
            
            
            
            
            
            
            
            
            
        } catch ( std::exception const & ex ) {
            std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
            std::cout << "error (below) on scaffold " << scaff_fname << " (will continue with others, if any)" << std::endl;
            std::cout << ex.what() << std::endl;
            std::cout << "scene residue numering (may help debug):" << std::endl;
            for( int i = 1; i <= scaffold_res.size(); ++i )
                std::cout << "scene res numbering: " << i-1 << " " << scaffold_sequence_glob0.at(scaffold_res[i]-1) << " pose number: " << scaffold_res[i] << std::endl;
        } catch ( ... ) {
            std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
            std::cout << "unknown error on scaffold " << scaff_fname << ", will continue with others, if any." << std::endl;
            std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
        }
        
    } // end of scaffold looping
    
    dokout.close();
    
    #ifdef USE_OPENMP
        omp_destroy_lock( &cout_lock );
        omp_destroy_lock( &dump_lock );
    #endif
    
    std::cout << "rifdock_DONE" << std::endl;
    
    return 0;
}
























