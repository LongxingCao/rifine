// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:

// INC

#define GLOBAL_VARIABLES_ARE_BAD
#include <rif_dock_test.hh>
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
	
	#include <riflib/rifdock_subroutines/HSearchFactory.hh>	

	#include <riflib/rifdock_subroutines/hack_pack.hh>
	#include <riflib/rifdock_subroutines/rosetta_rescore.hh>
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

		std::string const rif_type = get_rif_type_from_file( opt.rif_files.back() );
		BOOST_FOREACH( std::string fn, opt.rif_files ){
			std::string rif_type2 = get_rif_type_from_file( fn );
			runtime_assert_msg( rif_type==rif_type2, "mismatched rif types, expect: " + rif_type + " got: " + rif_type2 + " for " + fn );
		}
		std::cout << "read RIF type: " << rif_type << std::endl;

		cout << "Search Resls: " << opt.resl0;
			std::vector<float> RESLS(1,opt.resl0);
			for( int i = 1; i < opt.rif_files.size(); ++i ){
				RESLS.push_back( RESLS.back()/2.0 );
				cout << " " << RESLS.back();
			}
			cout << endl;

		std::cout << "opt.rosetta_score_fraction: " << opt.rosetta_score_fraction << std::endl;
		std::cout << "opt.rosetta_score_then_min_below_thresh: " << opt.rosetta_score_then_min_below_thresh << std::endl;
		std::cout << "opt.rosetta_score_at_least: " << opt.rosetta_score_at_least << std::endl;
		std::cout << "opt.rosetta_score_at_most: " << opt.rosetta_score_at_most << std::endl;
		std::cout << "opt.rosetta_min_fraction: " << opt.rosetta_min_fraction << std::endl;
		std::cout << "opt.rosetta_min_targetbb: " << opt.rosetta_min_targetbb << std::endl;
		std::cout << "opt.rosetta_min_allbb: " << opt.rosetta_min_allbb << std::endl;
		std::cout << "opt.rosetta_score_cut: " << opt.rosetta_score_cut << std::endl;
		std::cout << "opt.require_satisfaction: " << opt.require_satisfaction << std::endl;

		std::cout << "//////////////////////////// end options /////////////////////////////////" << std::endl;



		// for( int iscaff = 0; iscaff < opt.scaffold_fnames.size(); ++iscaff )
		// {
		// 	std::string scaff_fname = opt.scaffold_fnames.at(iscaff);
		// 	std::cout << scaff_fname << std::endl;
		// 	core::pose::Pose scaffold;
		// 	utility::vector1<core::Size> scaffold_res;
		// 	core::import_pose::pose_from_file(scaffold, scaff_fname);
		// 	scaff_fname = utility::file::file_basename(utility::file_basename(scaff_fname));
		// 	scaffold.dump_pdb(scaff_fname+"_0.pdb");
		// 	scaffold_res = devel::scheme::get_designable_positions_best_guess( scaffold, opt.dont_use_scaffold_loops );
		// 	::devel::scheme::pose_to_ala( scaffold, scaffold_res );
		// 	scaffold.dump_pdb(scaff_fname+"_1.pdb");
		// }
		// utility_exit_with_message("test_scaff sel");


		////////////////////////////// should be no more use of options at this point! ///////////////////////////


		double time_rif=0, time_pck=0, time_ros=0;

		std::mt19937 rng( 0);//std::random_device{}() );


		{
			std::string dokfile_fname_orig = opt.dokfile_fname;
			int i = 2;
			while( utility::file::file_exists(opt.dokfile_fname) ){
				opt.dokfile_fname = dokfile_fname_orig + "." + str(i);
				++i;
			}
			if( i != 2)
				std::cout << "WARNING!" << dokfile_fname_orig << " already exists, using "
			              << opt.dokfile_fname << " instead!" << std::endl;
             else
             	std::cout << "output scores to " << opt.dokfile_fname << std::endl;
		}
		utility::io::ozstream dokout( opt.dokfile_fname );


		devel::scheme::RifFactoryConfig rif_factory_config;
		rif_factory_config.rif_type = rif_type;
		shared_ptr<RifFactory> rif_factory = ::devel::scheme::create_rif_factory( rif_factory_config );


		// shared_ptr<RifFactory> rif_factory;
		// if ( opt.scaffold_provider_type == "SingleFile" ) {
		// 	rif_factory = ::devel::scheme::create_rif_factory<SingleFileScaffoldProvider>(rif_factory_config);
		// } else if ( opt.scaffold_provider_type == "Alex" ) {
		// 	rif_factory = ::devel::scheme::create_rif_factory<MorphingScaffoldProvider>(rif_factory_config);
		// } else {
		// 	utility_exit_with_message( "rif_dock_test: unknown scaffold provider type "+opt.scaffold_provider_type );
		// }




	print_header( "create rotamer index" );
		
		std::cout << "Loading " << opt.rot_spec_fname << "..." << std::endl;
		std::string rot_index_spec_file = opt.rot_spec_fname;
		shared_ptr< RotamerIndex > rot_index_p = ::devel::scheme::get_rotamer_index( rot_index_spec_file );
		RotamerIndex & rot_index( *rot_index_p );


		// {
		// 	utility::io::ozstream out("test.rotidx.gz",std::ios_base::binary);
		// 	rot_index.save(out);
		// 	out.close();
		// }
		// {
		// 	RotamerIndex ri2;
		// 	utility::io::izstream in("test.rotidx.gz",std::ios_base::binary);
		// 	ri2.load(in);
		// 	in.close();

		// 	std::cout << rot_index << std::endl;
		// 	std::cout << std::endl;
		// 	std::cout << ri2 << std::endl;

		// 	runtime_assert( ri2 == rot_index );
		// 	utility_exit_with_message("test rot index load/save");
		// }

		// std::cout << "================ RotamerIndex ===================" << std::endl;
		// std::cout << rot_index.size() << " " << rot_index.n_primary_rotamers() << std::endl;
		// std::cout << rot_index << std::endl;
		// {
		// 	utility::io::ozstream out("rot_index.pdb");
		// 	rot_index.dump_pdb( out );
		// 	utility_exit_with_message("ortsdn");
		// }
		// std::cout << "=================================================" << std::endl;

		RotamerRFOpts rotrfopts;
		rotrfopts.oversample     = opt.rotrf_oversample;
		rotrfopts.field_resl     = opt.rotrf_resl;
		rotrfopts.field_spread   = opt.rotrf_spread;
		rotrfopts.data_dir       = opt.rotrf_cache_dir;
		rotrfopts.scale_atr      = opt.rotrf_scale_atr;
		::devel::scheme::RotamerRFTablesManager rotrf_table_manager( rot_index_p, rotrfopts );
		// rotrf_table_manager.preinit_all();

		MakeTwobodyOpts make2bopts;
		// hacked by brian             VVVV
		make2bopts.onebody_threshold = 30.0;
		make2bopts.distance_cut = 15.0;
		make2bopts.hbond_weight = packopts.hbond_weight;
		make2bopts.favorable_2body_multiplier = opt.favorable_2body_multiplier;




	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	print_header( "read and prepare target structure" ); //////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	core::pose::Pose target;
	std::vector<SimpleAtom> target_simple_atoms;
	utility::vector1<core::Size> target_res;
	std::vector<HBondRay> target_donors, target_acceptors;
	float rif_radius=0.0, target_redundancy_filter_rg=0.0;
	{
		core::import_pose::pose_from_file( target, opt.target_pdb );

		if( opt.use_scaffold_bounding_grids ){
			for( int ir = 1; ir <= target.size(); ++ir ){
				utility::vector1<core::Size> resids(1,ir); // 1-index numbering
				std::vector<SchemeAtom> atoms;
				devel::scheme::get_scheme_atoms( target, resids, atoms );
				int restype = rot_index.chem_index_.resname2num( target.residue(ir).name3() ); // for UNK will be -1
				for( int ia = 0; ia < atoms.size(); ++ia){
					SchemeAtom const & a( atoms[ia] );
					runtime_assert( a.type() > 0 );
					if( a.type() >= 21 ) continue;
					SimpleAtom sa( a.position(), a.type(), restype, ia );
					target_simple_atoms.push_back(sa);
				}
			}
			std::cout << "target_simple_atoms.size() " << target_simple_atoms.size() << std::endl;
		}
		target_res = devel::scheme::get_res( opt.target_res_fname , target, /*nocgp*/false );
		get_rg_radius( target, target_redundancy_filter_rg, rif_radius, target_res, true ); // allatom for target
		rif_radius += 7.0; // hacky guess
		::devel::scheme::HBRayOpts hbopt;


		hbopt.satisfied_atoms = ::devel::scheme::get_satisfied_atoms(target);

		// utility_exit_with_message("MAKE SURE LKBALL STUFF ISN'T FUCKING UP!!!");


		BOOST_FOREACH( core::Size ir, target_res ){
			::devel::scheme::get_donor_rays   ( target, ir, hbopt, target_donors );
			::devel::scheme::get_acceptor_rays( target, ir, hbopt, target_acceptors );
		}
		std::cout << "target_donors.size() " << target_donors.size() << " target_acceptors.size() " << target_acceptors.size() << std::endl;
		// {
		// 	{
		// 		std::vector<HBondRay> tmpdon, tmpacc, tmpacclk;
		// 		::devel::scheme::HBRayOpts hbopt;
		// 		hbopt.lkball = false;
		// 		hbopt.withbb = true;
		// 		hbopt.add_acceptor_mid = true;
		// 		hbopt.satisfied_atoms = ::devel::scheme::get_satisfied_atoms(target);
		// 		BOOST_FOREACH( core::Size ir, target_res ){
		// 			::devel::scheme::get_donor_rays   ( target, ir, hbopt, tmpdon );
		// 			::devel::scheme::get_acceptor_rays( target, ir, hbopt, tmpacc );
		// 		}
		// 		hbopt.lkball = true;
		// 		BOOST_FOREACH( core::Size ir, target_res ){
		// 			::devel::scheme::get_acceptor_rays( target, ir, hbopt, tmpacclk );
		// 		}
		// 		target.dump_pdb("target.pdb");
		// 		utility::io::ozstream donout(utility::file_basename(opt.target_pdb)+"_donors.pdb");
		// 		::devel::scheme::dump_hbond_rays( donout, tmpdon, true );
		// 		donout.close();
		// 		utility::io::ozstream accout(utility::file_basename(opt.target_pdb)+"_orb_acceptors.pdb");
		// 		::devel::scheme::dump_hbond_rays( accout, tmpacc, false );
		// 		accout.close();
		// 		utility::io::ozstream accoutlk(utility::file_basename(opt.target_pdb)+"_lkb_acceptors.pdb");
		// 		::devel::scheme::dump_hbond_rays( accoutlk, tmpacclk, false );
		// 		accoutlk.close();
		// 		utility_exit_with_message("testing lkball replace orbs");
		// 	}
		// }
	}
	std::vector< VoxelArrayPtr > target_field_by_atype;
	std::vector< std::vector< VoxelArrayPtr > > target_bounding_by_atype;
	{
		target_bounding_by_atype.resize( RESLS.size() );
		devel::scheme::RosettaFieldOptions rfopts;
		rfopts.field_resl = opt.target_rf_resl;
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


		if( true ){
			// std::cout << "using target bounding grids, generating (or loading) them" << std::endl;
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
			if( opt.downscale_atr_by_hierarchy ){
				std::cout << "downscale_atr_by_hierarchy on target bounding steric grids" << std::endl;
				// std::cout << "  zeroing atr component of target bounding steric grids" << std::endl;
				for( int iresl = 0; iresl < RESLS.size(); ++iresl ){
					float correction = 1.0/RESLS[iresl];
					if( correction >= 1.0 ) break;
					BOOST_FOREACH( VoxelArrayPtr vap, target_bounding_by_atype[iresl] ){
						if( vap != nullptr ){
							std::exception_ptr exception = nullptr;
							#ifdef USE_OPENMP
							#pragma omp parallel for schedule(dynamic,64)
							#endif
							for( int k = 0; k < vap->num_elements(); ++k ){
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
								}
							}
							if( exception ) std::rethrow_exception(exception);
						}
					}
				}
			}
		}
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
		for( int i_readmap = 0; i_readmap < opt.rif_files.size(); ++i_readmap ){
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
				#ifdef USE_OPENMP
				#pragma omp critical
				#endif
				std::cout << "loaded RIF score for resl " << F(7,3,RESLS[i_readmap])
				          << " raw cart_resl: " << F(7,3,rif_ptr->cart_resl() )
				          << " raw ang_resl: " << F(7,3,rif_ptr->ang_resl() ) << std::endl;

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
	}


    if( 0 == opt.scaffold_fnames.size() ){
        std::cout << "WARNING: NO SCAFFOLDS!!!!!!" << std::endl;
    }

	for( int iscaff = 0; iscaff < opt.scaffold_fnames.size(); ++iscaff )
	{
		std::string scaff_fname = opt.scaffold_fnames.at(iscaff);
		std::vector<std::string> scaffold_sequence_glob0;				// Scaffold sequence in name3 space
		utility::vector1<core::Size> scaffold_res;//, scaffold_res_all; // Seqposs of residues to design, default whole scaffold
		try {

			runtime_assert( rot_index_p );
			std::string scafftag = utility::file_basename( utility::file::file_basename( scaff_fname ) );

			std::cout << "/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl;
			std::cout << "/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl;
			std::cout << "//////   begin scaffold " << scafftag << " " << iscaff << " of " << opt.scaffold_fnames.size() << std::endl;
			std::cout << "/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl;
			std::cout << "/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl;




			ScaffoldProviderOP scaffold_provider = get_scaffold_provider(
				iscaff,
				rot_index_p,
				opt,
				make2bopts,
				rotrf_table_manager);

			ScaffoldIndex rep_si = scaffold_provider->get_representative_scaffold_index();
			ScaffoldDataCacheOP rep_data_cache = scaffold_provider->get_data_cache_slow( rep_si );
			assert(rep_data_cache);

			// debugging info
			scaffold_sequence_glob0 = *(rep_data_cache->scaffold_sequence_glob0_p);
			scaffold_res = *(rep_data_cache->scaffold_res_p);

			// needed for scene
			float rep_scaff_radius = rep_data_cache->scaff_radius;

			// needed for cout
			float rep_scaff_redundancy_filter_rg = rep_data_cache->scaff_redundancy_filter_rg;
			Eigen::Vector3f rep_scaffold_center = rep_data_cache->scaffold_center;

			float rep_redundancy_filter_rg = std::min( rep_scaff_redundancy_filter_rg, target_redundancy_filter_rg );
			std::cout << "using redundancy_filter_rg: ~" << rep_redundancy_filter_rg << std::endl;


			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			print_header( "setup scene from scaffold and target" );
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


			// SOMETHING WRONG, SCORES OFF BY A LITTLE
			// setup objectives, moved into scaffold loop to guarantee clean slate for each scaff...
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



			ScenePtr scene_minimal( scene_prototype->clone_deep() );
			scene_minimal->add_actor( 0, VoxelActor(target_bounding_by_atype) );


			// utility_exit_with_message("FOO");

			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			print_header( "setup director based on scaffold and target sizes" ); //////////////////////////////////////////////////////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			DirectorBase director; {
				F3 target_center = pose_center(target);
				float const body_radius = std::min( rep_scaff_radius, rif_radius );
				double const cart_grid = opt.resl0*opt.hsearch_scale_factor/sqrt(3); // 1.5 is a big hack here.... 2 would be more "correct"
				double const hackysin = std::min( 1.0, opt.resl0*opt.hsearch_scale_factor/2.0/ body_radius );
				runtime_assert( hackysin > 0.0 );
				double const rot_resl_deg0 = asin( hackysin ) * 180.0 / M_PI;
				int nside = std::ceil( opt.search_diameter / cart_grid );
				std::cout << "search dia.    : " <<  opt.search_diameter << std::endl;
				std::cout << "nside          : " << nside        << std::endl;
				std::cout << "resl0:           " << opt.resl0 << std::endl;
				std::cout << "body_radius:     " << body_radius << std::endl;
				std::cout << "rif_radius:      " << rif_radius << std::endl;
				std::cout << "scaffold_radius: " << rep_scaff_radius << std::endl;
				std::cout << "cart_grid:       " << cart_grid  << std::endl;
				std::cout << "rot_resl_deg0:   " << rot_resl_deg0 << std::endl;
				I3 nc( nside, nside, nside );
				F3 lb = target_center + F3( -cart_grid*nside/2.0, -cart_grid*nside/2.0, -cart_grid*nside/2.0 );
				F3 ub = target_center + F3(  cart_grid*nside/2.0,  cart_grid*nside/2.0,  cart_grid*nside/2.0 );
				std::cout << "cart grid ub " << ub << std::endl;
				std::cout << "cart grid lb " << lb << std::endl;
				std::cout << "(ub-lb/nc) = " << ((ub-lb)/nc.template cast<float>()) << std::endl;
				std::cout << "cartcen to corner (cart. covering radius): " << sqrt(3.0)*cart_grid/2.0 << std::endl;
				shared_ptr<RifDockNestDirector> nest_director = make_shared<RifDockNestDirector>( rot_resl_deg0, lb, ub, nc, 1 );
				std::cout << "NestDirector:" << endl << *nest_director << endl;
				std::cout << "nest size0:    " << nest_director->size(0, RifDockIndex()).nest_index << std::endl;
				std::cout << "size of search space: ~" << float(nest_director->size(0, RifDockIndex()).nest_index)*1024.0*1024.0*1024.0 << " grid points" << std::endl;

				shared_ptr<RifDockScaffoldDirector> scaff_director = make_shared<RifDockScaffoldDirector>(scaffold_provider, 1 );

				std::vector<DirectorBase> director_list;
				director_list.push_back( nest_director );  // Nest director must come first!!!!
				director_list.push_back( scaff_director );

				director = make_shared<RifDockDirector>(director_list);
			}

			// Longxing
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			print_header( "perform test with scaffold in original position" ); //////////////////////////////////////////////////////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    		rep_data_cache->setup_onebody_tables( rot_index_p, opt);
			cout << std::endl;
			cout << "scores for scaffold in original position: " << std::endl;
			{
				EigenXform x(EigenXform::Identity());
				x.translation() = rep_scaffold_center;
				director->set_scene( RifDockIndex(0, rep_si), 0, *scene_minimal);
				scene_minimal->set_position(1,x);
				for(int i = 0; i < RESLS.size(); ++i){
					std::vector<float> sc = objectives[i]->scores(*scene_minimal);
					cout << "input bounding score " << i << " " << F(7,3,RESLS[i]) << " "
					     << F( 7, 3, sc[0]+sc[1] ) << " "
					     << F( 7, 3, sc[0]       ) << " "
					     << F( 7, 3, sc[1]       ) << endl;
				}

			}

			std::vector< ScenePtr > scene_pt( omp_max_threads_1() );
			BOOST_FOREACH( ScenePtr & s, scene_pt ) s = scene_minimal->clone_deep();

			RifDockData rdd {
						opt,
						RESLS,
						director,
						scene_pt,
						scene_minimal,
						target_simple_atoms,
						target_field_by_atype,
						&target_bounding_by_atype,
						&target_donors,
 						&target_acceptors,
 						target_redundancy_filter_rg,
 						target,
 						rot_index_p,
 						rotrf_table_manager,
 						objectives,
 						packing_objective,
 						packopts,
 						rif_ptrs,
 						rso_config,
 						rif_factory,
 						scaffold_provider
			};




			std::vector< SearchPointWithRots > packed_results;
			int64_t non0_space_size = 0;
			int64_t npack = 0;
			int64_t total_search_effort = 0;
			{
		        std::chrono::time_point<std::chrono::high_resolution_clock> start_rif = std::chrono::high_resolution_clock::now();

				///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				print_header( "perform hierarchical search" ); ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			    shared_ptr< std::vector< SearchPointWithRots > > hsearch_results_p; 

				{
					HSearchData data {
						total_search_effort,
						non0_space_size,
					};


					bool hsearch_success = call_hsearch_protocol( rdd, data, hsearch_results_p );
					if ( ! hsearch_success ) continue;
				}

				std::chrono::duration<double> elapsed_seconds_rif = std::chrono::high_resolution_clock::now()-start_rif;
				time_rif += elapsed_seconds_rif.count();


				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				//////////////////////////////////////////////         HACK PACK           /////////////////////////////////////////////////////////////
				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		        std::chrono::time_point<std::chrono::high_resolution_clock> start_pack = std::chrono::high_resolution_clock::now();
				

				scaffold_provider->set_fa_mode(true); // for legacy reasons this gets set here
		        if (opt.hack_pack) {
		        	hack_pack( hsearch_results_p, packed_results, rdd, total_search_effort, npack );
		        } else {
		        	packed_results = *hsearch_results_p;
		        }


				std::chrono::duration<double> elapsed_seconds_pack = std::chrono::high_resolution_clock::now()-start_pack;
				time_pck += elapsed_seconds_pack.count();
			}
			// std::vector< SearchPointWithRots > & packed_results = *packed_results_p;


			bool const do_rosetta_score = opt.rosetta_score_fraction > 0 || opt.rosetta_score_then_min_below_thresh > -9e8 || opt.rosetta_score_at_least > 0;

			if( do_rosetta_score && opt.hack_pack ){

				rosetta_rescore( packed_results, rdd, total_search_effort, time_ros );

			}


			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			print_header( "compile and filter results" ); ///////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			
			std::vector< RifDockResult > selected_results, allresults;
			compile_and_filter_results( packed_results, selected_results, allresults, rdd, dump_lock );



			std::cout << "allresults.size(): " << allresults.size() << " selected_results.size(): " << selected_results.size() << std::endl;

			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			print_header( "timing info" ); //////////////////////////////////////////////////////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			std::cout<<"total RIF     time: "<<KMGT(time_rif)<<" fraction: "<<time_rif/(time_rif+time_pck+time_ros)<<std::endl;
			std::cout<<"total Pack    time: "<<KMGT(time_pck)<<" fraction: "<<time_pck/(time_rif+time_pck+time_ros)<<std::endl;
			std::cout<<"total Rosetta time: "<<KMGT(time_ros)<<" fraction: "<<time_ros/(time_rif+time_pck+time_ros)<<std::endl;			

			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			print_header( "output results" ); //////////////////////////////////////////////////////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


			BOOST_FOREACH( RifDockResult const & r, allresults ){
				; // nothing with all results ATM
			}

			output_results(selected_results, rdd, dokout, npack);
			
			

		} catch( std::exception const & ex ) {
			std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
			std::cout << "error (below) on scaffold " << scaff_fname << " (will continue with others, if any)" << std::endl;
			std::cout << ex.what() << std::endl;
			std::cout << "scene residue numering (may help debug):" << std::endl;
			for( int i = 1; i <= scaffold_res.size(); ++i ){
				std::cout << "scene res numbering: " << i-1 << " " << scaffold_sequence_glob0.at(scaffold_res[i]-1) << " pose number: " << scaffold_res[i] << std::endl;
			}
			std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
		} catch ( ... ) {
			std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
			std::cout << "unknown error on scaffold " << scaff_fname << ", will continue with others, if any." << std::endl;
			std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
		}


	} // end scaffold loop


	dokout.close();






	#ifdef USE_OPENMP
		omp_destroy_lock( &cout_lock );
		omp_destroy_lock( &dump_lock );
	#endif

	std::cout << "rif_dock_test_DONE" << std::endl;

	return 0;
 }
