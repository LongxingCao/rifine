// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols

#include <riflib/rotamer_energy_tables.hh>
#include <riflib/util.hh>


#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/scoring/ScoreFunction.hh>

#include <utility/file/file_sys_util.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

#include <ObjexxFCL/format.hh>

#include <scheme/actor/Atom.hh>
#include <scheme/actor/BackboneActor.hh>
#include <scheme/chemical/RotamerIndex.hh>
#include <scheme/rosetta/score/RosettaField.hh>

#include <riflib/rotamer_energy_tables.hh>
#include <riflib/EtableParams_init.hh>
#include <riflib/RotamerGenerator.hh>
#include <riflib/util.hh>
#include <riflib/rosetta_field.hh>

#include <boost/multi_array.hpp>

#include <exception>
#include <stdexcept>

namespace devel {
namespace scheme {

using ObjexxFCL::format::I;

void get_onebody_rotamer_energies(
	core::pose::Pose const & scaffold,
	utility::vector1<core::Size> const & scaffold_res,
	devel::scheme::RotamerIndex const & rot_index,
	std::vector<std::vector<float> > & scaffold_onebody_rotamer_energies,
	std::vector<std::string> const & cachepath,
	std::string const & cachefile,
	bool replace_with_ala
){
	utility::io::izstream in;
	std::string cachefile_found = devel::scheme::open_for_read_on_path( cachepath, cachefile, in );
	if( cachefile.size() && cachefile_found.size() ){
		std::cout << "reading onebody energies from: " << cachefile << std::endl;
		// utility::io::izstream in( cachefile );
		size_t s1,s2;
		in.read((char*)&s1,sizeof(size_t));
		in.read((char*)&s2,sizeof(size_t));
		runtime_assert( 0 < s1 && s1 < 99999 );
		runtime_assert( 0 < s2 && s2 < 99999 );
		// std::cout << "nres: " << s1 << " nrot: " << s2 << std::endl;
		scaffold_onebody_rotamer_energies.resize( s1, std::vector<float>(s2) );
		for( size_t i = 0; i < s1; ++i ){
			for( size_t j = 0; j < s2; ++j ){
				in.read((char*)(&scaffold_onebody_rotamer_energies[i][j]),sizeof(float));
				runtime_assert( in.good() );
			}
		}
		std::string test;
		runtime_assert_msg( !(in >> test), "something left in buffer from: "+cachefile_found ); // nothing left
		in.close();
	} else {
		devel::scheme::compute_onebody_rotamer_energies(
			scaffold,
			scaffold_res,
			rot_index,
			scaffold_onebody_rotamer_energies,
					// this doesn't seem to work the way I expected...
					// false // do not mutate all to ALA!
			replace_with_ala
		);
		if( cachefile.size() ){
			std::cout << "saving onebody energies to: " << cachefile << std::endl;
			utility::io::ozstream out;//( cachefile );
			std::string writefile = open_for_write_on_path( cachepath, cachefile, out, true );
			runtime_assert( writefile.size() );
			size_t s1 = scaffold_onebody_rotamer_energies.size();
			size_t s2 = scaffold_onebody_rotamer_energies.front().size();
			out.write((char*)&s1,sizeof(size_t));
			out.write((char*)&s2,sizeof(size_t));
			for( size_t i = 0; i < s1; ++i ){
				for( size_t j = 0; j < s2; ++j ){
					out.write((char*)(&scaffold_onebody_rotamer_energies[i][j]),sizeof(float));
				}
			}
			out.close();
		}
	}

}


// this is the slow function, brian implemented a newer one, which is faster using low level code...
void
compute_onebody_rotamer_energies(
	core::pose::Pose const & pose,
	utility::vector1<core::Size> const & scaffold_res,
	RotamerIndex const & rot_index,
	std::vector<std::vector< float > > & onebody_rotamer_energies,
	bool replace_with_ala
){
	using devel::scheme::str;
	using devel::scheme::omp_max_threads_1;
	using devel::scheme::omp_thread_num_1;
	core::chemical::ResidueTypeSetCAP rts = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");
	core::conformation::ResidueOP ala = core::conformation::ResidueFactory::create_residue( rts.lock()->name_map("ALA") );


	std::vector<core::scoring::ScoreFunctionOP> score_func_per_thread(omp_max_threads_1());
	for( auto & score_func : score_func_per_thread ){
		// score_func = core::scoring::ScoreFunctionFactory::create_score_function( "talaris2014" );
		// score_func->set_etable( "FA_STANDARD_SOFT" );
		// // score_func->set_weight( core::scoring::fa_rep, score_func->get_weight(core::scoring::fa_rep)*0.67 );

		
		score_func = core::scoring::get_score_function();
		score_func->set_weight( core::scoring::fa_dun, score_func->get_weight(core::scoring::fa_dun)*0.67 );

		core::scoring::methods::EnergyMethodOptions opts = score_func->energy_method_options();
		core::scoring::hbonds::HBondOptions hopts = opts.hbond_options();
		hopts.use_hb_env_dep( false );
		opts.hbond_options( hopts );
		score_func->set_energy_method_options( opts );
	}

	// make all ala or gly
	core::pose::Pose bbone(pose);
	if( replace_with_ala ) pose_to_ala( bbone );
	double const base_score = score_func_per_thread.front()->score( bbone );
	// bbone.dump_pdb("bbone_test.pdb");


	onebody_rotamer_energies.resize( bbone.size() );
	std::vector<core::pose::Pose> pose_per_thread( omp_max_threads_1(), bbone );
	std::cout << "compute_onebody_rotamer_energies " << bbone.size() << "/" << rot_index.size();
	std::exception_ptr exception = nullptr;
	#ifdef USE_OPENMP
	#pragma omp parallel for schedule(dynamic,1)
	#endif
	for( int ir = 1; ir <= bbone.size(); ++ir ){
		if( exception ) continue;
		if( std::find(scaffold_res.begin(), scaffold_res.end(), ir) == scaffold_res.end() ){
			onebody_rotamer_energies[ir-1].resize( rot_index.size(), 999999.0 );
			continue;
		}
		try {
			core::pose::Pose & work_pose( pose_per_thread[ omp_thread_num_1()-1 ] );
			core::scoring::ScoreFunctionOP score_func = score_func_per_thread[ omp_thread_num_1()-1 ];
			onebody_rotamer_energies[ir-1].resize( rot_index.size(), 999999.0 );
			if( ! work_pose.residue(ir).is_protein()   ) continue;
			if(   work_pose.residue(ir).name3()=="GLY" ) continue;
			if(   work_pose.residue(ir).name3()=="PRO" ) continue;
			#ifdef USE_OPENMP
			#pragma omp critical
			#endif
			std::cout << (100.0*ir)/work_pose.size() << "% "; std::cout.flush();
			for( int jr = 0; jr < rot_index.size(); ++jr ){
				core::conformation::ResidueOP rot = core::conformation::ResidueFactory::create_residue( rts.lock()->name_map( rot_index.resname(jr) ) );
				work_pose.replace_residue( ir, *rot, true );
				for( int k = 0; k < rot_index.nchi(jr); ++k ){
					work_pose.set_chi( k+1, ir, rot_index.chi( jr, k ) );
				}
				onebody_rotamer_energies[ir-1][jr] = score_func->score( work_pose ) - base_score;
				// std::cout << "fa_dun " << ir << " " << jr << " "<< work_pose.energies().residue_total_energies(ir)[core::scoring::fa_dun] << std::endl;
				work_pose.replace_residue( ir, *ala, true );
				// if( jr > 2	 ) break;
			}
		} catch( ... ) {
			#ifdef USE_OPENMP
			#pragma omp critical
			#endif
			exception = std::current_exception();
		}
	}
	if( exception ) std::rethrow_exception(exception);

	std::cout << "compute_onebody_rotamer_energies_DONE" << std::endl;

	// for( int i = 0; i < onebody_rotamer_energies.size(); ++i ){
	// 	std::cout << "OBE " << i;
	// 	for( int j = 0; j < onebody_rotamer_energies[i].size(); ++j ){
	// 		std::cout << " " << onebody_rotamer_energies[i][j];
	// 	}
	// 	std::cout << std::endl;
	// }


}

/* 
    
    
    void
    compute_onebody_rotamer_energies(
                                     core::pose::Pose const & pose,
                                     utility::vector1<core::Size> const & scaffold_res,
                                     RotamerIndex const & rot_index,
                                     std::vector<std::vector< float > > & onebody_rotamer_energies,
                                     bool replace_with_ala
                                     ){
        using devel::scheme::str;
        using devel::scheme::omp_max_threads_1;
        using devel::scheme::omp_thread_num_1;
        
        // make all ala or gly
        core::pose::Pose bbone(pose);
        if( replace_with_ala ) pose_to_ala( bbone );
        std::vector<core::pose::Pose> pose_per_thread( omp_max_threads_1(), bbone );
        
        std::vector<core::scoring::ScoreFunctionOP> score_func_per_thread(omp_max_threads_1());
        utility::vector1<bool> true_vect(bbone.size(), true);
        
        //runtime_assert(rot_index.per_thread_rotamers_.size() > 0);
        int ala_rot = rot_index.ala_rot();
        
        core::scoring::get_score_function();	// get one first to prevent multi-thread crash
        #ifdef USE_OPENMP
        #pragma omp parallel for schedule(dynamic,1)
        #endif
        for( int i = 0; i < score_func_per_thread.size(); i++ ){
            
            // We specifically don't use the thread number here because we're filling in the per-thread stuff
            core::scoring::ScoreFunctionOP & score_func = score_func_per_thread[i];
            core::pose::Pose & work_pose = pose_per_thread[i];
            
            score_func = core::scoring::get_score_function();
            score_func->set_weight( core::scoring::fa_dun, score_func->get_weight(core::scoring::fa_dun)*0.67 );
            
            core::scoring::methods::EnergyMethodOptions opts = score_func->energy_method_options();
            core::scoring::hbonds::HBondOptions hopts = opts.hbond_options();
            hopts.use_hb_env_dep( false );
            hopts.decompose_bb_hb_into_pair_energies(true);
            opts.hbond_options( hopts );
            score_func->set_energy_method_options( opts );
            
            // all work poses are scored here
            score_func->score( work_pose );
            score_func->setup_for_packing( work_pose, true_vect, true_vect );
            
        }
        
        
        utility::vector1<core::Real> per_res_nbr_dist( bbone.size(), 0 );
        float max_res_radius = 0;
        for ( core::Size ir = 1; ir <= bbone.size(); ir++ ) {
            float radius = bbone.residue(ir).nbr_radius();
            per_res_nbr_dist[ir] = radius;
            max_res_radius = std::max<float>( max_res_radius, radius );
        }
        
        float max_rotamer_radius = rot_index.get_max_nbr_radius();
        float max_interaction_radius = score_func_per_thread.front()->info()->max_atomic_interaction_distance();
        float max_radius = max_rotamer_radius + max_res_radius + max_interaction_radius;
        
        // Copied from core/pack/packer_neighbors.cc
        core::conformation::PointGraphOP point_graph( new core::conformation::PointGraph );
        core::conformation::residue_point_graph_from_conformation( bbone.conformation(), *point_graph );
        core::conformation::find_neighbors<core::conformation::PointGraphVertexData,core::conformation::PointGraphEdgeData>(
                                                                                                                            point_graph, max_radius );
        
        utility::graph::Graph neighbor_graph( bbone.size() );
        
        for ( int ir = 1; ir <= bbone.size(); ir++ ) {
            for ( auto
                 iter = point_graph->get_vertex( ir ).const_upper_edge_list_begin(),
                 iter_end = point_graph->get_vertex( ir ).const_upper_edge_list_end();
                 iter != iter_end; ++iter ) {
                
                int other_pos = iter->upper_vertex();
                core::Real other_dist = per_res_nbr_dist[other_pos];
                core::Real dist = other_dist + max_rotamer_radius + max_interaction_radius;
                
                if ( iter->data().dsq() < dist*dist ) {
                    neighbor_graph.add_edge( ir, other_pos );
                }
            }
        }
        
        
        onebody_rotamer_energies.resize( bbone.size() );
        std::cout << "compute_onebody_rotamer_energies " << bbone.size() << "/" << rot_index.size() << " ";
        std::exception_ptr exception = nullptr;
        #ifdef USE_OPENMP
        #pragma omp parallel for schedule(dynamic,1)
        #endif
        for( int ir = 1; ir <= bbone.size(); ++ir ){
            if( exception ) continue;
            if( std::find(scaffold_res.begin(), scaffold_res.end(), ir) == scaffold_res.end() ){
                onebody_rotamer_energies[ir-1].resize( rot_index.size(), 12345.0 );
                continue;
            }
            try {
                core::pose::Pose & work_pose( pose_per_thread[ omp_thread_num_1()-1 ] );
                core::scoring::ScoreFunctionOP score_func = score_func_per_thread[ omp_thread_num_1()-1 ];
                onebody_rotamer_energies[ir-1].resize( rot_index.size(), 12345.0 );
                if( ! work_pose.residue(ir).is_protein()   ) continue;
                if(   work_pose.residue(ir).name3()=="GLY" ) continue;
                if(   work_pose.residue(ir).name3()=="PRO" ) continue;
                #ifdef USE_OPENMP
                #pragma omp critical
                #endif
                {
                    std::cout << (100.0*ir)/work_pose.size() << "% "; std::cout.flush();
                }
                
                
                // Let's make a rotset for the current position
                core::pack::rotamer_set::RotamerSet_ rotset;
                rotset.set_resid(ir);
                
                core::conformation::ResidueOP original_rot = work_pose.residue( ir ).clone();
                for ( int irot = 0; irot < rot_index.size(); irot++ ) {
                    core::conformation::ResidueOP pt_rot = rot_index.get_per_thread_rotamer( omp_thread_num(), irot );
                    work_pose.replace_residue( ir, *pt_rot, true );	// I give up, there's just so much you have to update without replace residue
                    rotset.add_rotamer( work_pose.residue(ir) );	// not cloning because yolo
                }
                work_pose.replace_residue( ir, *original_rot, true );
                // Add the current rotamer to the end to get the base score
                rotset.add_rotamer( *original_rot );
                int original_rot_index = rotset.num_rotamers();
                
                score_func->prepare_rotamers_for_packing(work_pose, rotset);
                
                // Now make the energy containers that rosetta wants
                //   ala    rotamers
                ObjexxFCL::FArray2D< core::PackerEnergy > pair_energy_table( 1, rotset.num_rotamers(), 0.0 );
                utility::vector1< core::PackerEnergy > onebody_energies( rotset.num_rotamers(), 0 );
                
                // First do the actual one-body energies
                for ( int irot = 1; irot <= rotset.num_rotamers(); irot++ ) {
                    core::scoring::EnergyMap emap;
                    score_func->eval_ci_1b( *(rotset.rotamer( irot )), work_pose, emap );
                    score_func->eval_cd_1b( *(rotset.rotamer( irot )), work_pose, emap );
                    onebody_energies[ irot ] += static_cast< core::PackerEnergy > (score_func->weights().dot( emap ));
                }
                score_func->evaluate_rotamer_intrares_energies( rotset, work_pose, onebody_energies );
                
                // Loop over interacting positions
                for ( utility::graph::Graph::EdgeListConstIter
                     iter = neighbor_graph.get_node( ir )->const_edge_list_begin(),
                     iter_end = neighbor_graph.get_node( ir )->const_edge_list_end();
                     iter != iter_end; ++iter ) {
                    
                    int seqpos = (*iter)->get_other_ind( ir );
                    
                    core::pack::rotamer_set::RotamerSet_ other_rotset;
                    other_rotset.set_resid(seqpos);
                    other_rotset.add_rotamer( *(work_pose.residue(seqpos).clone()));
                    score_func->prepare_rotamers_for_packing(work_pose, other_rotset);
                    
                    score_func->evaluate_rotamer_pair_energies( rotset, other_rotset, work_pose, pair_energy_table );
                    
                }
                
                
                // Long range interactions have their own graph structure
                for ( auto
                     lr_iter = score_func->long_range_energies_begin(),
                     lr_end = score_func->long_range_energies_end();
                     lr_iter != lr_end; ++lr_iter ) {
                    
                    core::scoring::LREnergyContainerCOP lrec = work_pose.energies().long_range_container( (*lr_iter)->long_range_type() );
                    if ( !lrec || lrec->empty() ) continue; // only score non-emtpy energies.
                    
                    for ( core::scoring::ResidueNeighborConstIteratorOP
                         rni = lrec->const_neighbor_iterator_begin(ir),
                         rniend = lrec->const_neighbor_iterator_end(ir);
                         (*rni) != (*rniend); ++(*rni) ) {
                        
                        core::Size seqpos = rni->neighbor_id();
                        runtime_assert( seqpos != ir );
                        
                        
                        core::pack::rotamer_set::RotamerSet_ other_rotset;
                        other_rotset.set_resid(seqpos);
                        other_rotset.add_rotamer( *(work_pose.residue(seqpos).clone()));
                        score_func->prepare_rotamers_for_packing(work_pose, other_rotset);
                        
                        (*lr_iter)->evaluate_rotamer_pair_energies( rotset, other_rotset, work_pose, *score_func, 
                                                                   score_func->weights(), pair_energy_table );
                        
                        
                    }
                }
                
                // Now we have the scores, lets fill the table
                
                // Everything is relative to the original, so get that score first
                float base_score = onebody_energies[original_rot_index] + pair_energy_table( 1, original_rot_index);
                
                for ( int irot = 0; irot < rot_index.size(); irot++ ) {
                    float raw_score = onebody_energies[irot+1] + pair_energy_table( 1, irot+1);
                    onebody_rotamer_energies[ir-1][irot] = raw_score - base_score;
                }
                
                
                
            } catch( ... ) {
                #ifdef USE_OPENMP
                #pragma omp critical
                #endif
                exception = std::current_exception();
            }
        }
        if( exception ) std::rethrow_exception(exception);
        
        std::cout << "compute_onebody_rotamer_energies_DONE" << std::endl;
        
        // for( int i = 0; i < onebody_rotamer_energies.size(); ++i ){
        // 	std::cout << "OBE " << i;
        // 	for( int j = 0; j < onebody_rotamer_energies[i].size(); ++j ){
        // 		std::cout << " " << onebody_rotamer_energies[i][j];
        // 	}
        // 	std::cout << std::endl;
        // }
        
        
    }

*/

void get_per_rotamer_rf_tables_one(
	devel::scheme::RotamerIndex const & rot_index,
	int irot,
	RotamerRFOpts opts,
	int const N_ATYPE,
	std::vector< std::vector<
		::scheme::shared_ptr<
			::scheme::objective::voxel::VoxelArray< 3, float >
		>
	> >	& fields_by_atype_per_rot,
	bool verbose
){
	typedef ::scheme::objective::voxel::FieldCache3D< float > FieldCache;
	typedef ::scheme::objective::voxel::BoundingFieldCache3D< float, ::scheme::objective::voxel::AggMin > BoundingFieldCache;
	typedef ::scheme::objective::voxel::VoxelArray< 3, float > VoxelArray;
	typedef ::scheme::actor::Atom< Eigen::Vector3f > SchemeAtom;
	typedef ::scheme::util::SimpleArray<3,float> F3;

	using ::devel::scheme::str;
	using ::devel::scheme::omp_thread_num_1;

	runtime_assert_msg( N_ATYPE==21, "N_ATYPE assumed to be 21" );

	std::vector< shared_ptr< VoxelArray > > & field_by_atype = fields_by_atype_per_rot[irot];
	// part 1 of total sync hack... assume ready iff size == N_ATYPE+1
	field_by_atype.resize(N_ATYPE+2,nullptr);

	if( rot_index.structural_parent(irot) != irot ){
		field_by_atype.resize(N_ATYPE+1,nullptr);
		return; // will use parent fields
	}

	std::string const & resname = rot_index.rotamers_[irot].resname_;
	// #ifdef USE_OPENMP
	// #pragma omp critical
	// #endif
	// cout << "build rotamer RF for " << resname << " " << irot << endl;
	std::vector<SchemeAtom> atoms;
	// use only heavy atoms beyond the CB (which is #3 here)
	for( int i = 4; i < rot_index.rotamers_[irot].nheavyatoms; ++i ){
		if( rot_index.rotamers_[irot].atoms_[i].type() > N_ATYPE ) continue;
		atoms.push_back( rot_index.rotamers_[irot].atoms_[i] );
		// cout << i << " " << rot_index.rotamers_[irot].atoms_[i] << " " << rot_index.rotamers_[irot].atoms_[i].type() << endl;
	}
	if( atoms.size() == 0 ){
		// if no sidechain atoms, leave empty and continue
		field_by_atype.resize(N_ATYPE+1,nullptr);
		return;
	}

	::scheme::rosetta::score::RosettaField< SchemeAtom, devel::scheme::EtableParamsInit > rosetta_field( atoms );

	F3 lb(9e9,9e9,9e9),ub(-9e9,-9e9,-9e9);
	uint64_t target_atoms_hash = 0;
	for( auto const & a : atoms ){
		lb = lb.min(a.position());
		ub = ub.max(a.position());
	}
	// cout << "rosetta_field lb " << lb << endl;
	// cout << "rosetta_field ub " << ub << endl;

	std::string cachefile = opts.data_dir
	                 +"/__rotamer_"+str(irot,3)+resname
	                 +"__resl"      +boost::lexical_cast<std::string>(opts.field_resl)
	                 +"__spread"    +boost::lexical_cast<std::string>(opts.field_spread)
		             +"__oversamp"  +boost::lexical_cast<std::string>(opts.oversample);
    if( opts.scale_atr != 1.0f ){
    	cachefile = cachefile +"__scaleatr" + boost::lexical_cast<std::string>(opts.scale_atr);
    }
	cachefile = cachefile + "__21atypes.rosetta_field.gz";

	// must figure out how to save/load bounding grid without need for original

	// #ifdef USE_OPENMP
	// #pragma omp parallel for schedule(dynamic,1)
	// #endif

	if( utility::file::file_exists(cachefile) ){
		if( verbose ){
			#ifdef USE_OPENMP
			#pragma omp critical
			#endif
			std::cout << "thread " << I(3,omp_thread_num_1()) << " init  rot_rf_table CACHE AT " << cachefile << std::endl;
		}
		utility::io::izstream in( cachefile, std::ios::binary );
		for( int itype = 1; itype <= N_ATYPE; ++itype ){
			field_by_atype[itype] = ::scheme::make_shared< VoxelArray >( lb-6.0f-opts.field_spread, ub+6.0f+opts.field_spread, opts.field_resl );
			field_by_atype[itype]->load(in);
			// std::cout << *field_by_atype[itype] << std::endl;
		}
		in.close();
	} else {
		#ifdef USE_OPENMP
		#pragma omp critical
		#endif
		std::cout << "thread " << I(3,omp_thread_num_1()) << " init  rot_rf_table CACHE TO " << cachefile << std::endl;
		utility::io::ozstream out( cachefile , std::ios::binary );
		for( int itype = 1; itype <= N_ATYPE; ++itype ){
			// EtableParamsInit uses scheme atom numbering
			::scheme::rosetta::score::RosettaFieldAtype< SchemeAtom, devel::scheme::EtableParamsInit > rfa( rosetta_field, itype );
			::scheme::shared_ptr<VoxelArray> tmp = ::scheme::make_shared<FieldCache>( rfa, lb-6.0f, ub+6.0f, opts.field_resl, "", false, opts.oversample );
			if( opts.field_spread < 0.01 ){
				field_by_atype[itype] = tmp;
			} else {
				field_by_atype[itype] = ::scheme::make_shared< BoundingFieldCache >( *tmp, opts.field_spread, opts.field_resl );
			}
			if( opts.scale_atr != 1.0f ){
				size_t const nelem = field_by_atype[itype]->num_elements();
				for( size_t k = 0; k < nelem; ++k ){
					if( field_by_atype[itype]->data()[k] < 0.0 ){
						field_by_atype[itype]->data()[k] *= opts.scale_atr;
					}
				}
			}
			runtime_assert( utility::file::create_directory_recursive(opts.data_dir) );
			field_by_atype[itype]->save( out );
			// std::cout << "    " << I(3,irot) << " " << I(2,itype) << " " << *field_by_atype[itype] << std::endl;
		}
		out.close();
	}

	// part 2 of total sync hack (see above)
	#ifdef USE_OPENMP
	#pragma omp critical
	field_by_atype.resize(N_ATYPE+1);
	#endif

}


void get_per_rotamer_rf_tables(
	devel::scheme::RotamerIndex const & rot_index,
	RotamerRFOpts opts,
	int const N_ATYPE,
	std::vector< std::vector<
		::scheme::shared_ptr<
			::scheme::objective::voxel::VoxelArray< 3, float >
		>
	> >	& fields_by_atype_per_rot
){
	typedef ::scheme::objective::voxel::FieldCache3D< float > FieldCache;
	typedef ::scheme::objective::voxel::BoundingFieldCache3D< float, ::scheme::objective::voxel::AggMin > BoundingFieldCache;
	typedef ::scheme::objective::voxel::VoxelArray< 3, float > VoxelArray;
	typedef ::scheme::actor::Atom< Eigen::Vector3f > SchemeAtom;
	typedef ::scheme::util::SimpleArray<3,float> F3;

	using ::devel::scheme::str;
	using ::devel::scheme::omp_thread_num_1;

	runtime_assert_msg( N_ATYPE==21, "N_ATYPE assumed to be 21" );

	std::exception_ptr exception = nullptr;
	#ifdef USE_OPENMP
	#pragma omp parallel for schedule(dynamic,1)
	#endif
	for( int irot = 0; irot < rot_index.size(); ++irot ){
		if( exception ) continue;
		try{
			bool verbose = irot==0;
			get_per_rotamer_rf_tables_one( rot_index, irot, opts, N_ATYPE, fields_by_atype_per_rot, verbose );
		} catch( ... ) {
			#ifdef USE_OPENMP
			#pragma omp critical
			#endif
			exception = std::current_exception();
		}
	}
	if( exception ) std::rethrow_exception(exception);

	// fill in pointers to duplicate tables for rotamers with redundand heavy atom positions
	for( int irot = 0; irot < rot_index.size(); ++irot ){
		if( rot_index.structural_parent(irot) != irot ){
			for( int itype = 1; itype <= N_ATYPE; ++itype ){
				fields_by_atype_per_rot[ irot ][ itype ] = fields_by_atype_per_rot[ rot_index.structural_parent(irot) ][ itype ];
			}
		}
	}
	// for( int irot = 0; irot < rot_index.size(); ++irot ){ // ala and gly are null
	// 	for( int itype = 1; itype <= N_ATYPE; ++itype ){
	// 		cout << irot << " " << itype << endl;
	// 		runtime_assert( fields_by_atype_per_rot[ irot ][ itype ] );
	// 	}
	// }


}






void
make_twobody_tables(
	core::pose::Pose const & scaffold,
	devel::scheme::RotamerIndex const & rot_index,
	std::vector<std::vector<float> > const & onebody_energies,
	RotamerRFTablesManager & rotrfmanager,
	MakeTwobodyOpts opts,
	::scheme::objective::storage::TwoBodyTable<float> & twob
){
	// typedef ::scheme::objective::voxel::VoxelArray< 3, float > VoxelArray;
	// typedef Eigen::Transform<float,3,Eigen::AffineCompact> EigenXform;
	typedef ::scheme::actor::BackboneActor<EigenXform> BackboneActor;
	typedef ::scheme::actor::Atom< Eigen::Vector3f > SchemeAtom;

	runtime_assert( onebody_energies.size() == scaffold.size() );
 	// boost::multi_array<float,2> obe( boost::extents[scaffold.size()][rot_index.size()] );
	runtime_assert( twob.onebody_.shape()[0] == scaffold.size() );
	runtime_assert( twob.onebody_.shape()[1] == rot_index.size() );
	{
		runtime_assert( onebody_energies.size() == scaffold.size() );
		for( int i = 0; i < scaffold.size(); ++i ){
			runtime_assert( onebody_energies[i].size() == rot_index.size() );
			for( int j = 0; j < rot_index.size(); ++j ){
				twob.onebody_[i][j] = onebody_energies[i][j];
			}
		}
	}
	int avg = 0;
	for( int i = 0; i < scaffold.size(); ++i ){
		int count = 0;
		for( int j = 0; j < rot_index.size(); ++j ){
			count += ( twob.onebody_[i][j] < 0.0 );
		}
		// cout << i+1 << " " << count << endl;
		avg += count;
	}
	// cout << "avg count = " << (float)avg/scaffold.size() << endl;

	twob.init_onebody_filter( opts.onebody_threshold );

	double const dthresh2 = opts.distance_cut * opts.distance_cut;

	std::exception_ptr exception = nullptr;
	#ifdef USE_OPENMP
	#pragma omp parallel for schedule(dynamic,1)
	#endif
	for( int ir = 0; ir < scaffold.size(); ++ir ){
		if( exception ) continue;
		try {
			if( !scaffold.residue(ir+1).is_protein() ) continue;
			BackboneActor bbi( scaffold.residue(ir+1).xyz("N"), scaffold.residue(ir+1).xyz("CA"), scaffold.residue(ir+1).xyz("C") );

			// #ifdef USE_OPENMP
			// #pragma omp critical
			// #endif
			// std::cout << "compute twobody for ir-* " << ir << std::endl;

			for( int jr = 0; jr < ir; ++jr ){
				if( !scaffold.residue(jr+1).is_protein() ) continue;

				double dis2 = scaffold.residue(ir+1).xyz("CA").distance_squared( scaffold.residue(jr+1).xyz("CA") );
				// double dthr = scaffold.residue(ir).nbr_radius() +                   scaffold.residue(jr+1).nbr_radius();

				if( dis2 > dthresh2 ){
					// twob.clear_twobody(ir,jr); // no need
					continue;
				}

				BackboneActor bbj( scaffold.residue(jr+1).xyz("N"), scaffold.residue(jr+1).xyz("CA"), scaffold.residue(jr+1).xyz("C") );

				twob.init_twobody(ir,jr);

				EigenXform X2i = bbi.position().inverse() * bbj.position();
				EigenXform X2j = bbj.position().inverse() * bbi.position();
				auto const & to_sp( rot_index.to_structural_parent_frame_ );

				float minscore=9e9, maxscore=-9e9;
				std::vector<int> randirotsel( twob.nsel_[ir] );
				for( int q=0; q<randirotsel.size(); ++q ) randirotsel[q] = randirotsel.size()-1-q;
				// this is fucked somehow... random would be better, reversed helps some
				// random_permutation( randirotsel, numeric::random::rg() );
				for( int irotsel_irand = 0; irotsel_irand < twob.nsel_[ir]; ++irotsel_irand ){
					int irotsel = randirotsel[ irotsel_irand ];
					int irot = twob.sel2all_[ir][irotsel];
					runtime_assert( irot >= 0 );
					for( int jrotsel = 0; jrotsel < twob.nsel_[jr]; ++jrotsel ){
						int jrot = twob.sel2all_[jr][jrotsel];
						runtime_assert( jrot >= 0 );

						twob.twobody_[ir][jr][irotsel][jrotsel] = 0.0; // set to 0

						// get lj, sol
						if( rot_index.nheavyatoms(irot) > rot_index.nheavyatoms(jrot) ){ // irot is bigger
							if( rotrfmanager.get_rotamer_rf_tables(irot)[ 1 ] ){
								for( int ja = 4; ja < rot_index.nheavyatoms(jrot); ++ja ){ // use only heavy atoms beyond the CB (which is #3 here)
									int jatype = rot_index.rotamers_[ jrot ].atoms_[ja].type();
									runtime_assert( jatype > 0 && jatype < 22 );
									Eigen::Vector3f pos_ja = to_sp.at(irot) * X2i * rot_index.rotamers_[ jrot ].atoms_[ja].position();
									// runtime_assert( rotrfmanager.get_rotamer_rf_tables(irot).size() );
									float const atomscore = rotrfmanager.get_rotamer_rf_tables(irot).at( jatype )->at( pos_ja );
									runtime_assert_msg( atomscore < 9999.0, "very high atomscore" );
									twob.twobody_[ir][jr][irotsel][jrotsel] += atomscore;
								}
							} else {
								if( rot_index.resname(irot)!="ALA"&&rot_index.resname(irot)!="GLY" ){
									utility_exit_with_message( "no rotrf table for "+str(irot)+" / "+ str(ir)+rot_index.resname(irot)
									    + " other is" + str(jr)+rot_index.resname(jrot) );
								}
							}
						} else {
							if( rotrfmanager.get_rotamer_rf_tables(jrot)[ 1 ] ){
								for( int ia = 4; ia < rot_index.nheavyatoms(irot); ++ia ){ // use only heavy atoms beyond the CB (which is #3 here)
									int iatype = rot_index.rotamers_[ irot ].atoms_[ia].type();
									if( iatype > 21 ){
										std::cout << iatype << " " << irot << " " << ia << " " << rot_index.rotamers_[irot].atoms_[ia].data().atomname << " "
										          << rot_index.rotamers_[irot].resname_ << " " << rot_index.nheavyatoms(irot) << std::endl;
									}
									runtime_assert( iatype > 0 && iatype < 22 );
									Eigen::Vector3f pos_ia = to_sp.at(jrot) * X2j * rot_index.rotamers_[ irot ].atoms_[ia].position();
									// runtime_assert( rotrfmanager.get_rotamer_rf_tables(jrot).size() );
									float const atomscore = rotrfmanager.get_rotamer_rf_tables(jrot).at( iatype )->at( pos_ia );
									runtime_assert_msg( atomscore < 9999.0, "very high atomscore" );
									twob.twobody_[ir][jr][irotsel][jrotsel] += atomscore;
								}
							} else {
								if( rot_index.resname(jrot)!="ALA"&&rot_index.resname(jrot)!="GLY" ){
									utility_exit_with_message( "no rotrf table for "+str(jrot)+" / "+str(jr)+rot_index.resname(jrot)
									    + " other is" + str(ir)+rot_index.resname(irot) );
								}
							}
						}

						// this is basically a copy of what's in ScoreRotamerVsTarget, without the multidentate stuff
						if( rot_index.rotamer(irot).acceptors_.size() > 0 ||
							rot_index.rotamer(irot).donors_   .size() > 0 )
						{
							float hbscore = 0.0;
							for( int i_hr_rot_acc = 0; i_hr_rot_acc < rot_index.rotamer(irot).acceptors_.size(); ++i_hr_rot_acc )
							{
								HBondRay hr_rot_acc = rot_index.rotamer(irot).acceptors_[i_hr_rot_acc];
								Eigen::Vector3f dirpos = hr_rot_acc.horb_cen + hr_rot_acc.direction;
								hr_rot_acc.horb_cen  = X2j * hr_rot_acc.horb_cen;
								hr_rot_acc.direction = X2j * dirpos - hr_rot_acc.horb_cen;
								for( int i_hr_tgt_don = 0; i_hr_tgt_don < rot_index.rotamer(jrot).donors_.size(); ++i_hr_tgt_don )
								{
									HBondRay const & hr_tgt_don = rot_index.rotamer(jrot).donors_[i_hr_tgt_don];
									float const thishb = score_hbond_rays( hr_tgt_don, hr_rot_acc );
									hbscore += thishb * opts.hbond_weight;
								}
							}
							for( int i_hr_rot_don = 0; i_hr_rot_don < rot_index.rotamer(irot).donors_.size(); ++i_hr_rot_don )
							{
								HBondRay hr_rot_don = rot_index.rotamer(irot).donors_[i_hr_rot_don];
								Eigen::Vector3f dirpos = hr_rot_don.horb_cen + hr_rot_don.direction;
								hr_rot_don.horb_cen  = X2j * hr_rot_don.horb_cen;
								hr_rot_don.direction = X2j * dirpos - hr_rot_don.horb_cen;
								for( int i_hr_tgt_acc = 0; i_hr_tgt_acc < rot_index.rotamer(jrot).acceptors_.size(); ++i_hr_tgt_acc )
								{
									HBondRay const & hr_tgt_acc = rot_index.rotamer(jrot).acceptors_[i_hr_tgt_acc];
									float const thishb = score_hbond_rays( hr_rot_don, hr_tgt_acc );
									hbscore += thishb * opts.hbond_weight;
								}
							}
							twob.twobody_[ir][jr][irotsel][jrotsel] += hbscore;
						}



						if( twob.twobody_[ir][jr][irotsel][jrotsel] > 12345.0 ){
							twob.twobody_[ir][jr][irotsel][jrotsel] = 12345.0;
						}

						minscore = std::min( minscore, twob.twobody_[ir][jr][irotsel][jrotsel] );
						maxscore = std::max( maxscore, twob.twobody_[ir][jr][irotsel][jrotsel] );
					}
				}

				if( minscore > -0.01 && maxscore < 0.01 ){
					twob.clear_twobody( ir, jr );
				} else {
					// using namespace ObjexxFCL::format;
					// #pragma omp critical
					// std::cout << I(3,ir) << " " << I(3,jr) << " " << F(9,1,maxscore) << std::endl;
				}

			}
		} catch( ... ) {
			#ifdef USE_OPENMP
			#pragma omp critical
			#endif
			exception = std::current_exception();
		}
	}
	if( exception ) std::rethrow_exception(exception);

}

void
get_twobody_tables(
	std::vector<std::string> const & cachepath,
	std::string const & cachefile,
	std::string & description, // in/out
	core::pose::Pose const & scaffold,
	devel::scheme::RotamerIndex const & rot_index,
	std::vector<std::vector<float> > const & onebody_energies,
	RotamerRFTablesManager & rotrfmanager,
	MakeTwobodyOpts opts,
	::scheme::objective::storage::TwoBodyTable<float> & twob
){
	utility::io::izstream in;
	std::string cachefile_found;
	if( cachefile.size() ) cachefile_found = devel::scheme::open_for_read_on_path( cachepath, cachefile, in );
	if( cachefile.size() && cachefile_found.size() ){
		std::cout << "reading twobody energies from: " << cachefile_found << std::endl;
		twob.load( in, description );
		in.close();
	} else {
		twob.init( scaffold.size(), rot_index.size() );
		make_twobody_tables( scaffold, rot_index, onebody_energies, rotrfmanager, opts, twob );
		if( cachefile.size() ) std::cout << "created twobody energies and saving to: " << cachefile << std::endl;
		if( description=="" ) description = "No description, Will sucks. Complain to willsheffler@gmail.com\n";
		utility::io::ozstream out;//( cachefile );
		devel::scheme::open_for_write_on_path( cachepath, cachefile, out, true );
		twob.save( out, description );
		out.close();
	}


	if ( opts.favorable_2body_multiplier != 1 ) {
		for ( uint64_t i = 0; i < twob.twobody_.size(); i++ ) {
			for ( uint64_t j = 0; j < twob.twobody_[i].size(); j++ ) {
				for ( uint64_t k = 0; k < twob.twobody_[i][j].size(); k++ ) {
					for ( uint64_t l = 0; l < twob.twobody_[i][j][k].size(); l++ ) {
						float val = twob.twobody_[i][j][k][l];
						if ( val < 0 ) {
							twob.twobody_[i][j][k][l] = val * opts.favorable_2body_multiplier;
						}
					}
				}
			}
		}
	}


}



}}

