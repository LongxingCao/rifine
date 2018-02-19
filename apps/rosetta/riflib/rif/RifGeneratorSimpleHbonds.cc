// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.



#include <riflib/rif/RifGeneratorSimpleHbonds.hh>




	#include <ObjexxFCL/FArray3D.hh>
	#include <ObjexxFCL/format.hh>

	#include <boost/assign/std/vector.hpp>
	#include <boost/foreach.hpp>
	#include <boost/lexical_cast.hpp>
	#include <boost/math/special_functions/sign.hpp>
	#include <boost/random/mersenne_twister.hpp>
	#include <boost/random/uniform_real.hpp>

	#include <core/id/AtomID.hh>
	#include <core/pose/Pose.hh>
        #include <core/pose/variant_util.hh>
	#include <core/scoring/motif/util.hh>

	#include <devel/init.hh>
	#include <riflib/HBondedPairGenerator.hh>
	#include <riflib/RotamerGenerator.hh>
	#include <riflib/util.hh>
	#include <riflib/rif/make_hbond_geometries.hh>

	#include <map>

	#include <scheme/actor/Atom.hh>
	#include <scheme/actor/BackboneActor.hh>
	#include <scheme/actor/VoxelActor.hh>
	#include <scheme/chemical/RotamerIndex.hh>
	#include <scheme/kinematics/Director.hh>
	#include <scheme/kinematics/Scene.hh>
	#include <scheme/nest/pmap/OriTransMap.hh>
	#include <scheme/objective/ObjectiveFunction.hh>
	#include <scheme/objective/hash/XformMap.hh>
	#include <scheme/objective/storage/RotamerScores.hh>
	#include <scheme/objective/voxel/FieldCache.hh>
	// #include <scheme/objective/voxel/VoxelArray.hh>
	#include <scheme/rosetta/score/RosettaField.hh>
	#include <scheme/util/StoragePolicy.hh>

	#include <utility/file/file_sys_util.hh>
	#include <utility/io/izstream.hh>
	#include <utility/io/ozstream.hh>
	#include <utility/tools/make_vector1.hh>


	#include <boost/random/mersenne_twister.hpp>
	#include <boost/random/uniform_real.hpp>

namespace devel {
namespace scheme {
namespace rif {


struct HBJob {
	std::string don, acc, don_or_acc;
	int nrots;
	int ires;
	bool operator<( HBJob const & other ) const { return nrots > other.nrots; }
};


	void
	RifGeneratorSimpleHbonds::generate_rif(
		RifAccumulatorP accumulator,
		RifGenParamsP params
	){

		core::pose::Pose const & target = *params->target;
		std::string const & target_tag = params->target_tag;
		utility::vector1<int> target_res = params->target_res;
		shared_ptr<RotamerIndex const> rot_index_p = params->rot_index_p;
		std::vector<std::string> const & cache_data_path = params->cache_data_path;
		std::vector< VoxelArray* > & field_by_atype = params->field_by_atype;

		// typedef Eigen::Transform<float,3,Eigen::AffineCompact> EigenXform;

		using core::id::AtomID;
		using std::cout;
		using std::endl;
		using namespace devel::scheme;
		typedef numeric::xyzVector<core::Real> Vec;
		typedef numeric::xyzVector<float> Vecf;
		typedef numeric::xyzMatrix<core::Real> Mat;
		typedef numeric::xyzTransform<core::Real> Xform;
		using ObjexxFCL::format::F;
		using ObjexxFCL::format::I;
		using devel::scheme::rif::RelRotPos;

		omp_lock_t cout_lock, io_lock, pose_lock, hacky_rpms_lock, hbond_geoms_cache_lock;
		omp_init_lock( &cout_lock );
		omp_init_lock( &io_lock );
		omp_init_lock( &pose_lock );
		omp_init_lock( &hbond_geoms_cache_lock );

		for( auto dir : cache_data_path ){
			std::cout << "RifGeneratorSimpleHbonds, cache data path: " << dir << std::endl;
		}


		RotamerIndex const & rot_index( *rot_index_p );

		std::vector<boost::random::mt19937> rngs;
		for( int i = 0; i < omp_max_threads_1(); ++i ){
			rngs.push_back( boost::random::mt19937( (unsigned int)time(0) + i) );
		}
		boost::uniform_real<> uniform;

		std::cout << "========================== DO HBOND ==================================" << std::endl;

		utility::vector1<std::string> donresn, donresn_std;
		utility::vector1<std::string> accresn, accresn_std;
		// utility::vector1<core::Size> target_res;
		{
			donresn_std.push_back( "ARG" );
			donresn_std.push_back( "ASN" );
			donresn_std.push_back( "GLN" );
			donresn_std.push_back( "HIS" );
			donresn_std.push_back( "HIS_D" );
			donresn_std.push_back( "LYS" );
			donresn_std.push_back( "SER" );
			donresn_std.push_back( "THR" );
			donresn_std.push_back( "TRP" );
			donresn_std.push_back( "TYR" );
            donresn_std.push_back( "ADE" );
            donresn_std.push_back( "CYT" );
            donresn_std.push_back( "GUA" );
            donresn_std.push_back( "THY" );
            donresn_std.push_back( "RAD" );
            donresn_std.push_back( "RCY" );
            donresn_std.push_back( "RGU" );
            donresn_std.push_back( "URA" );

			accresn_std.push_back( "ASN" );
			accresn_std.push_back( "ASP" );
			accresn_std.push_back( "GLN" );
			accresn_std.push_back( "GLU" );
			accresn_std.push_back( "HIS" );
			accresn_std.push_back( "HIS_D" );
			accresn_std.push_back( "SER" );
			accresn_std.push_back( "THR" );
			accresn_std.push_back( "TYR" );
            accresn_std.push_back( "ADE" );
            accresn_std.push_back( "CYT" );
            accresn_std.push_back( "GUA" );
            accresn_std.push_back( "THY" );
            accresn_std.push_back( "RAD" );
            accresn_std.push_back( "RCY" );
            accresn_std.push_back( "RGU" );
            accresn_std.push_back( "URA" );

			// if( donresn_user.size()==0 ) donresn_user = donresn_std;
			// if( accresn_user.size()==0 ) accresn_user = accresn_std;
			for( auto s : donresn_user ){
				if( std::find( donresn_std.begin(), donresn_std.end(), s ) == donresn_std.end() ){
					donresn_std.push_back(s);
					std::cout << "adding don res " << s << " to std set" << std::endl;
				}
			}
			for( auto s : accresn_user ){
				if( std::find( accresn_std.begin(), accresn_std.end(), s ) == accresn_std.end() ){
					accresn_std.push_back(s);
					std::cout << "adding acc res " << s << " to std set" << std::endl;
				}
			}
		}
		std::vector< ::scheme::chemical::HBondRay > target_donors, target_acceptors;
		for( auto ir : target_res ){
			::devel::scheme::get_donor_rays   ( target, ir, params->hbopt, target_donors );
			::devel::scheme::get_acceptor_rays( target, ir, params->hbopt, target_acceptors );
		}
		{
			if( target_donors.size() ){
				utility::io::ozstream donout(params->output_prefix+"donors.pdb.gz");
				::devel::scheme::dump_hbond_rays( donout, target_donors, true );
				donout.close();
			}
			if( target_acceptors.size() ){
				utility::io::ozstream accout(params->output_prefix+"acceptors.pdb.gz");
				::devel::scheme::dump_hbond_rays( accout, target_acceptors, false );
				accout.close();
			}
			// utility_exit_with_message("debug hbonds");
		}
		std::cout << "target_donors.size() " << target_donors.size() << " target_acceptors.size() " << target_acceptors.size() << std::endl;
		int n_sat_groups = 0;
		if( accumulator->rif()->has_sat_data_slots() ) n_sat_groups = target_donors.size() + target_acceptors.size();
		std::vector<utility::io::ozstream*> rif_hbond_vis_out_satgroups(n_sat_groups,nullptr);
		std::vector<std::vector<utility::io::ozstream*>> rif_hbond_vis_out_double_satgroups(
		                                                     n_sat_groups, std::vector<utility::io::ozstream*>(n_sat_groups,nullptr));
		std::cout << "n_sat_groups = " << n_sat_groups << std::endl;

	// target.dump_pdb("test.pdb");
	// utility_exit_with_message("DEBUG EXTRA ACCEPTORS!!!");
		devel::scheme::ScoreRotamerVsTarget<
				VoxelArrayPtr, ::scheme::chemical::HBondRay, ::devel::scheme::RotamerIndex
			> rot_tgt_scorer;
		{
			rot_tgt_scorer.rot_index_p_ = rot_index_p;
			rot_tgt_scorer.target_field_by_atype_ = field_by_atype;
			rot_tgt_scorer.target_donors_ = target_donors;
			rot_tgt_scorer.target_acceptors_ = target_acceptors;
			rot_tgt_scorer.hbond_weight_ = opts.hbond_weight;
			rot_tgt_scorer.upweight_multi_hbond_ = opts.upweight_multi_hbond;
			rot_tgt_scorer.upweight_iface_ = 1.0;

		}

		// int npairs = donresn.size()*accresn_user.size()  +  accresn.size()*donresn_user.size() ;

		// build up a joblist for hbond rif gen
		std::vector<HBJob> hb_jobs;
		core::chemical::ResidueTypeSetCAP rts = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");

		for(int ires = 1; ires <= target_res.size(); ++ires){
			int const ir = target_res[ires];
			std::string resn = target.residue(ir).name();
			std::cout << "RifGenSimpleHbonds checking res " << resn << std::endl;
			HBJob j;
			j.ires = ir;
			for( int iacc = 1; iacc <= accresn_user.size(); ++iacc ){
				core::chemical::ResidueType const & rtype = rts.lock()->name_map(accresn_user[iacc]);
				if( !rtype.has("N") || !rtype.has("CA") || !rtype.has("C") ){
					std::cout << "not putting " << accresn_user[iacc] << " into rif, no N,CA,C" << std::endl;
					continue;
				}
				j.don = "GLY";
				j.acc = accresn_user[iacc];
				if( std::find( accresn_std.begin(), accresn_std.end(), j.acc ) == accresn_std.end() ) continue; // no non-standard res in RIF
				j.don_or_acc = "DON_";
				std::pair<size_t,size_t> b = rot_index.index_bounds(j.acc.substr(0,3));
				j.nrots = b.second-b.first;
				//if( target.residue(ir).is_protein() && target.residue(ir).has("H") ) hb_jobs.push_back( j );
				if( target.residue(ir).has("H") ) hb_jobs.push_back( j );
				if( std::find(donresn_std.begin(),donresn_std.end(),resn)!=donresn_std.end() ){ // is donor
					j.don = resn;
					hb_jobs.push_back( j );
				}
			}
			for( int idon = 1; idon <= donresn_user.size(); ++idon ){
				core::chemical::ResidueType const & rtype = rts.lock()->name_map(donresn_user[idon]);
				if( !rtype.has("N") || !rtype.has("CA") || !rtype.has("C") ){
					std::cout << "not putting " << donresn_user[idon] << " into rif, no N,CA,C" << std::endl;
					continue;
				}
				j.don = donresn_user[idon];
				if( std::find( donresn_std.begin(), donresn_std.end(), j.don ) == donresn_std.end() ) continue; // no non-standard res in RIF
				j.acc = "GLY";
				j.don_or_acc = "ACC_";
				std::pair<size_t,size_t> b = rot_index.index_bounds(j.don.substr(0,3));
				j.nrots = b.second-b.first;
				//if( target.residue(ir).is_protein() && target.residue(ir).has("O") ) hb_jobs.push_back( j );
				if( target.residue(ir).has("O") ) hb_jobs.push_back( j );
				if( std::find(accresn_std.begin(),accresn_std.end(),resn)!=accresn_std.end() ){ // is acceptor
					j.acc = resn;
					hb_jobs.push_back( j );
				}
			}
		}
		runtime_assert_msg( hb_jobs.size() , "no hbond jobs generated!" );
		std::sort( hb_jobs.begin(), hb_jobs.end() );
		// for( int i = 0; i < hb_jobs.size(); ++i ){
		// 	HBJob j = hb_jobs[i];
		// 	std::cout << "HBJob " << I(4,i) << I(3,j.ires) << " " << target.residue(j.ires).name() << " "
		// 	          << j.don_or_acc << " " << j.don << " " << j.acc << " " << j.nrots << std::endl;
		// }
		// utility_exit_with_message("check HBJob LIST");

		std::map< std::string, utility::vector1< RelRotPos > * > hbond_geoms_cache;
		std::map< std::string, omp_lock_t > hbond_io_locks;
		for( int ihbjob = 0; ihbjob < hb_jobs.size(); ++ihbjob ){
			std::string don = hb_jobs[ihbjob].don;
			std::string acc = hb_jobs[ihbjob].acc;
			std::string don_or_acc =  hb_jobs[ihbjob].don_or_acc;
			int ir =  hb_jobs[ihbjob].ires;
			std::string hbgeomtag = don_or_acc + don + "-" + acc;
			bool override = target.size()==1 && ( target.residue(1).name()==don || target.residue(1).name()==acc );
			override |= !target.residue(ir).is_protein();
			if( override ) hbgeomtag += "__" + target_tag;
			hbond_geoms_cache[ hbgeomtag ] = nullptr;
			omp_lock_t tmplock;
			hbond_io_locks[ hbgeomtag ] = tmplock;
			omp_init_lock( & hbond_io_locks[ hbgeomtag ] );
		}

		// make sure hbond geometries exist and are loaded
		std::cout << "preloading / generating hbond_geometry files" << std::endl;
		#ifdef USE_OPENMP
		#pragma omp parallel for schedule(dynamic,1)
		#endif
		for( int ihbjob = 0; ihbjob < hb_jobs.size(); ++ihbjob ){
			std::string don = hb_jobs[ihbjob].don;
			std::string acc = hb_jobs[ihbjob].acc;
			std::string don_or_acc =  hb_jobs[ihbjob].don_or_acc;
			int nrots = hb_jobs[ihbjob].nrots;
			int ir =  hb_jobs[ihbjob].ires;
			std::string hbgeomtag = don_or_acc + don + "-" + acc;

            // std::cout << hbgeomtag << std::endl;
            // continue;

			std::map< std::string, core::pose::Pose > hbgeom_exemplars_rtype_override;
			bool override = target.size()==1 && ( target.residue(1).name()==don || target.residue(1).name()==acc );
			override |= !target.residue(ir).is_protein();
			if( override ){
				std::cout << "add " << target.residue(ir).name() << " to hbgeom_exemplars_rtype_override" << std::endl;
				omp_set_lock(&pose_lock);
				core::pose::Pose tmp0;
				tmp0.append_residue_by_jump(target.residue(ir),1);
				hbgeom_exemplars_rtype_override[ target.residue(ir).name() ] = tmp0;
				auto & tmp = hbgeom_exemplars_rtype_override[ target.residue(ir).name() ];
				runtime_assert( tmp.size() == 1 );
				if( tmp.residue(1).is_lower_terminus() ) core::pose::remove_lower_terminus_type_from_pose_residue( tmp, 1 );
				if( tmp.residue(1).is_upper_terminus() ) core::pose::remove_upper_terminus_type_from_pose_residue( tmp, 1 );
				using core::chemical::VIRTUAL_DNA_PHOSPHATE;
				if( tmp.residue(1).has_variant_type(VIRTUAL_DNA_PHOSPHATE) ){
					core::pose::remove_variant_type_from_pose_residue(tmp, VIRTUAL_DNA_PHOSPHATE, 1);
				}
				omp_unset_lock(&pose_lock);
				hbgeomtag += "__" + target_tag;
			}

			// seems init of hbond_geoms_cache is sloppy... check if exists and init if not
			if( hbond_geoms_cache.find(hbgeomtag) == hbond_geoms_cache.end() ){
				utility_exit_with_message("hbond_geoms_cache key not found "  +hbgeomtag);
			}

			// load hbond geom data if needed
			if( ! hbond_geoms_cache[hbgeomtag] ){

				// likewise, init of locks is sloppy... check and init lock if necessary
				if( hbond_io_locks.find(hbgeomtag) == hbond_io_locks.end() ){
					#ifdef USE_OPENMP
					#pragma omp critical
					#endif
					omp_init_lock( &hbond_io_locks[hbgeomtag] ); // inserts the lock and then inits
					// std::cout << "hbond_io_locks missing " << hbgeomtag << std::endl;
					// utility_exit_with_message("hbond_io_locks not initialized properly");
				}

				std::string cachefile = "__HBOND_GEOMS";
					cachefile += "__maxtip" + boost::lexical_cast<std::string>( opts.tip_tol_deg    ) ;
					cachefile += "__resl"   + boost::lexical_cast<std::string>( opts.rot_samp_resl  ) ;
					cachefile += "__range"  + boost::lexical_cast<std::string>( opts.rot_samp_range ) ;
					cachefile += "__ex1_0";
					cachefile += "__ex2_0";
					cachefile += "__ex3_0";
					cachefile += "__ex4_0";
					cachefile += "__nrots"  + boost::lexical_cast<std::string>( nrots ) ;
					cachefile += "__" + hbgeomtag + ".rel_rot_pos.gz";

				auto cache = new utility::vector1< RelRotPos >;
				omp_set_lock( &hbond_io_locks[hbgeomtag] ); // make sure nobody tries to use this while filling in...
				// runtime_assert( hbond_geoms_cache.find(hbgeomtag) != hbond_geoms_cache.end() )
				hbond_geoms_cache[hbgeomtag] = cache;
				omp_unset_lock( & hbond_geoms_cache_lock );

				bool failed_to_read = true;
				utility::io::izstream instream;
				std::string cachefile_found = devel::scheme::open_for_read_on_path( cache_data_path, cachefile, instream );
				if( cachefile_found.size() ){

					if( ihbjob==0 ){
						omp_set_lock(&cout_lock);
						cout << "load hbgeom " << cachefile_found << endl;
						cout << "            (will not log rest)" << endl;
						omp_unset_lock(&cout_lock);
					} else {
						std::cout << "*"; std::cout.flush();
					}

					size_t n;
					runtime_assert( instream.good() );
					instream.read( (char*)(&n), sizeof(size_t) );
					cache->resize( n );
					for(size_t i = 0; i < n; ++i){
						if( !instream.good() ) break;
						RelRotPos r;
						instream.read( (char*)(&r), sizeof(RelRotPos) );
						cache->at(i+1) = r;
					}
					instream.close();
					runtime_assert( instream.good() );
					failed_to_read = cache->size() != n;
				}

				if( failed_to_read ){

					utility::vector1< RelRotPos > & hbond_geoms( *hbond_geoms_cache[hbgeomtag] );

					omp_set_lock(&cout_lock);
					cout << "GENERATING HBOND GEOMETRIES pair " << ihbjob+1 << " of " << hb_jobs.size()
					     << " : " << don << "/" << acc << " -- FIX_" << don_or_acc << endl;
					omp_unset_lock(&cout_lock);

					devel::scheme::rif::MakeHbondGeomOpts mhbopts;
					mhbopts.tip_tol_deg    = opts.tip_tol_deg;
					mhbopts.rot_samp_resl  = opts.rot_samp_resl;
					mhbopts.rot_samp_range = opts.rot_samp_range;
					devel::scheme::rif::make_hbond_geometries(
						rot_index,
						don,
						acc,
						don_or_acc=="DON_",
						don_or_acc=="ACC_",
						hbgeom_exemplars_rtype_override,
						hbond_geoms,
						mhbopts
					);

					// fixed now... cachefile will have target_tag iff any exemplar
					// if( hbgeom_exemplars_rtype_override.size() != 0 ){
						// std::cout << "WARNING: storing exemplar to cache!!!" << std::endl;
					// }

					size_t n = hbond_geoms.size();
					utility::io::ozstream out;
					std::string cachefile_found = devel::scheme::open_for_write_on_path( cache_data_path, cachefile, out, true );
					if( cachefile_found.size() ){
						omp_set_lock(&cout_lock);
							cout << "SAVING " << KMGT(n) << " HBOND GEOMETRIES TO " << cachefile_found << endl;
						omp_unset_lock(&cout_lock);
						runtime_assert( out.good() );
						out.write( (char*)(&n), sizeof(size_t) );
						cout << n << endl;
						for(size_t i = 0; i < n; ++i){
							out.write( (char*)( &hbond_geoms[i+1] ), sizeof(RelRotPos) );
						}
						out.close();
					} else {
						std::cout << "WARNING: can't save HBOND GEOMETRIES for " << cachefile << ", they will be regenerated every time!" << std::endl;
					}
				}

				omp_unset_lock( &hbond_io_locks[hbgeomtag] ); // now is ready


			}
			if( ! hbond_geoms_cache[hbgeomtag] ){
				utility_exit_with_message( "hbond_geoms_cache missing for " + hbgeomtag );
			}
		}
		std::cout << endl;

		for( int ihbjob = 0; ihbjob < hb_jobs.size(); ++ihbjob ){
			std::string don = hb_jobs[ihbjob].don;
			std::string acc = hb_jobs[ihbjob].acc;
			std::string don_or_acc =  hb_jobs[ihbjob].don_or_acc;
			int nrots = hb_jobs[ihbjob].nrots;
			int ir =  hb_jobs[ihbjob].ires;
			std::string hbgeomtag = don_or_acc + don + "-" + acc;
			bool override = target.size()==1 && ( target.residue(1).name()==don || target.residue(1).name()==acc );
			override |= !target.residue(ir).is_protein();
			if( override ) hbgeomtag += "__" + target_tag;

			utility::io::ozstream *rif_hbond_vis_out = nullptr;

			omp_set_lock( &hbond_io_locks[hbgeomtag] );
			if( ! hbond_geoms_cache[hbgeomtag] ){
				utility_exit_with_message( "hbond_geoms_cache missing for " + hbgeomtag );
			}
			utility::vector1< RelRotPos > const & hbond_geoms( *hbond_geoms_cache[hbgeomtag] );
			omp_unset_lock( &hbond_io_locks[hbgeomtag] );

			// loop over hbond geometries, then loop over residues which might have those geoms
			// do in this order to reduce the numker of geom datasets in memory at once
			// std::vector< std::vector< std::pair<uint64_t,std::pair<float,int32_t> > >  > to_insert( omp_max_threads_1() );
			// std::vector< std::pair<uint64_t,std::pair<float,int32_t> > > to_insert;


			std::string resn = target.residue(ir).name();
			if( don_or_acc == "DON_" ) assert( don == resn || don == "GLY" );
			if( don_or_acc == "ACC_" ) assert( acc == resn || acc == "GLY" );

			std::string anchor_atom1 = target.residue(ir).atom_name( target.residue(ir).nheavyatoms()-2 );
			std::string anchor_atom2 = target.residue(ir).atom_name( target.residue(ir).nheavyatoms()-1 );
			std::string anchor_atom3 = target.residue(ir).atom_name( target.residue(ir).nheavyatoms()-0 );
			if( resn == "SER" ){
				anchor_atom1 = "HG"; // special case for ser because so small, use HG, CB, OG
			}
			if( don_or_acc == "DON_" && don == "GLY" ){
				runtime_assert( acc != "PRO" );
				anchor_atom1 = "H";
				anchor_atom2 = "N";
				anchor_atom3 = "CA";
			}
			if( don_or_acc == "ACC_" && acc == "GLY" ){
				runtime_assert( don != "PRO" );
				anchor_atom1 = "O";
				anchor_atom2 = "C";
				anchor_atom3 = "CA";
			}
			if( don_or_acc == "ACC_" && ( acc=="ADX" || acc=="CYX" || acc=="GUX" || acc=="THX" ) ){
				anchor_atom1 = "OP1";
				anchor_atom2 = "P";
				anchor_atom3 = "OP2";
			}

			// std::cout << "HBOND ALIGN STUB " << anchor_atom1 << " " << anchor_atom2 << " " << anchor_atom3 << std::endl;
			Xform target_frame( target.xyz(AtomID(target.residue(ir).atom_index(anchor_atom1),ir)),
			                    target.xyz(AtomID(target.residue(ir).atom_index(anchor_atom2),ir)),
			                    target.xyz(AtomID(target.residue(ir).atom_index(anchor_atom3),ir)) );

			// #pragma omp critical
			// std::cout << std::endl << "FRAME ALIGN: " << AtomID(target.residue(ir).atom_index(anchor_atom1),ir) << " "
			// 							 << AtomID(target.residue(ir).atom_index(anchor_atom2),ir) << " "
			// 							 << AtomID(target.residue(ir).atom_index(anchor_atom3),ir) << "   "
			// 							 << anchor_atom1 << " " << target.residue(ir).xyz( anchor_atom1 ) << "   "
			// 							 << anchor_atom2 << " " << target.residue(ir).xyz( anchor_atom2 ) << "   "
			// 							 << anchor_atom3 << " " << target.residue(ir).xyz( anchor_atom3 ) << "   "
			// 							 << std::endl
			// 							 << std::endl;

			// omp_set_lock(&cout_lock);
			// std::cout << "    Thread " << I(3,omp_get_thread_num()) << " checking residue " << I(4,ir) << " for " << hbgeomtag << ", N= "
			// 			 << KMGT(hbond_geoms.size()) << "     " << I(4,ihbjob) << " of " << hb_jobs.size() <<  std::endl;
			// omp_unset_lock(&cout_lock);

			int rrpcount = 0;
			float range = opts.hbond_cart_sample_hack_range;
			// if( hbond_geoms.size() > 2*1000*1000 ) range /= 2.0; // hack to reduce non-ideal sampling for high-count hbonds
			// if( acc=="TYR" || don=="TYR" ) range *= 1.5; // TYR is long and rigid....
			float const range_nsamp = std::max( 0.0001f, std::ceil( range / (float)opts.hbond_cart_sample_hack_resl ) );
			if( range == 0.0 ) range = 0.001;

			std::exception_ptr exception = nullptr;
			#ifdef USE_OPENMP
			#pragma omp parallel for schedule(dynamic,64)
			#endif
			for( int i = 1; i <= hbond_geoms.size(); ++i ){
				if(exception) continue;
				try {
					RelRotPos const & rel_rot_pos( hbond_geoms[i] );
					// if( ++rrpcount%100000==0 ){ cout << (float)rrpcount/hbond_geoms.size()*100.0 << "%(t"<<omp_thread_num_1()<<") "; cout.flush(); }

					int irot = rel_rot_pos.rotamer;
					float hb_pair_score = rel_rot_pos.score;
					int hbonding_atom_num = rel_rot_pos.hbonding_atom_num;
					std::vector<SchemeAtom> const & res_atoms( rot_index.atoms(irot) );

					Xform const res_atomsbbf( Vec(res_atoms[0].position()[0],res_atoms[0].position()[1],res_atoms[0].position()[2]) ,
					                          Vec(res_atoms[1].position()[0],res_atoms[1].position()[1],res_atoms[1].position()[2]) ,
					                          Vec(res_atoms[2].position()[0],res_atoms[2].position()[1],res_atoms[2].position()[2]) );
					if( rel_rot_pos.n .distance_squared( rel_rot_pos.ca ) < 1.0 ||
					    rel_rot_pos.ca.distance_squared( rel_rot_pos.c  ) < 1.0 ||
					    rel_rot_pos.c .distance_squared( rel_rot_pos.n  ) < 1.0  )
					{
						std::cout << "bad data in rel_rot_pos, sourced from: " << hbgeomtag << std::endl;
						continue;
					}
					Xform const hbondbbf     ( rel_rot_pos.n    , rel_rot_pos.ca   , rel_rot_pos.c     );
					Xform const hbonder_frame( rel_rot_pos.stub1, rel_rot_pos.stub2, rel_rot_pos.stub3 );

					// ~res_atomsbbf to rotamer local bb coords
					//	hbondbbf to hbonder bb frame
					Xform const xalign = ~hbonder_frame * target_frame * hbondbbf * ~res_atomsbbf;

					// std::cout << ~hbonder_frame * target_frame << endl;
					Vec hbpos1(
						res_atoms[rel_rot_pos.hbonding_atom_num].position()[0],
						res_atoms[rel_rot_pos.hbonding_atom_num].position()[1],
						res_atoms[rel_rot_pos.hbonding_atom_num].position()[2]
					);
					Vec hbpos2 = target.residue(ir).xyz( rel_rot_pos.hbonding_atom_num_other+1 );
					if( opts.debug && fabs( hbpos2.distance(xalign*hbpos1) - 2.8 ) > 0.5 ) {
						{
							omp_set_lock(&cout_lock);
								cout << hbpos2.distance(xalign*hbpos1) << " "
									 << rel_rot_pos.hbonding_atom_num_other << " "
									 << target.residue(ir).atom_name(rel_rot_pos.hbonding_atom_num_other+1)  << " "
									 << rot_index.rotamers_[irot].resname_  << " "
									 << rot_index.chem_index_.atom_data( rot_index.rotamers_[irot].resname_ , rel_rot_pos.hbonding_atom_num ).atomname
									 << endl;
							omp_unset_lock(&cout_lock);
							omp_set_lock(&io_lock);
								utility::io::ozstream out("test.pdb");
								for( auto a : res_atoms ){
									Vec tmp( a.position()[0], a.position()[1], a.position()[2] );
									tmp = xalign*tmp;
									a.set_position( tmp );
									::scheme::actor::write_pdb( out, a, rot_index.chem_index_ );
								}
								out.close();
							omp_unset_lock(&io_lock);
							utility_exit_with_message("debug hbond distance check");
						}
					}


					for( float dx = -range; dx <= range+0.000001; dx += range/range_nsamp ){
					for( float dy = -range; dy <= range+0.000001; dy += range/range_nsamp ){
					for( float dz = -range; dz <= range+0.000001; dz += range/range_nsamp ){

						// if( 0.375 < fabs( 2.7 - hbpos2.distance( xalign * ( hbpos1+Vec(dx,dy,dz) ) ) ) ) continue;
						using devel::scheme::score_hbond_rays;

						Vec const N  = xalign * Vec(res_atoms[0].position()[0]+dx,res_atoms[0].position()[1]+dy,res_atoms[0].position()[2]+dz);
						Vec const CA = xalign * Vec(res_atoms[1].position()[0]+dx,res_atoms[1].position()[1]+dy,res_atoms[1].position()[2]+dz);
						Vec const C  = xalign * Vec(res_atoms[2].position()[0]+dx,res_atoms[2].position()[1]+dy,res_atoms[2].position()[2]+dz);
						::scheme::actor::BackboneActor<EigenXform> bbactor( N, CA, C );

						// float positioned_rotamer_score = dx*dx+dy*dy+dz*dz;
						// for(int ia = 0; ia < res_atoms.size(); ++ia){
						// 	int at = res_atoms[ia].type();
						// 	if( ia == rel_rot_pos.hbonding_atom_num || at > 21 ) continue; // skip hbonding atom for fa-rep issues && hydrogen
						// 	// #pragma omp critical
						// 	// if( ia == hbonding_atom_num ) cout << ia << " " << res_atoms[ia].type() << " " << res_atoms[ia].data().atomname << endl;
						// 	// if( ia > 2 && at > 6 ) continue;
						// 	Vec v( res_atoms[ia].position()[0]+dx, res_atoms[ia].position()[1]+dy, res_atoms[ia].position()[2]+dz );
						// 	v = xalign * v;
						// 	positioned_rotamer_score += field_by_atype[at]->at( v[0], v[1], v[2] );
						// }
						// positioned_rotamer_score += hb_pair_score;


						int sat1=-1, sat2=-1;
						float positioned_rotamer_score = rot_tgt_scorer.score_rotamer_v_target_sat( irot, bbactor.position_, sat1, sat2, 10.0, 0 );
						if( positioned_rotamer_score > opts.score_threshold ) continue;

						if( n_sat_groups > 0 ){
							runtime_assert( sat1 < n_sat_groups && sat2 < n_sat_groups );
							// if( sat1 < 0 || sat2 >= 0 ){
							// 	#pragma omp critical
							// 	{
							// 		std::cout << "bad_sat score: " << positioned_rotamer_score << " "
							//             << opts.score_threshold << " " << sat1 << " " << sat2 << std::endl;
							// 		utility::io::ozstream out("bad_sat.pdb");
							// 		for( auto a : res_atoms ){
							// 			Vec tmp( a.position()[0]+dx, a.position()[1]+dy, a.position()[2]+dz );
							// 			tmp = xalign*tmp;
							// 			a.set_position( tmp );
							// 			::scheme::actor::write_pdb( out, a, rot_index.chem_index_ );
							// 		}
							// 		out.close();
							// 		utility_exit_with_message("why is sat1 < 0????");
							// 	}
							// }
						}

						accumulator->insert( bbactor.position_, positioned_rotamer_score, irot, sat1, sat2 );


						// // shitty test output
						if( opts.dump_fraction > 0 ){
							double const runif = uniform(rngs[omp_thread_num_1()-1]);
							float dump_chance = opts.dump_fraction;
							// if( positioned_rotamer_score < -2.0 ) dump_chance *= 2.0;
							// if( positioned_rotamer_score < -3.0 ) dump_chance *= 2.0;
							// if( positioned_rotamer_score < -4.0 ) dump_chance *= 2.0;
							// if( positioned_rotamer_score < -5.0 ) dump_chance *= 2.0;
							// if( positioned_rotamer_score < -6.0 ) dump_chance *= 2.0;
							// if( n_sat_groups > 0 && sat1 >= 0 && sat2 >= 0 ) dump_chance *= 10.0;
							if( runif < dump_chance ){
								// std::cout << "do test output" << std::endl;
								omp_set_lock(&io_lock);
									if( rif_hbond_vis_out == nullptr ){
										std::string outfilename = params->output_prefix+"RifGen_Hbond_vis_"+I(3,ir)+hbgeomtag+".pdb.gz";
										// std::cout << "init1 " << outfilename << " " << runif << " " << opts.dump_fraction << std::endl;
										rif_hbond_vis_out = new utility::io::ozstream( outfilename );
									}
									*rif_hbond_vis_out << "MODEL " << irot << endl;
									for( auto a : res_atoms ){
										Vec tmp( a.position()[0]+dx, a.position()[1]+dy, a.position()[2]+dz );
										tmp = xalign*tmp;
										a.set_position( tmp );
										::scheme::actor::write_pdb( *rif_hbond_vis_out, a, rot_index.chem_index_ );
									}
									*rif_hbond_vis_out << "ENDMDL" << endl;
									if( n_sat_groups > 0 ){
										if( sat1 >= 0 && sat2 >= 0 ){
											if( sat1 < sat2 ) std::swap(sat1,sat2);
											if( rif_hbond_vis_out_double_satgroups[sat1][sat2] == nullptr ){
												std::string outfilename = params->output_prefix+"RifGen_Hbond_vis_doublesat"+I(4,sat1)+I(4,sat2)+".pdb.gz";
												// std::cout << "init1 " << outfilename << " " << runif << " " << opts.dump_fraction << std::endl;
												rif_hbond_vis_out_double_satgroups[sat1][sat2] = new utility::io::ozstream( outfilename );
											}
											*rif_hbond_vis_out_double_satgroups[sat1][sat2] << "MODEL " << irot << endl;
											for( auto a: res_atoms ){
												Vec tmp( a.position()[0]+dx, a.position()[1]+dy, a.position()[2]+dz );
												tmp = xalign*tmp;
												a.set_position( tmp );
												::scheme::actor::write_pdb( *rif_hbond_vis_out_double_satgroups[sat1][sat2], a, rot_index.chem_index_ );
											}
											*rif_hbond_vis_out_double_satgroups[sat1][sat2] << "ENDMDL" << endl;
										} else {
											for( int i12 = 0; i12 < 2; ++i12 ){
												int sat = i12 ? sat1 : sat2;
												if( sat >= 0 ){
													runtime_assert( sat < n_sat_groups );
													if( rif_hbond_vis_out_satgroups[sat] == nullptr ){
														std::string outfilename = params->output_prefix+"RifGen_Hbond_vis_sat"+I(4,sat)+".pdb.gz";
														// std::cout << "init1 " << outfilename << " " << runif << " " << opts.dump_fraction << std::endl;
														rif_hbond_vis_out_satgroups[sat] = new utility::io::ozstream( outfilename );
													}
													*rif_hbond_vis_out_satgroups[sat] << "MODEL " << irot << endl;
													for( auto a: res_atoms ){
														Vec tmp( a.position()[0]+dx, a.position()[1]+dy, a.position()[2]+dz );
														tmp = xalign*tmp;
														a.set_position( tmp );
														::scheme::actor::write_pdb( *rif_hbond_vis_out_satgroups[sat], a, rot_index.chem_index_ );
													}
													*rif_hbond_vis_out_satgroups[sat] << "ENDMDL" << endl;
												}
											}
										}
									}
								omp_unset_lock(&io_lock);
							}
						}
						// for( int itest = 0; itest < test_bbs.size(); ++itest ){
						// 	float const rms2 = (( test_bbs[itest][0] - Eigen::Vector3f(N [0],N [1],N [2]) ).squaredNorm() +
						// 							  ( test_bbs[itest][1] - Eigen::Vector3f(CA[0],CA[1],CA[2]) ).squaredNorm() +
						// 							  ( test_bbs[itest][2] - Eigen::Vector3f(C [0],C [1],C [2]) ).squaredNorm() ) / 3.0;
						// 	// cout << sqrt(rms) << " "; cout.flush();
						// 	if( rms2 < test_rms2_cut ){
						// 		{
						// 			if( rif_hbond_vis_out == nullptr ){
						// 				std::cout << "init2 rif_hbond_vis_"+I(3,ir)+hbgeomtag+".pdb" << std::endl;
						// 				rif_hbond_vis_out = new utility::io::ozstream("rif_hbond_vis_"+I(3,ir)+hbgeomtag+".pdb");
						// 			}
						// 			omp_set_lock(&io_lock);
						// 			// std::cout << "CLOSE ROTAMER " << irot << " hbond_atom_num " << 
						//                rel_rot_pos.hbonding_atom_num << " score: " << rel_rot_pos.score << std::endl;
						// 			*rif_hbond_vis_out << "MODEL " << boost::lexical_cast<std::string>(key) << endl;
						// 			for( auto a: res_atoms ){
						// 				Vec tmp( a.position()[0]+dx, a.position()[1]+dy, a.position()[2]+dz );
						// 				tmp = xalign*tmp;
						// 				a.set_position( tmp );
						// 				::scheme::actor::write_pdb( *rif_hbond_vis_out, a, rot_index.chem_index_ );
						// 			}
						// 			*rif_hbond_vis_out << "ENDMDL" << endl;
						// 			// rif_hbond_vis_out.close();
						// 			// utility_exit_with_message("foo");
						// 			omp_unset_lock(&io_lock);
						// 		}
						// 	}
						// }

					}}}


				} catch( ... ) {
					#ifdef USE_OPENMP
					#pragma omp critical
					#endif
					exception = std::current_exception();
				}
			}
			if( exception ) std::rethrow_exception(exception);


			// omp_set_lock(&cout_lock);
			#ifdef USE_OPENMP
			#pragma omp critical
			#endif
			{
				std::cout << "RifGenSimpHB: " << I(4,ir) << " " << ObjexxFCL::format::LJ(20,hbgeomtag) << " done, ngeom: "
				          << KMGT(hbond_geoms.size()) << " job " << I(4,ihbjob) << " of " << hb_jobs.size() << "  ";
				accumulator->checkpoint( cout );
				if( ihbjob%1 == 0 ){
					accumulator->report( cout );
				}
			}
			// omp_unset_lock(&cout_lock);

			if( rif_hbond_vis_out ){
				rif_hbond_vis_out->close();
				delete rif_hbond_vis_out;
			}


		}

		// cout << "RIF so far: " << " non0 in RIF: " << KMGT(rif.size()) << " N_motifs_found: "
		//      << KMGT(N_motifs_found) << " coverage: " << (double)N_motifs_found/rif.size() << endl;


		// cleanup hbond_io_locks
		for( auto & i : hbond_io_locks ) omp_destroy_lock( &(i.second) );
		for( auto & i : hbond_geoms_cache ) delete i.second;


		for( auto ozp: rif_hbond_vis_out_satgroups ){
			if( ozp != nullptr ){
				ozp->close(); // example pdb out file
				delete ozp;
			}
		}
		for( auto vozp: rif_hbond_vis_out_double_satgroups ){
			for( auto ozp : vozp ){
				if( ozp != nullptr ){
					ozp->close(); // example pdb out file
					delete ozp;
				}
			}
		}


		omp_destroy_lock( & cout_lock ) ;
		omp_destroy_lock( & io_lock );
		omp_destroy_lock( & pose_lock );
		omp_destroy_lock( & hbond_geoms_cache_lock );

		//change me back to exit if want cactually check simple hbonds
		if( accumulator->n_motifs_found() == 0 ){
			utility_exit_with_message("no hbonds found, something is wrong");
			//change me back to exit if want cactually check simple hbonds
			
		}
	}


}
}
}

