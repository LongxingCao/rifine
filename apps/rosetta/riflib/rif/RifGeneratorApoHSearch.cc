// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.



#include <riflib/rif/RifGeneratorApoHSearch.hh>


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
	#include <core/scoring/motif/util.hh>

	#include <devel/init.hh>
	#include <riflib/RotamerGenerator.hh>
	#include <riflib/util.hh>

	#include <map>

	#include <scheme/actor/Atom.hh>
	#include <scheme/actor/BackboneActor.hh>
	#include <scheme/actor/VoxelActor.hh>
	#include <scheme/chemical/RotamerIndex.hh>
	#include <scheme/kinematics/Director.hh>
	#include <scheme/kinematics/Scene.hh>
	#include <scheme/nest/pmap/OriTransMap.hh>
	#include <scheme/objective/ObjectiveFunction.hh>
	// #include <scheme/objective/hash/XformMap.hh>
	// #include <scheme/objective/storage/RotamerScores.hh>
	#include <scheme/objective/voxel/FieldCache.hh>
	// #include <scheme/objective/voxel/VoxelArray.hh>
	#include <scheme/rosetta/score/RosettaField.hh>
	#include <scheme/util/StoragePolicy.hh>
	#include <scheme/chemical/stub.hh>
	#include <scheme/io/dump_pdb_atom.hh>

	#include <utility/file/file_sys_util.hh>
	#include <utility/io/izstream.hh>
	#include <utility/io/ozstream.hh>
	#include <utility/tools/make_vector1.hh>

	#include <parallel/algorithm>

	#include <boost/random/mersenne_twister.hpp>
	#include <boost/random/uniform_real.hpp>


namespace devel {
namespace scheme {
namespace rif {


	#pragma pack (push, 4) // allows size to be 12 rather than 16
	struct SearchPoint {
		float score;
		uint64_t index;
		SearchPoint() : score(9e9), index(0) {}
		SearchPoint(uint64_t i) : score(9e9), index(i) {}
		bool operator < (SearchPoint const & o) const {
			return score < o.score;
		}
	};
	#pragma pack (pop)
	struct CMPSCORE {
		bool operator()( SearchPoint const & a, SearchPoint const & b ) const {
			return a.score < b.score;
		}
	};
	struct CMPINDEX {
		bool operator()( SearchPoint const & a, SearchPoint const & b ) const {
			return a.index < b.index;
		}
	};

	template<class Scene, class Director, class RotamerIndex>
	void dump_scene(
		Director const & director,
		Scene & tmp_scene,
		RotamerIndex const & rot_index,
		uint64_t index,
		int resl,
		std::ostream & out
	){
		director.set_scene( index, resl, tmp_scene );
		out << "MODEL" << std::endl;
		write_pdb( out, tmp_scene, rot_index.chem_index_ );
		out << "ENDMDL" << std::endl;

	}


	void
	RifGeneratorApoHSearch::generate_rif(
		RifAccumulatorP accumulator,
		RifGenParamsP params
	){

		core::pose::Pose const & target = *params->target;
		std::string const & target_tag = params->target_tag;
		utility::vector1<int> target_res = params->target_res;
		shared_ptr<RotamerIndex const> rot_index_p = params->rot_index_p;
		std::vector<std::string> const & cache_data_path = params->cache_data_path;
		std::vector< VoxelArray* > & field_by_atype = params->field_by_atype;


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
		using devel::scheme::KMGT;
		using ::scheme::chemical::make_stub;
		using ::Eigen::Vector3f;

		typedef ::scheme::util::SimpleArray<3,float> F3;
		typedef ::scheme::util::SimpleArray<3,int> I3;

		// typedef ::scheme::objective::storage::RotamerScores< 12 > XMapVal;
		// typedef ::scheme::objective::hash::XformMap< EigenXform
		// 									, XMapVal
		// 									, ::scheme::objective::hash::XformHash_bt24_BCC6
		// 								> XMap;
		// cout << "size of hashed value: " << sizeof(typename XMap::Map::value_type) << endl;
		// load 1.25, 0-162m, 1-208m, 2-294m, 3-453m
		// load 3.67, 0-56, 1-69m


		typedef ::scheme::actor::      Atom< Eigen::Vector3f > SchemeAtom;
		typedef ::scheme::actor::SimpleAtom< Eigen::Vector3f > SimpleAtom;
		// runtime_assert( sizeof(SchemeAtom) == 16 );
		// typedef scheme::rosetta::score::RosettaField<SchemeAtom,devel::scheme::EtableParamsInit> RosettaField;

		typedef ::scheme::objective::voxel::VoxelArray<3,float> VoxelArray;
		typedef ::scheme::objective::voxel::BoundingFieldCache3D<float> BoundingGrid;



		typedef ::scheme::nest::NEST<
										  6,
										  EigenXform,
										  ::scheme::nest::pmap::OriTransMap,
										  ::scheme::util::StoreNothing, // do not store a transform in the Nest
										  uint64_t,
										  float,
										  false // do not inherit from NestBase
										 > Nest;

		typedef ::scheme::kinematics::NestDirector< Nest, uint64_t > Director;


		typedef SimpleAtom SceneAtom;
		typedef ::scheme::actor::VoxelActor<EigenXform,float> VoxelActor;

		typedef ::scheme::actor::Score_Voxel_vs_Atom< VoxelActor, SceneAtom > VoxelScore;
		typedef ::scheme::objective::ObjectiveFunction<
		boost::mpl::vector<
			VoxelScore
		>,
		int // Config type, just resl
		> Objective;


		typedef boost::mpl::vector<
			SceneAtom,
			VoxelActor
		> Actors;
		typedef ::scheme::kinematics::impl::Conformation<Actors> Conformation;
		typedef ::scheme::kinematics::Scene<Conformation,EigenXform> Scene;
		// typedef scheme::shared_ptr<scheme::kinematics::SceneBase<EigenXform> > SceneP;


		int64_t const DIM = 6;
		int64_t const DIMPOW2 = 1<<DIM;
		int64_t const beam_size = int64_t( opts.beam_size_M * 1000000.0 / DIMPOW2 ) * DIMPOW2;


		omp_lock_t cout_lock, io_lock, accum_lock;
		omp_init_lock( &cout_lock );
		omp_init_lock( &io_lock );
		omp_init_lock( &accum_lock );

		std::vector<boost::random::mt19937> rngs;
		for( int i = 0; i < omp_max_threads_1(); ++i ){
			rngs.push_back( boost::random::mt19937( (unsigned int)time(0) + i) );
		}
		boost::uniform_real<> uniform;


		std::vector<int> rots;
		for( auto resn : apores ){
			if( resn.size() > 3 ){
				int irot = boost::lexical_cast<int>(resn.substr(3,resn.size()-3));
				std::string threeletter = resn.substr(0,3);
				std::cerr << rot_index_p->rotamers_[irot].resname_ << std::endl;
				if( rot_index_p->rotamers_[irot].resname_ == threeletter ){
					if( rot_index_p->is_primary(irot) && rot_index_p->is_structural_primary(irot) ){
						rots.push_back( irot );
					} else {
						utility_exit_with_message(
							"rotamer number " + str(irot) + " is not primary rotamer of " + resn);
					}
				} else {
					utility_exit_with_message("residue name3 dosen't match: "+resn);
				}
			} else {
				std::pair<int,int> bounds = rot_index_p->index_bounds( resn );
				for( int irot = bounds.first; irot < bounds.second; ++irot ){
					if( rot_index_p->is_primary(irot) && rot_index_p->is_structural_primary(irot) ){
						rots.push_back( irot );
					}
				}
			}
		}

		std::map<std::string,float> abs_score_cut_by_res;
		abs_score_cut_by_res["ALA"] = -0.3;
		abs_score_cut_by_res["PHE"] = -2.8;
		abs_score_cut_by_res["ILE"] = -1.7;
		abs_score_cut_by_res["LEU"] = -1.5;
		abs_score_cut_by_res["MET"] = -1.7;
		abs_score_cut_by_res["VAL"] = -1.0;
		abs_score_cut_by_res["TRP"] = -3.4;
		abs_score_cut_by_res["LYS"] = -2.6;

		// #ifdef USE_OPENMP
		// #pragma omp parallel for schedule(dynamic,1)
		// #endif
		for( int ijob = 0; ijob < rots.size(); ++ijob ){

			int irot = rots[ijob];
			std::string resn = rot_index_p->rotamers_[irot].resname_;
			runtime_assert_msg( abs_score_cut_by_res.find(resn) != abs_score_cut_by_res.end(), "unsupported res "+resn );
			float const abs_score_cut_by_res_thisres = abs_score_cut_by_res[resn] * opts.score_cut_adjust;

			// utility::io::ozstream rif_apo_vis_out("rif_apo_vis_"+resn+str(irot)+".pdb");

			{
				omp_set_lock(&cout_lock);
				// cout << "========================================================================================================" << endl;
				cout << "================== ApoHSearch rotamer " << irot << " " << resn << " chis: ";
				for( int i = 0; i < rot_index_p->rotamers_[irot].chi_.size(); ++i ) cout << " " << rot_index_p->rotamers_[irot].chi_[i];
				cout << " ================== Progress: " << ijob+1 << " of " << rots.size() << " " << ijob*1.f/rots.size()*100.0f << "\% ==================" << endl;
				// cout << "========================================================================================================" << endl;
				omp_unset_lock(&cout_lock);
			}

			Scene scene_proto(2);
			float rotamer_radius = 0;
			Eigen::Vector3f rotamer_center(0,0,0);
			{
				int firstatom = 3; // CB
				if( rot_index_p->nchi(irot) == 0 ) firstatom = 0;
				if( rot_index_p->nchi(irot) == 1 ) firstatom = 1;
				// these things are not takin into account below when the N/CA/C score is added in
				// possible double repulsive for the bb atoms... probably no big deal
				int natom = 0;
				for( int ia = firstatom; ia < rot_index_p->nheavyatoms(irot); ++ia){
					if( ia == 2 && rot_index_p->nchi(irot) > 0 )  continue;
					rotamer_center += rot_index_p->atom(irot,ia).position();
					natom++;
				}
				rotamer_center /= float(natom);
				BOOST_FOREACH( SchemeAtom const & a, rot_index_p->rotamers_[irot].atoms_ )
					rotamer_radius = std::max( (a.position()-rotamer_center).norm(), rotamer_radius );
				for( int ia = firstatom; ia < rot_index_p->rotamers_[irot].atoms_.size(); ++ia){
					if( ia == 2 && rot_index_p->nchi(irot) > 0 )  continue;
					SchemeAtom const & a( rot_index_p->rotamers_[irot].atoms_[ia] );
					if( a.type() >= 21 ) continue;
					runtime_assert( rot_index_p->chem_index_.resname2num_.find(resn) != rot_index_p->chem_index_.resname2num_.end() );
					int restype = rot_index_p->chem_index_.resname2num_.find(resn)->second;
					int atype = a.type();
					if( atype > 6 ) { // polar sidechain or bb atom
						atype = 20; // Score_Voxel_vs_Atom will treat as replsive-only
					}
					SceneAtom sa( a.position()-rotamer_center, atype, restype, ia );
					runtime_assert( rot_index_p->chem_index_.atom_data( restype, ia ) == a.data() );
					scene_proto.add_actor(1,sa);
				}
				scene_proto.add_actor( 0, VoxelActor( bounding_by_atype ) );
			}


			struct RotChild {
				int rotid;
				Eigen::Vector3f Ncen, CAcen, Ccen;
			};
			std::vector<RotChild> inv_rotamer_backbones;

			std::cout << "RifGeneratorApoHSearch: add child rotamers:";
			for( size_t crot = 0; crot < rot_index_p->size(); ++crot ){
				if( rot_index_p->structural_parent_of_.at(crot) == irot && rot_index_p->is_primary(crot) ){
					RotChild child;
					child.rotid = crot;
					auto x = rot_index_p->to_structural_parent_frame_.at(crot);
					child.Ncen  = x * rot_index_p->atom(crot,0).position() - rotamer_center;
					child.CAcen = x * rot_index_p->atom(crot,1).position() - rotamer_center;
					child.Ccen  = x * rot_index_p->atom(crot,2).position() - rotamer_center;
					auto testca  = x * rot_index_p->atom(crot,3).position() - rotamer_center;
					auto primaryca = rot_index_p->atom(irot,3).position()-rotamer_center;
					runtime_assert( (primaryca - testca).norm() < 0.001 );
					inv_rotamer_backbones.push_back(child);
					std::cout << " " << resn << crot;
				}
			}
			std::cout << std::endl;

			float const half_tgt_resl = RESLS.front()/2.0;
			float rot_resl_deg;
				if( half_tgt_resl/2.0 > rotamer_radius ) rot_resl_deg = 180.0;
				else 			                         rot_resl_deg = 2.0 * 180.0 / M_PI * std::asin( half_tgt_resl/2.0/rotamer_radius );
			// rot_resl_deg = std::min( rot_resl_deg, 30.0f );

			F3 lb0( bounding_by_atype[0][1]->lb_ + rotamer_radius );
			F3 ub0( bounding_by_atype[0][1]->ub_ + rotamer_radius );
			I3 nc( std::ceil( (ub0[0]-lb0[0])/half_tgt_resl*sqrt(3.0)/2.0 ),
					 std::ceil( (ub0[1]-lb0[1])/half_tgt_resl*sqrt(3.0)/2.0 ),
					 std::ceil( (ub0[2]-lb0[2])/half_tgt_resl*sqrt(3.0)/2.0 )   );
			F3 lb = ( lb0 + ub0 - nc.template cast<float>() * half_tgt_resl/sqrt(3)*2.0 )/2.0;
			F3 ub = ( lb0 + ub0 + nc.template cast<float>() * half_tgt_resl/sqrt(3)*2.0 )/2.0;
			Director d( rot_resl_deg, lb, ub, nc, 1 );
			std::cout << "NEST info base resl: " << (ub-lb)/nc.template cast<float>() << " " << rot_resl_deg << std::endl;
			// {
				// cout << "NEST RAD " << rotamer_radius << endl;
				// cout << "NEST ROT " << rot_resl_deg << endl;
				// cout << "NEST LB0 " << lb0 << endl;
				// cout << "NEST UB0 " << ub0 << endl;
				// cout << "NEST NC  " << nc << endl;
				// cout << "NEST lb  " << lb << endl;
				// cout << "NEST ub  " << ub << endl;
				// cout << "NEST sp  " << (ub-lb)/nc.template cast<float>() << endl;
				// cout << "NUM   ORI CELLS " << KMGT( d.nest_.ori_map_.num_cells()) << " " << d.nest_.ori_map_.nside_ << endl;
				// cout << "NUM TRANS CELLS " << KMGT( d.nest_.trans_map_.num_cells()) << endl;
				// cout << "NUM   TOT CELLS " << KMGT( d.nest_.size(0)) << endl;
			// }



			std::vector< Scene > scene_per_thread( omp_max_threads_1() );
			for( auto & s : scene_per_thread ) s = scene_proto;

			Objective objective;

			std::vector< std::vector< SearchPoint > > samples( RESLS.size() );
				samples[0].resize( d.nest_.size(0) );
				for( uint64_t i = 0; i < d.nest_.size(0); ++i )	samples[0][i] = SearchPoint( i );

			for( int r = 0; r < RESLS.size()-1; ++r){
				if( 0 == samples[r].size() ) break;
				cout << "Hstage: " << r << " resl: " << F(4,2,RESLS[r]) << " nsamp: " << KMGT(samples[r].size()) << " ";
				int64_t const out_interval = samples[r].size()/50;
				std::exception_ptr exception = nullptr;
				#ifdef USE_OPENMP
				#pragma omp parallel for schedule(dynamic,8192)
				#endif
				for( int64_t i = 0; i < samples[r].size(); ++i ){
					if( exception ) continue;
					try {
						if( i%out_interval==0 ){
							cout << '*'; cout.flush();// (float)i/samples[r].size()*100.0 << "% "; cout.flush();
						}
						// uint64_t i = numeric::random::uniform()*d.nest_.size(0);
						uint64_t const isamp = samples[r][i].index;
						Scene & tscene( scene_per_thread[omp_get_thread_num()] );
						d.set_scene( isamp, r, tscene );
						// this is necssary, lots seem to have the same score
						samples[r][i].score = /*samples[r][i].rank =*/ objective( tscene, r ).template get<VoxelScore>();// - numeric::random::uniform()/1000.0;
					} catch( ... ) {
						#ifdef USE_OPENMP
						#pragma omp critical
						#endif
						exception = std::current_exception();
					}
				}
				if( exception ) std::rethrow_exception(exception);

				SearchPoint max_pt, min_pt;
				int64_t len = samples[r].size();
				if( samples[r].size() > beam_size/DIMPOW2 ){
					__gnu_parallel::nth_element( samples[r].begin(), samples[r].begin()+beam_size/DIMPOW2, samples[r].end() );
					len = beam_size/DIMPOW2;
					min_pt = *__gnu_parallel::min_element( samples[r].begin(), samples[r].begin()+len );
					max_pt = *(samples[r].begin()+beam_size/DIMPOW2);
				} else {
					min_pt = *__gnu_parallel::min_element( samples[r].begin(), samples[r].end() );
					max_pt = *__gnu_parallel::max_element( samples[r].begin(), samples[r].end() );
				}

				float const hsearch_score_cut = std::min( opts.abs_score_cut, abs_score_cut_by_res_thisres );
				omp_set_lock(&cout_lock);
				cout << " branching: " << F(9,6,min_pt.score) << " to " << F(9,6, std::min(hsearch_score_cut,max_pt.score)) << endl;
				omp_unset_lock(&cout_lock);

				// this hackyness is necessary.. don't want to explicidly build final samples vector... too big
				if( r+2 >= samples.size() ) break;

				// cout << "populating new sample array" << endl; // 26 sec for 2G... can speed up with threads?
				std::vector< std::vector< SearchPoint > > samples_thread(omp_max_threads_1());

				// #ifdef USE_OPENMP
				// #pragma omp parallel for schedule(dynamic,1024)
				// #endif
				// for( int64_t i = 0; i < len; ++i ){
				// 	if( samples[r][i].score > hsearch_score_cut ) continue;
				// 	uint64_t isamp0 = samples[r][i].index;
				// 	for( uint64_t j = 0; j < DIMPOW2; ++j ){
				// 		uint64_t isamp = isamp0 * DIMPOW2 + j;
				// 		samples_thread[omp_get_thread_num()].push_back( SearchPoint(isamp) );
				// 	}
				// }
				// for( int i = 0; i < samples_thread.size(); ++i ){
				// 	  std::copy ( samples_thread[i].begin(), samples_thread[i].end(), std::back_inserter( samples[r+1] ) );
				// }
				// samples[r+1].reserve( beam_size );
				for( int64_t i = 0; i < len; ++i ){
					if( samples[r][i].score > hsearch_score_cut ) continue;
					uint64_t isamp0 = samples[r][i].index;
					for( uint64_t j = 0; j < DIMPOW2; ++j ){
						uint64_t isamp = isamp0 * DIMPOW2 + j;
						samples[r+1].push_back( SearchPoint(isamp) );
					}
				}
				// cout << "done populating new sample array, clearing" << endl;
				samples[r].clear();

			}

			float const final_score_cut = std::min( opts.abs_score_cut, abs_score_cut_by_res_thisres );
			uint64_t num_final_samples = 0;{
				std::vector<uint64_t> tmp(omp_max_threads_1(),0);
				#ifdef USE_OPENMP
				#pragma omp parallel for schedule(dynamic,8192)
				#endif
				for( int64_t i = 0; i < samples[RESLS.size()-2].size(); ++i ){
					tmp[omp_get_thread_num()] += ( samples[RESLS.size()-2][i].score < final_score_cut );
				}
				for( int i = 0; i < tmp.size(); ++i ) num_final_samples += tmp[i];
			}


			// final
				float score_weight = 1.0;
				if( opts.downweight_hydrophobics ){
					score_weight = 0.8;
					if( resn == "TRP" ) score_weight = 0.5;
					if( resn == "PHE" ) score_weight = 0.6;
					if( resn == "TYR" ) score_weight = 0.6;
					if( resn == "MET" ) score_weight = 0.7;
				}

				typedef std::tuple<float,EigenXform,int> TestHit;
				std::vector<TestHit> test_hits;
				int r = RESLS.size()-1;
				cout << "Hstage: " << r << " resl: " << F(4,2,RESLS.back()) << " nsamp: " << KMGT(num_final_samples*DIMPOW2) << " ";
				int64_t const out_interval = samples[r-1].size()/50;
				float min_score = 9e9;
				std::vector<double> avg_scores( omp_max_threads_1(), 0.0 );
				std::vector<uint64_t> avg_scores_count( omp_max_threads_1(), 0 );
				std::exception_ptr exception = nullptr;

				int64_t block_size = 8192;
				for( int64_t iblock = 0; iblock < samples[r-1].size()/block_size; ++iblock )
				{
					int64_t block_begin = iblock * block_size;
					int64_t block_end = std::min<int64_t>( samples[r-1].size(), block_begin+block_size );

					#ifdef USE_OPENMP
					#pragma omp parallel for schedule(dynamic,1)
					#endif
					for( int64_t i = block_begin; i < block_end; ++i ){
						if(exception) continue;
						try{
							if( i%out_interval==0 ){
								cout << '*'; cout.flush();// (float)i/samples[r].size()*100.0 << "% "; cout.flush();
							}
							if( samples[r-1][i].score > final_score_cut ) continue;
							uint64_t isamp0 = samples[r-1][i].index;
							for( uint64_t j = 0; j < DIMPOW2; ++j ){
								uint64_t isamp = isamp0 * DIMPOW2 + j;
								Scene & tscene( scene_per_thread[omp_get_thread_num()] );
								d.set_scene( isamp, r, tscene );
								float score0 = objective( tscene, r ).template get<VoxelScore>();// - numeric::random::uniform()/1000.0;
								if(score0 < 0){
									avg_scores[ omp_get_thread_num() ] += score0;
									avg_scores_count[ omp_get_thread_num() ]++;
								}

								// BB atoms are repulsive-only, so this is ok
								if( score0 > final_score_cut ) continue;

								for( auto const & child : inv_rotamer_backbones )
								{
									int crot = child.rotid;
									Vector3f Nchild  = tscene.position(1) * child.Ncen;
									Vector3f CAchild = tscene.position(1) * child.CAcen;
									Vector3f Cchild  = tscene.position(1) * child.Ccen;
									::scheme::actor::BackboneActor<EigenXform> bbactor_child( Nchild, CAchild , Cchild );

									// if(runiftest < 0.0001) {
										// 	utility::io::ozstream out("test"+str(++count)+".pdb");


										// 	EigenXform x = tscene.position(1);
										// 	EigenXform y = EigenXform::Identity();
										// 	y.translation() = -rotamer_center;
										// 	EigenXform z = rot_index_p->to_structural_parent_frame_.at(crot);
										// 	rot_index_p->dump_pdb( out, crot, x*y*z );

										// 	dump_scene( d, tscene, *rot_index_p, isamp, RESLS.size()-1, out );

										// 	out << "MODEL" << std::endl;
										// 	::scheme::io::dump_pdb_atom_resname_atomname( out, "TST", "  N ",  Nchild );
										// 	::scheme::io::dump_pdb_atom_resname_atomname( out, "TST", " CA ", CAchild );
										// 	::scheme::io::dump_pdb_atom_resname_atomname( out, "TST", "  C ",  Cchild );
										// 	out << "ENDMDL" << std::endl;
										// 	out.close();
										// }

									float score = score0;
									// treat all bb atoms as 20 by convention, bb is repl-only
									score += std::max(0.0f,bounding_by_atype.back().at(20)->at( Nchild ));
									score += std::max(0.0f,bounding_by_atype.back().at(20)->at( CAchild ));
									score += std::max(0.0f,bounding_by_atype.back().at(20)->at( Cchild ));

									if( score < min_score ){
										#pragma omp critical
										if( score < min_score ){
											min_score = score;
										}
									}
									if( score > final_score_cut ) continue;

									accumulator->insert( bbactor_child.position_, score_weight*score, crot );

									if( opts.dump_fraction > 0 ){
										double const runif = uniform(rngs[omp_thread_num_1()-1]);
										if( runif < opts.dump_fraction ){
											omp_set_lock(&io_lock);
												test_hits.push_back( std::make_tuple(score,tscene.position(1),crot) );
											omp_unset_lock(&io_lock);
										}
									}

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

					// this is why blocks... can't condense while other threads running
					if( accumulator->need_to_condense() ){
						accumulator->checkpoint( cout );
					}

				} // end blocks

				float avg_score = 0.0;
				uint64_t tot_samp = 0;
				for( int i = 0; i < avg_scores.size(); ++i){
					avg_score += avg_scores[i] / avg_scores_count[i];
					tot_samp += avg_scores_count[i];
				}
				avg_score /= avg_scores.size();

				std::cout << endl << "SCOREINFO " << resn << " min: " << F(7,3,min_score)
				          << " cut: " << F(3,1,abs_score_cut_by_res_thisres) << " avg: " << F(7,3,avg_score) << " nsamp " << KMGT(tot_samp) << endl;


			accumulator->checkpoint( cout );
			accumulator->report( cout );

			if( test_hits.size() ){
				std::sort(test_hits.begin(), test_hits.end(),
					[](TestHit a, TestHit b) { return std::get<0>(a) < std::get<0>(b); });
				utility::io::ozstream out( params->output_prefix+"RifGen_Apo_test_hits_"+resn+boost::lexical_cast<std::string>(irot)+".pdb.gz" );
				for( auto h : test_hits ){
					float score;
					int irot;
					EigenXform x;
					std::tie(score,x,irot) = h;
					EigenXform y = EigenXform::Identity();
					y.translation() = -rotamer_center;
					EigenXform z = rot_index_p->to_structural_parent_frame_.at(irot);
					rot_index_p->dump_pdb( out, irot, x*y*z );
				}
				out.close();
			}
			test_hits.clear();

			samples.back().clear();

			// rif_apo_vis_out.close();

		}


		omp_destroy_lock( & cout_lock ) ;
		omp_destroy_lock( & io_lock );
		omp_destroy_lock( & accum_lock );


	}


}
}
}

