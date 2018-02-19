// #include <core/scoring/motif/motif_hash_stuff.hh>
// #include <core/import_pose/import_pose.hh>
#include <basic/options/option_macros.hh>

// #include <core/chemical/ChemicalManager.hh>
// #include <core/pose/PDBInfo.hh>

// #include <core/scoring/ScoreFunctionFactory.hh>
// #include <core/pose/Pose.hh>
// #include <core/io/pdb/pose_io.hh>
// #include <core/import_pose/import_pose.hh>
// #include <core/id/AtomID.hh>
#include <devel/init.hh>

#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/file/file_sys_util.hh>
#include <ObjexxFCL/format.hh>


#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <map>

#include <riflib/RifFactory.hh>

#include <riflib/util.hh>



using std::cout;
using std::endl;
using devel::scheme::KMGT;


		// Float lever_radius,
		// Float lever_dis,
		// XformHashNeighbors<Hasher> & nbcache

OPT_1GRP_KEY( StringVector, scheme, base_xmap_files             )
	OPT_1GRP_KEY( Real              , scheme, base_xmap_cart_resl )
	OPT_1GRP_KEY( Real              , scheme, base_xmap_ang_resl )
	OPT_1GRP_KEY( RealVector        , scheme, hash_cart_resls        )
	OPT_1GRP_KEY( RealVector        , scheme, hash_cart_bounds       )
	OPT_1GRP_KEY( RealVector        , scheme, hash_ang_resls         )
	OPT_1GRP_KEY( RealVector        , scheme, lever_radii      )
	OPT_1GRP_KEY( RealVector        , scheme, lever_bounds     )
	OPT_1GRP_KEY( String            , scheme, nbcache_path )
	OPT_1GRP_KEY( Boolean           , scheme, verbose          )
	OPT_1GRP_KEY( Integer           , scheme, nbcache_rotation_sample_factor          )
	OPT_1GRP_KEY( Boolean           , scheme, make_hacky_grids )

	void REGISTER_OPTIONS() {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		NEW_OPT( scheme::base_xmap_files, "xmap files to make 'coarse' grids for"     , utility::vector1<std::string>() );
		NEW_OPT( scheme::base_xmap_cart_resl, "cart resl for input xmap", 0 );
		NEW_OPT( scheme::base_xmap_ang_resl, "ang resl for input xmap", 0 );
		NEW_OPT( scheme::hash_cart_resls, "cartesian resolution(s) of hash table(s)"      , utility::vector1<double>() );
		NEW_OPT( scheme::hash_cart_bounds, "bound on cartesian coordinates"               , utility::vector1<double>() );
		NEW_OPT( scheme::hash_ang_resls,  "ang reslolution(s) of hash table(s) in degrees", utility::vector1<double>() );
		NEW_OPT( scheme::lever_radii      , ""                                      , utility::vector1<double>() );
		NEW_OPT( scheme::lever_bounds     , ""                                      , utility::vector1<double>() );
		NEW_OPT( scheme::nbcache_path , ""                                      , ".");
		NEW_OPT( scheme::verbose, "", false );
		NEW_OPT( scheme::nbcache_rotation_sample_factor, "", 1000 );
		NEW_OPT( scheme::make_hacky_grids, "", true );
	}


void
make_bounding_grids(){

	using namespace basic::options;
	using ObjexxFCL::format::F;
	namespace sopt = basic::options::OptionKeys::scheme;
	using namespace devel::scheme;


	runtime_assert( option[sopt::lever_radii     ]().size() == option[sopt::lever_bounds     ]().size() );
	runtime_assert( option[sopt::lever_radii     ]().size() == option[sopt::hash_ang_resls   ]().size() );
	runtime_assert( option[sopt::lever_radii     ]().size() == option[sopt::hash_cart_resls  ]().size() );
	runtime_assert( option[sopt::lever_radii     ]().size() == option[sopt::hash_cart_bounds ]().size() );

	runtime_assert( option[sopt::base_xmap_files]().size() == 1 );

	std::string base_xmap_file = option[sopt::base_xmap_files]().front();
	runtime_assert( utility::file::file_exists( base_xmap_file ) );
	runtime_assert( base_xmap_file.substr(base_xmap_file.size()-8) == ".xmap.gz" ||
	                base_xmap_file.substr(base_xmap_file.size()-7) == ".rif.gz"  ||
	                base_xmap_file.substr(base_xmap_file.size()-5) == ".xmap"  ||
	                base_xmap_file.substr(base_xmap_file.size()-4) == ".rif"
	                );
	cout << "reading " << base_xmap_file << endl;

	std::string rif_type = get_rif_type_from_file( base_xmap_file );
	std::cout << "======================================================================" << std::endl;
	std::cout << "read RIF type: " << rif_type << std::endl;
	std::cout << "======================================================================" << std::endl;
	devel::scheme::RifFactoryConfig rif_factory_config;
	rif_factory_config.rif_type = rif_type;
	shared_ptr<RifFactory> rif_factory = ::devel::scheme::create_rif_factory( rif_factory_config );

	std::string ref_description;
	RifConstPtr ref_rif = rif_factory->create_rif_from_file( base_xmap_file, ref_description );
	std::cout << "read description: " << std::endl << ref_description << std::endl;
	std::cout << "======================================================================" << std::endl;

	std::cout << "RIF size: " << KMGT( ref_rif->mem_use() ) << " load: " << ref_rif->load_factor() << std::endl;

	int nbase = ref_rif->size();
	std::cout <<"read full xmap, size: " << KMGT(nbase) << std::endl;


	for( int ibound = 1; ibound <= option[sopt::lever_bounds]().size(); ++ibound ){

		double const lever_radius      = option[sopt::lever_radii      ]().at( ibound );
		double const lever_bound       = option[sopt::lever_bounds     ]().at( ibound );
		double const hash_ang_resl     = option[sopt::hash_ang_resls   ]().at( ibound );
		double const hash_cart_resl    = option[sopt::hash_cart_resls  ]().at( ibound );
		double const hash_cart_bound   = option[sopt::hash_cart_bounds ]().at( ibound );

		double const cart_bound = lever_bound;
		double const  ang_bound_rad = lever_bound / lever_radius;
		double const  ang_bound = ang_bound_rad * 180.0 / M_PI;
		cout << "========================================= make bounding grids for resl. " << lever_bound << " ==============================================" << endl;
		cout << "cart_bound: " << cart_bound << ", ang_bound: " << ang_bound << endl;
		cout << "cart_hash_resl: " << hash_cart_resl << ", hash_ang_resl: " << hash_ang_resl << std::endl;

		RifPtr new_rif = rif_factory->create_rif_from_rif( ref_rif, hash_cart_resl, hash_ang_resl, hash_cart_bound );
		std::cout << "new map size " << KMGT(new_rif->size()) << " ratio: " <<  (float)new_rif->size() / (float)ref_rif->size() << std::endl;

		{
			cout << "before resize, bounding grid size: " << KMGT( new_rif->mem_use() )
			     << " load: " << new_rif->load_factor()
				 << ", sizeof(value_type) " << new_rif->sizeof_value_type() << endl;
			 cout << "NOT RESIZING bounding grid" << endl;
			// bounding_xmapfull.map_.max_load_factor(0.8);
			// bounding_xmapfull.map_.min_load_factor(0.4);
			// bounding_xmapfull.map_.resize(0);
			// cout << "after  resize, bounding grid size: " << KMGT( bounding_xmapfull.mem_use() )
			     // << " load: " << bounding_xmapfull.map_.size()*1.f/bounding_xmapfull.map_.bucket_count() << endl;
			std::string digits = boost::lexical_cast<std::string>(lever_bound);
			if( digits.size() == 1 ) digits = "0" + digits;
			std::string tag = "_BOUNDING_";
			if( option[sopt::make_hacky_grids]() ) tag += "HACK_";
			else                                   tag += "RANDHACK_";
			std::ostringstream oss_description;
			oss_description << "==== bounding xmap ====" << endl;
			if( option[sopt::make_hacky_grids]() ){
				oss_description << "!!!!!!!!!!!!!!!!!!!!!!!!!!!! USING_HACKY_GRIDS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
			}
			oss_description << "hash_cart_bound: " << hash_cart_bound << endl;
			oss_description << "lever_radius:  " << lever_radius << endl;
			oss_description << "lever_bound:   " << lever_bound << endl;
			oss_description << "==== source oss_description ====\n" << ref_description;
			std::string description = oss_description.str();

			utility::io::ozstream out( base_xmap_file+tag+"RIF_"+digits + ".xmap.gz" );
			new_rif->save( out, description );
			out.close();
		}



	}

}


int main(int argc, char *argv[])
{

	using namespace ::devel::scheme;

	REGISTER_OPTIONS();
	devel::init(argc,argv);

	make_bounding_grids();

	return 0;
}






/*

template<class T>
float fullval_to_numeric( T const & t ){
	return t;
}

template< int N >
float fullval_to_numeric( ::scheme::objective::storage::RotamerScores<N> const & rs ){
	return rs.score(0);
}

////////////// depricated... left for possible future reference
template<
	class EigenXform,
	class XMapValFull,
	template<class Xform> class HasherFull ,
	template<class Xform> class HasherBounding
>
void
make_bounding_grids_failed(){

	int const ArrayBits1 = 0;
	int const ArrayBits2 = 0;
	typedef float XMapValBounding;

	typedef ::scheme::objective::hash::XformMap< EigenXform, XMapValFull, HasherFull > XMapFull;

	using namespace basic::options;
	using namespace core::scoring::motif;
	namespace sopt = basic::options::OptionKeys::scheme;
	using ObjexxFCL::format::F;
	using devel::scheme::omp_max_threads;
	using devel::scheme::omp_thread_num;

	boost::random::mt19937 rng((unsigned int)time(0) + 736495684);
	int const num_threads = omp_max_threads();

	typedef scheme::objective::hash::XformMap< EigenXform, XMapValBounding, HasherBounding > XMapBounding;
	typedef scheme::objective::hash::XformHashNeighbors<typename XMapBounding::Hasher> XHNB;

	runtime_assert( option[sopt::lever_radii     ]().size() == option[sopt::lever_bounds     ]().size() );
	runtime_assert( option[sopt::lever_radii     ]().size() == option[sopt::hash_ang_resls   ]().size() );
	runtime_assert( option[sopt::lever_radii     ]().size() == option[sopt::hash_cart_resls  ]().size() );
	runtime_assert( option[sopt::lever_radii     ]().size() == option[sopt::hash_cart_bounds ]().size() );

	XMapFull full_xmap; // read resls from file unless artificial
	std::string full_xmap_description;
	std::string base_xmap_file = "artificial_test";
	if( option[sopt::base_xmap_files]().size() > 0 ){
		base_xmap_file = option[sopt::base_xmap_files]().front();
		cout << "reading " << base_xmap_file << endl;
		runtime_assert( option[sopt::base_xmap_files]().size() == 1 );
		runtime_assert( utility::file::file_exists( base_xmap_file ) )
		runtime_assert( base_xmap_file.substr(base_xmap_file.size()-8) == ".xmap.gz" ||
		                base_xmap_file.substr(base_xmap_file.size()-7) == ".rif.gz"  ||
		                base_xmap_file.substr(base_xmap_file.size()-5) == ".xmap"  ||
		                base_xmap_file.substr(base_xmap_file.size()-4) == ".rif"
		                );
		utility::io::izstream in( base_xmap_file );
		runtime_assert( full_xmap.load( in, full_xmap_description ) );
		in.close();
		std::cout << "RIF size: " << KMGT( full_xmap.mem_use() ) << " load: " << full_xmap.map_.size()*1.f/full_xmap.map_.bucket_count() << std::endl;
	} else if( option[sopt::test_structure].user() ) {
		cout << "make artificial test RIF from: " << option[sopt::test_structure]() << endl;
		full_xmap.init( option[sopt::base_xmap_cart_resl](), option[sopt::base_xmap_ang_resl]() );
		core::pose::Pose pose;
		core::import_pose::pose_from_pdb( pose, option[sopt::test_structure]() );
		for( int itest = 0; itest < option[sopt::num_test_points](); ++itest ){
			{
				numeric::xyzVector<double> cen   = protocols::sic_dock::center_of_geom( pose );
				numeric::xyzVector<double> trans = numeric::random::random_vector_unit_cube() * 20.0 - 10.0;
				numeric::xyzMatrix<double> rot   = numeric::random::random_rotation();
				protocols::sic_dock::trans_pose( pose, -cen );
				protocols::sic_dock::rot_pose( pose, rot );
				protocols::sic_dock::trans_pose( pose, trans );
			}
			std::vector<int> picked_res;
			for( int ir = 1; ir <= pose.size(); ++ir ) picked_res.push_back(ir);
			numeric::random::random_permutation( picked_res, numeric::random::rg() );
			std::vector<EigenXform> target_frames;
			float score = 2.0 - (double)itest / ((double)option[sopt::num_test_points]()-1);
			std::cout << itest << " " << score << std::endl;
			for( int i = 0; i < option[sopt::num_target_frames](); ++i ){
				int ir = picked_res[i];
				scheme::actor::BackboneActor<EigenXform> bbactor( pose.residue(ir).xyz("N"), pose.residue(ir).xyz("CA"), pose.residue(ir).xyz("C") );
				::scheme::objective::storage::RotamerScores< 12 > rs;
				rs.add_rotamer( 1, score );
				full_xmap.insert( bbactor.position(), rs );
			}
			pose.dump_pdb( "test"+ObjexxFCL::format::I(4,itest+1)+".pdb" );
		}
	} else {
		utility_exit_with_message( "must specify base_xmap_files OR test_structure" );
	}
	// int nbase = full_xmap.count_not( 0.0 );
	// int nbase = full_xmap.count_not( XMapValFull() );
	int nbase = full_xmap.size();
	std::cout <<"read full xmap, size: " << KMGT(nbase) << std::endl;


	for( int ibound = 1; ibound <= option[sopt::lever_bounds]().size(); ++ibound ){

		double const lever_radius      = option[sopt::lever_radii      ]().at( ibound );
		double const lever_bound       = option[sopt::lever_bounds     ]().at( ibound );
		double const hash_ang_resl     = option[sopt::hash_ang_resls   ]().at( ibound );
		double const hash_cart_resl    = option[sopt::hash_cart_resls  ]().at( ibound );
		double const hash_cart_bound   = option[sopt::hash_cart_bounds ]().at( ibound );
		std::string nbcache_file = option[sopt::nbcache_path]() + "/__scheme_make_bounding_grids_"
			+"cr"+F(5,2,hash_cart_resl)+"_"
			+"ar"+F(5,2,hash_ang_resl)+"_"
			+"lr"+F(5,2,lever_radius)+"_"
			+"lb"+F(5,2,lever_bound)+".nbcache.gz";


		double const cart_bound = lever_bound;
		double const  ang_bound_rad = lever_bound / lever_radius;
		double const  ang_bound = ang_bound_rad * 180.0 / M_PI;
		cout << "========================================= make bounding grids for resl. " << lever_bound << " ==============================================" << endl;
		cout << "cart_bound: " << cart_bound << ", ang_bound: " << ang_bound << endl;
		cout << "cart_hash_resl: " << hash_cart_resl << ", hash_ang_resl: " << hash_ang_resl << std::endl;

		std::cout << "insert all values into coarse_xmap "; std::cout.flush();
		// XMapBounding coarse_xmap( hash_cart_resl, hash_ang_resl, hash_cart_bound );

		XMapFull * coarse_xmapfull_p;

		if( option[sopt::make_hacky_grids]() || lever_bound > 3.0 ){

			coarse_xmapfull_p = new XMapFull( hash_cart_resl, hash_ang_resl, hash_cart_bound );
			int progress0 = 0;
			BOOST_FOREACH( typename XMapFull::Map::value_type const & v, full_xmap.map_ ){
				if( ++progress0 % std::max((size_t)1,(full_xmap.size()/100)) == 0 ){
					std::cout << '*'; std::cout.flush();
				}
				EigenXform x = full_xmap.hasher_.get_center( v.first );
				XMapValBounding val = fullval_to_numeric( v.second );
				if( val != 0.0 ){
					// coarse_xmap.insert( x, val );
					uint64_t k = coarse_xmapfull_p->hasher_.get_key(x);
					typename XMapFull::Map::iterator iter = coarse_xmapfull_p->map_.find(k);
					if( iter == coarse_xmapfull_p->map_.end() ){
						coarse_xmapfull_p->map_.insert( std::make_pair(k,v.second) );
					} else {
						iter->second.merge( v.second );
					}
					// bounding_xmap.insert_sphere( x, lever_radius, lever_bound, val, nbcache );
				}
			}
			std::cout << endl;
			std::cout << "coarse_map size " << KMGT(coarse_xmapfull_p->size()) << " ratio: " <<  (float)coarse_xmapfull_p->size() / (float)full_xmap.size() << std::endl;

		} else {
			coarse_xmapfull_p = & full_xmap;
		}

		// XMapBounding bounding_xmap( hash_cart_resl, hash_ang_resl, hash_cart_bound );
		XMapFull bounding_xmapfull( hash_cart_resl, hash_ang_resl, hash_cart_bound );

		if( option[sopt::make_hacky_grids]() ){

			// bounding_xmap = coarse_xmap;
			bounding_xmapfull = *coarse_xmapfull_p;
			std::cout << "make_hacky_grids, using coarse_xmap as bounding map" << std::endl;

		} else if( true ) {
			std::cout << "using better hack to make bounding grids" << std::endl;


			float cart_rand_spread = cart_bound - hash_cart_resl / 2.0; // no justification for this...
			float ang_rand_spread  =  ang_bound -  hash_ang_resl / 2.0;
			std::cout << "cart_rand_spread: " << cart_rand_spread << ", " << "ang_rand_spread: " << ang_rand_spread << endl;

			float cratio = cart_bound / hash_cart_resl;
			float aratio = ang_bound / hash_ang_resl;

			int Nrandxform = 100.0  *  cratio*cratio*cratio  *  aratio*aratio*aratio;
			std::cout << "N insertion samples: " << KMGT(Nrandxform) << std::endl;
			if( lever_radius < 3.0 ) Nrandxform /= 3.0;

			std::vector<EigenXform> rand_xforms(Nrandxform);
			rand_xforms[0] = EigenXform::Identity();
			for( int i = 1; i < Nrandxform; ++i ){
				::scheme::numeric::rand_xform_sphere( rng, rand_xforms[i], cart_rand_spread, float(ang_rand_spread * M_PI/180.0) );
			}

			std::vector<XMapFull> bounding_xmaps( omp_max_threads(), XMapFull( hash_cart_resl, hash_ang_resl, hash_cart_bound ) );

			std::vector< std::pair< typename XMapFull::Map::key_type, typename XMapFull::Map::data_type > > coarse_xmap_values;
			// std::vector< typename XMapFull::Map::value_type > coarse_xmap_values;
			coarse_xmap_values.reserve(coarse_xmapfull_p->size());
			BOOST_FOREACH( typename XMapFull::Map::value_type const & v, coarse_xmapfull_p->map_ ){
				coarse_xmap_values.push_back( v );
			}

			std::cout << "make grid with " << KMGT(Nrandxform) << " rand nbr samp. ";
			float avg_nbcount = 0;
			#ifdef USE_OPENMP
			#pragma omp parallel for schedule(dynamic,1)
			#endif
			for( int i = 0; i < coarse_xmap_values.size(); ++i ){
				int ithread = omp_thread_num();
				XMapFull & bounding_xmap( bounding_xmaps[ithread] );
				typename XMapFull::Map::data_type const val = coarse_xmap_values[i].second;
				EigenXform x0 = coarse_xmapfull_p->hasher_.get_center( coarse_xmap_values[i].first );
				// std::set<uint64_t> seenit;
				google::dense_hash_set<uint64_t> seenit;
				seenit.set_empty_key(0);
				for( int irand = 0; irand < Nrandxform; ++irand ){
					EigenXform const x = rand_xforms[irand] * x0;
					uint64_t const k = bounding_xmap.hasher_.get_key(x);
					if( seenit.find(k) == seenit.end() ){
						seenit.insert(k);
						typename XMapFull::Map::iterator iter = bounding_xmap.map_.find(k);
						if( iter == bounding_xmap.map_.end() ){
							bounding_xmap.map_.insert( std::make_pair( k, val ) );
						} else {
							iter->second.merge( val );
						}
					}
				}
				#pragma omp critical
				avg_nbcount += seenit.size();

				if( i % std::max((size_t)1,(coarse_xmapfull_p->size()/100)) == 0 ){
					std::cout << "*"; std::cout.flush();
				}
			}
			std::cout << std::endl;
			avg_nbcount /= coarse_xmap_values.size();
			std::cout << "avg_nbcount " << avg_nbcount << ", nsamp / nbinfound = " << Nrandxform / avg_nbcount << std::endl;

			if( lever_radius > 3.0 ) delete coarse_xmapfull_p;

// utility_exit_with_message("testing rand neighbor gen");

			int n = 1; while( n < bounding_xmaps.size() ) n *= 2;
			while( n > 1 ){
				std::cout << "condense per-thread maps, to merge: " << n << std::endl;
				#ifdef USE_OPENMP
				#pragma omp parallel for schedule(dynamic,1)
				#endif
				for( int i = 0; i < n; ++i ){
					if( i < n/2 && i + n/2 < bounding_xmaps.size() ){
						XMapFull       & merge_to   = bounding_xmaps[i];
						XMapFull const & merge_from = bounding_xmaps[n/2+i];
						BOOST_FOREACH( typename XMapFull::Map::value_type const & v, merge_from.map_ ){
							typename XMapFull::Map::iterator iter = merge_to.map_.find( v.first );
							if( iter == merge_to.map_.end() ){
								merge_to.map_.insert( v );
							} else {
								iter->second.merge( v.second );
							}
						}
						// std::cout << "merged " << i+n/2 << " into " << i << std::endl;
					}
				}
				n /= 2;
			}
			bounding_xmapfull = bounding_xmaps[0];

		} else { // the way that doesn't work yet, using nb cells from library funcs

	// 		std::cout << "now do insert_sphere for all " << KMGT( coarse_xmap.size() ) << " values in coarse xmap with " << num_threads << " threads" << std::endl;

	// 		std::vector<XMapBounding> bounding_xmaps( omp_max_threads(), XMapBounding( hash_cart_resl, hash_ang_resl, hash_cart_bound ) );
	// 		XHNB nbcache_tmp( cart_bound, ang_bound, bounding_xmaps.front().hasher_, option[sopt::nbcache_rotation_sample_factor]() );
	// 		if( utility::file::file_exists( nbcache_file ) ){
	// 			cout << "reading nbcache file: " << nbcache_file << endl;
	// 			utility::io::izstream in( nbcache_file );
	// 			nbcache_tmp.load(in);
	// 			in.close();
	// 		}
	// 		std::vector<XHNB> nbcaches( num_threads, nbcache_tmp );

	// 		std::vector< std::pair< typename XMapBounding::Map::key_type, typename XMapBounding::Map::data_type > > coarse_xmap_values;
	// 		coarse_xmap_values.reserve(coarse_xmap.size());
	// 		BOOST_FOREACH( typename XMapBounding::Map::value_type const & v, coarse_xmap.map_ ){
	// 			coarse_xmap_values.push_back(v);
	// 		}
	// 		std::vector<double> avg_nbcounts(omp_max_threads(),0.0);

	// 		#ifdef USE_OPENMP
	// 		#pragma omp parallel for schedule(dynamic,1)
	// 		#endif
	// 		for( int i = 0; i < coarse_xmap_values.size(); ++i ){
	// 			int ithread = omp_thread_num();
	// 			XMapBounding & bounding_xmap( bounding_xmaps[ithread] );
	// 			XHNB & nbcache( nbcaches[ithread] );
	// 			double avg_nbcount = 0;
	// 			XMapValBounding val = coarse_xmap_values[i].second;
	// 			EigenXform x = coarse_xmap.hasher_.get_center( coarse_xmap_values[i].first );
	// 			avg_nbcounts[ithread] += bounding_xmap.insert_sphere( x, lever_radius, lever_bound, val, nbcache );
	// 			if( i % std::max((size_t)1,(coarse_xmap.size()/100)) == 0 ){
	// 				std::cout << "progress " << (float)i / (float)coarse_xmap_values.size() *100.0 << "% " << std::endl;
	// 			}
	// 		}
	// 		double avg_nbcount = 0;
	// 		BOOST_FOREACH( double d, avg_nbcounts ) avg_nbcount += d;
	// 		avg_nbcount /= double( coarse_xmap_values.size() );


	// // test = [(str(i)+" ") for i in range(12)]
	// // print
	// // n = 1
	// // while n < len(test): n *= 2
	// // while True:
	// // 	print "=======",n,"======="
	// // 	for i in range(n):
	// // 		if i < n/2 and i+n/2 < len(test):
	// // 			print i,n/2+i
	// // 			test[i] += test[n/2+i]
	// // 	n /= 2

	// // 	if n == 1: break

	// // #for i in test: print i
	// 		// hacky parallel merge
	// 		int n = 1; while( n < bounding_xmaps.size() ) n *= 2;
	// 		while( n > 1 ){
	// 			std::cout << "condense per-thread maps, to merge: " << n << std::endl;
	// 			#ifdef USE_OPENMP
	// 			#pragma omp parallel for schedule(dynamic,1)
	// 			#endif
	// 			for( int i = 0; i < n; ++i ){
	// 				if( i < n/2 && i + n/2 < bounding_xmaps.size() ){
	// 					XMapBounding       & merge_to   = bounding_xmaps[i];
	// 					XMapBounding const & merge_from = bounding_xmaps[n/2+i];
	// 					BOOST_FOREACH( typename XMapBounding::Map::value_type const & v, merge_from.map_ ){
	// 						typename XMapBounding::Map::iterator iter = merge_to.map_.find( v.first );
	// 						if( iter == merge_to.map_.end() ){
	// 							merge_to.map_.insert( v );
	// 						} else {
	// 							iter->second = std::min( v.second, iter->second );
	// 						}
	// 					}
	// 					std::cout << "merged " << i+n/2 << " into " << i << std::endl;
	// 				}
	// 			}
	// 			n /= 2;
	// 		}
	// 		bounding_xmap = bounding_xmaps[0];

	// 		// std::cout << "condense per-thread maps" << std::endl;
	// 		// XMapBounding bounding_xmap( hash_cart_resl, hash_ang_resl, hash_cart_bound );
	// 		// BOOST_FOREACH( XMapBounding const & xm, bounding_xmaps ){
	// 		// 	BOOST_FOREACH( typename XMapBounding::Map::value_type const & v, xm.map_ ){
	// 		// 		typename XMapBounding::Map::iterator iter = bounding_xmap.map_.find( v.first );
	// 		// 		if( iter == bounding_xmap.map_.end() ){
	// 		// 			bounding_xmap.map_.insert( v );
	// 		// 		} else {
	// 		// 			iter->second = std::min( v.second, iter->second );
	// 		// 		}
	// 		// 	}
	// 		// }

	// 		std::cout << endl << "average nbcount: " << avg_nbcount << " hash mem use: " << KMGT(bounding_xmap.mem_use()) << std::endl;

	// 		std::cout << "save nbcache_file " << nbcache_file << std::endl;
	// 		{
	// 			BOOST_FOREACH( XHNB const & nbcache, nbcaches ){
	// 				nbcache_tmp.merge( nbcache );
	// 			}
	// 			utility::io::ozstream out( nbcache_file );
	// 			nbcache_tmp.save( out );
	// 			out.close();
	// 		}

		} // end if hacky


		{
			cout << "before resize, bounding grid size: " << KMGT( bounding_xmapfull.mem_use() )
			     << " load: " << bounding_xmapfull.map_.size()*1.f/bounding_xmapfull.map_.bucket_count()
				 << ", sizeof(value_type) " << sizeof(typename XMapFull::Map::value_type) << endl;
			 cout << "NOT RESIZING bounding grid" << endl;
			// bounding_xmapfull.map_.max_load_factor(0.8);
			// bounding_xmapfull.map_.min_load_factor(0.4);
			// bounding_xmapfull.map_.resize(0);
			// cout << "after  resize, bounding grid size: " << KMGT( bounding_xmapfull.mem_use() )
			     // << " load: " << bounding_xmapfull.map_.size()*1.f/bounding_xmapfull.map_.bucket_count() << endl;
			std::string digits = boost::lexical_cast<std::string>(lever_bound);
			if( digits.size() == 1 ) digits = "0" + digits;
			std::string tag = "_BOUNDING_";
			if( option[sopt::make_hacky_grids]() ) tag += "HACK_";
			else                                   tag += "RANDHACK_";
			std::ostringstream description;
			description << "==== bounding xmap ====" << endl;
			if( option[sopt::make_hacky_grids]() ){
				description << "!!!!!!!!!!!!!!!!!!!!!!!!!!!! USING_HACKY_GRIDS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
			}
			description << "hash_cart_bound: " << hash_cart_bound << endl;
			description << "lever_radius:  " << lever_radius << endl;
			description << "lever_bound:   " << lever_bound << endl;
			description << "==== source description ====\n" << full_xmap_description;

			// utility::io::ozstream out( utility::file_basename(base_xmap_file)+tag+digits + ".xmap.gz" );
			// bounding_xmapfull.save( out, description.str() );
			// out.close();

			utility::io::ozstream out2( base_xmap_file+tag+"RIF_"+digits + ".xmap.gz" );
			bounding_xmapfull.save( out2, description.str() );
			out2.close();
		}


		// int nbounding = bounding_xmapfull.count_not(0.0);
		// cout << "Bounding size cost " << (float)nbounding/nbase << endl;
		// float fracfill = (float)nbounding/bounding_xmapfull.map_.size() / (float)(1<<ArrayBits2);
		// float memper = (float)((1<<ArrayBits2)+8.0) / (1<<ArrayBits2);

		// cout << "array fill frac " << fracfill << " zorder array[" << (1<<ArrayBits2) << "] mem_saving = " << 9.0 / memper * fracfill << "-fold " << endl;
		// some very quick test data:
		//            1   2    4    8   16   32   64   128  256  512  1024
		// R2.0/1.0: 1.0 1.53 2.05 2.54 2.79 2.34 1.97 1.53 1.17 0.90 0.66
		// R3.0/1.0:               2.23 2.39 1.99      1.28           0.52
		// R1.0/0.5                3.20 3.60 3.05
		// R1.5/0.5                3.40 4.00 3.45
		// R2.5/0.5                3.11 3.68


	}

}


*/









