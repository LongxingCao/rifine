// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols

#include <riflib/rosetta_field.hh>

	#include <scheme/objective/hash/XformMap.hh>
	#include <scheme/actor/BackboneActor.hh>
	#include <scheme/chemical/RotamerIndex.hh>
	#include <core/conformation/Residue.hh>
	#include <scheme/rosetta/score/RosettaField.hh>
	#include <riflib/util.hh>
	#include <riflib/EtableParams_init.hh>

	#include <core/chemical/AtomType.hh>

	#include <utility/io/izstream.hh>
	#include <utility/io/ozstream.hh>
	#include <utility/file/file_sys_util.hh>

	#include <ObjexxFCL/format.hh>

	#include <core/types.hh>

	#include <core/pose/Pose.hh>

	#include <boost/foreach.hpp>

	#include <Eigen/Dense>

	#include <exception>
	#include <stdexcept>


namespace devel {
namespace scheme {

using std::endl;
using std::cout;

int N_ATYPE = 21;








void
parse_atomids(
	core::pose::Pose const & refpose,
	utility::vector1<int> const & atomreslist,
	std::vector<core::id::AtomID> & atomids_out,
	std::string printme
){
	runtime_assert_msg( atomreslist.size() % 2 == 0, "atomno/res list must be even length (atomno/rsd pairs)" );
	int ia;
	for( int iro = 1; iro <= atomreslist.size(); ++iro ){
		if( iro%2==1 ){
			 ia = atomreslist[iro];
		} else {
			int ir = atomreslist[iro];
			runtime_assert_msg( 0 < ir && ir <= refpose.size(), "residue index out of bounds: " + str(ir) + " nres: " + str( refpose.size()) );
			runtime_assert_msg( 0 < ia && ia <= refpose.residue(ir).natoms(), "atom index out of bounds: " + str(ia) + " natom: " + str( refpose.residue(ir).natoms()) + " res: " + str(ir) );
			atomids_out.push_back( core::id::AtomID(ia,ir) );
			if( printme.size() ) std::cout << "parse_atomids selected atom " << printme << ": " << ir << " " << refpose.residue(ir).name3() << " " << ia << " " << refpose.residue(ir).atom_name(ia) << std::endl;
		}
	}
}


void
get_scheme_atoms(
	core::pose::Pose const & target,
	utility::vector1<core::Size> const & target_res,
	std::vector< ::scheme::actor::Atom< Eigen::Vector3f > > & out,
	bool bbonly
){
	// if(target_res.size()==0) for(core::Size i=1; i<=target.size(); ++i) target_res.push_back(i);

		std::vector<int> atypemap = get_rif_atype_map();

		for(int ires = 1; ires <= target_res.size(); ++ires){
			int ir = target_res[ires];
			core::conformation::Residue const & r( target.residue(ir) );
			int Natoms = r.nheavyatoms();
			if( bbonly ){
				Natoms = r.last_backbone_atom();
				// std::cout << Natoms << " " << r.atom_name(Natoms) << std::endl;
				if( r.atom_name(Natoms+1) == " CB " ) ++Natoms;
				if( r.name3()=="PRO" ) Natoms = r.nheavyatoms();
				if( r.name3()=="CYD" ) Natoms = r.nheavyatoms();
			}
			for(int ia = 1; ia <= Natoms; ++ia){
				if( r.is_virtual(ia) ){
					std::cout << "skipping virtual atom: " << r.name() << " " << r.atom_name(ia) << std::endl;
					continue;
				}
				int at = atypemap[ r.atom_type_index(ia) ];
				if( at > N_ATYPE ){
					// utility_exit_with_message("heavy atom type > N_ATYPE: "+str(at)+" "+r.name()+" "+r.atom_name(ia) );
					std::cout << "WARNING: heavy atom type "<<r.atom_type_index(ia)<<" > N_ATYPE: "+str(at)+" "+r.name()+" "+r.atom_name(ia) 
					          << " will treat as carbon for sterics!" << std::endl;
					at = 5;
				}
				::scheme::actor::Atom< Eigen::Vector3f > a( r.xyz(ia), at );
				out.push_back(a);
			}
		}
}


void
get_scheme_atoms_cbonly(
	core::pose::Pose const & target,
	utility::vector1<core::Size> const & target_res,
	std::vector< ::scheme::actor::Atom< Eigen::Vector3f > > & out
){
	for(int ires = 1; ires <= target_res.size(); ++ires){
		int ir = target_res[ires];
		core::conformation::Residue const & r( target.residue(ir) );
		if( !r.is_protein() ) continue;
		::scheme::actor::Atom< Eigen::Vector3f > a;
		if( r.has("CB") ) a = ::scheme::actor::Atom< Eigen::Vector3f >( r.xyz("CB"), 5 ); // CB ONLY
		else              a = ::scheme::actor::Atom< Eigen::Vector3f >( r.xyz("CA"), 2 ); // CA ONLY
		out.push_back(a);
	}
}


void
get_scheme_atoms(
	core::pose::Pose const & target,
	std::vector< ::scheme::actor::Atom< Eigen::Vector3f > > & out,
	bool bbonly
){
	utility::vector1<core::Size> target_res;
	if(target_res.size()==0){
		for(core::Size i = 1; i <= target.size(); ++i){
			target_res.push_back(i);
		}
	}
	get_scheme_atoms( target, target_res, out, bbonly );
}

std::string
get_rosetta_fields(
	std::string const & target_fname,
	core::pose::Pose const & target,
	utility::vector1<core::Size> const & target_res,
	RosettaFieldOptions const & opts,
	std::vector<
		// shared_ptr<
			::scheme::objective::voxel::VoxelArray<3,float> *
		// >
	> & field_by_atype,
	bool verbose
){

	std::string target_tag = utility::file_basename( target_fname );
	std::string hashstr = ::devel::scheme::get_res_list_hash( target_res );

	// std::cout << target_res.size() << " " << target_res_hash << std::endl;
	// utility_exit_with_message("dbg tgt res hash");
	if( ! utility::file::file_exists( opts.data_dir) ){
		if( !opts.fail_if_no_cached_data ){
			utility::file::create_directory_recursive(opts.data_dir);
			if( ! utility::file::file_exists( opts.data_dir) ){
				utility_exit_with_message("missing data dir: '" + opts.data_dir + "'" );
			}
		}
	}

	std::string cache_prefix = opts.data_dir+"/__RF_"+target_tag
		             +"_trhash" + hashstr
	                 +"_resl"   + boost::lexical_cast<std::string>(opts.field_resl)
		             +"_osamp"  + boost::lexical_cast<std::string>(opts.oversample);
	if( opts.block_hbond_sites       ) cache_prefix += "_blockhb";
	if( opts.repulsive_only_boundary ) cache_prefix += "_replonlybdry";
	if( opts.repulsive_atoms.size() ){
		cache_prefix += "_repl";
		for( auto aid : opts.repulsive_atoms ){
			cache_prefix += "_"+str(aid.atomno())+"."+str(aid.rsd());
		}
	}

	return get_rosetta_fields_specified_cache_prefix(
		cache_prefix,
		target_fname,
		target,
		target_res,
		opts,
		field_by_atype,
		verbose
	);
}



std::string
get_rosetta_fields_specified_cache_prefix(
	std::string const & cache_prefix,
	std::string const & ,//target_fname,
	core::pose::Pose const & target,
	utility::vector1<core::Size> const & target_res,
	RosettaFieldOptions const & opts,
	std::vector<
		// shared_ptr<
			::scheme::objective::voxel::VoxelArray<3,float> *
		// >
	> & field_by_atype,
	bool verbose
){
  	using ObjexxFCL::format::I;

	field_by_atype.resize(N_ATYPE+1,nullptr);

	// typedef Eigen::Transform<double,3,Eigen::AffineCompact> EigenXform;
	// load 1.25, 0-162m, 1-208m, 2-294m, 3-453m
	// load 3.67, 0-56, 1-69m


	typedef ::scheme::actor::Atom< Eigen::Vector3f > SchemeAtom;
	// runtime_assert( sizeof(SchemeAtom) == 16 );
	typedef ::scheme::rosetta::score::RosettaField<SchemeAtom,devel::scheme::EtableParamsInit> RosettaField;

	typedef ::scheme::util::SimpleArray<3,float> F3;

	typedef ::scheme::objective::voxel::FieldCache3D<float> FieldCache;
	typedef ::scheme::objective::voxel::VoxelArray<3,float> VoxelArray;

	float field_resl = opts.field_resl;
	int oversample = opts.oversample;

		std::vector<int> atypemap = get_rif_atype_map();

		std::vector<SchemeAtom> target_atoms;
		{
			std::vector<SchemeAtom> selected_atoms;
			get_scheme_atoms( target, target_res, selected_atoms );
			// for(int ires = 1; ires <= target_res.size(); ++ires){
			// 	int ir = target_res[ires];
			// 	core::conformation::Residue const & r( target.residue(ir) );
			// 	for(int ia = 1; ia <= r.nheavyatoms(); ++ia){
			// 		if( r.is_virtual(ia) ) continue;
			// 		SchemeAtom a( r.xyz(ia), atypemap[ r.atom_type_index(ia) ] );
			// 		selected_atoms.push_back(a);
			// 	}
			// }
			for(int ir = 1; ir <= target.size(); ++ir){
				core::conformation::Residue const & r( target.residue(ir) );
				for(int ia = 1; ia <= r.nheavyatoms(); ++ia){
					if( r.is_virtual(ia) ) continue;
					int at = atypemap[ r.atom_type_index(ia) ];
					if( at > N_ATYPE ){
						// utility_exit_with_message("heavy atom type > N_ATYPE: "+str(at)+" "+r.name()+" "+r.atom_name(ia) );
						std::cout << "WARNING: heavy atom type"<<r.atom_type_index(ia)<<" > N_ATYPE: "+str(at)+" "+r.name()+" "+r.atom_name(ia) 
						          << " will treat as carbon for sterics!" << std::endl;
						at = 5;
					}

					bool is_selected = ( std::find(target_res.begin(),target_res.end(),ir) != target_res.end() );
					bool is_repulser = ( std::find(opts.repulsive_atoms.begin(),opts.repulsive_atoms.end(),core::id::AtomID(ia,ir)) != opts.repulsive_atoms.end() );
					if( is_repulser ) std::cout << "repulser atom: " << ir << " " << target.residue(ir).name3() << " " << ia << " " << target.residue(ir).atom_name(ia) << std::endl;

					if( is_selected ){

						int16_t atrepl( -12345 ); // can fit in int16_t
						SchemeAtom a( r.xyz(ia), is_repulser? atrepl : at );
						target_atoms.push_back(a);

					} else {

						int at_boundary = at;
						if( opts.repulsive_only_boundary ) at_boundary = -at;
						SchemeAtom a( r.xyz(ia), at_boundary ); // -at for neg-only
						bool close_to_selected = false;
						for( auto const & b : selected_atoms ){
							float const dx = b.position()[0]-a.position()[0];
							float const dy = b.position()[1]-a.position()[1];
							float const dz = b.position()[2]-a.position()[2];
							float const dis2 = dx*dx+dy*dy+dz*dz;
							if( dis2 < 225.0 ) close_to_selected = true;
						}
						if( close_to_selected ) target_atoms.push_back(a);
					}

				}

				if( opts.block_hbond_sites ){
					for(int i = 1; i <= r.Hpos_polar().size(); ++i){
						int ia = r.Hpos_polar()[i];
						// std::cout << "block_hbond_sites" << std::endl;
						int ib = r.atom_base(ia);
						target_atoms.push_back( SchemeAtom( r.xyz(ia), -5 ) );
						target_atoms.push_back( SchemeAtom( r.xyz(ia) + (r.xyz(ia)-r.xyz(ib)).normalized()*2.8, -5 ) );
					}
				}

			}

			if( target_atoms.size() == 0 ) utility_exit_with_message("no target atoms!");
		}
		RosettaField rosetta_field( target_atoms );

		F3 lb(9e9,9e9,9e9),ub(-9e9,-9e9,-9e9);
		uint64_t target_atoms_hash = 0;
		for( auto const & a : target_atoms ){
			lb = lb.min(a.position());
			ub = ub.max(a.position());
		}
		std::cout << "rosetta_field lb: " << lb << " ub: " << ub << " size(A): " << ub-lb << std::endl;

		std::exception_ptr exception = nullptr;
		#ifdef USE_OPENMP
		#pragma omp parallel for schedule(dynamic,1)
		#endif
		for( int itype = 1; itype <= N_ATYPE; ++itype ){
            if( exception ) continue;
            if( opts.one_atype_only && itype != opts.one_atype_only ) continue;
			try {
				std::string cachefile = cache_prefix +"__atype"+boost::lexical_cast<std::string>(itype)+".rosetta_field.gz";

				::scheme::rosetta::score::RosettaFieldAtype< SchemeAtom, devel::scheme::EtableParamsInit > rfa( rosetta_field, itype );

				if( utility::file::file_exists(cachefile) ){
					if( opts.generate_only ) continue;
					if( verbose ){
						#ifdef USE_OPENMP
						#pragma omp critical
						#endif
						std::cout<< "thread " << I(3,omp_thread_num_1()) << " init  rosetta_field " << I(2,itype) << " CACHE AT " << cachefile << std::endl;
					}
					// field_by_atype[itype] = boost::make_shared< FieldCache >( rfa, lb-6.0f, ub+6.0f, field_resl, "", true, oversample ); // no init
					field_by_atype[itype] = new FieldCache( rfa, lb-6.0f, ub+6.0f, field_resl, "", true, oversample ); // no init
					utility::io::izstream in( cachefile, std::ios::binary );
					field_by_atype[itype]->load(in);
					in.close();
				} else {
					if( opts.fail_if_no_cached_data ){
						#ifdef USE_OPENMP
						#pragma omp critical
						#endif
						{
							std::cout << "fail_if_no_cached_data set, and data not available: " << std::endl;
							std::cout << cachefile << std::endl;
							utility_exit_with_message("required data not available");
						}
					}
					#ifdef USE_OPENMP
					#pragma omp critical
					#endif
					std::cout << "thread " << I(3,omp_thread_num_1()) << " init  rosetta_field " << I(2,itype) << " CACHE TO " << cachefile << std::endl;
					// field_by_atype[itype] = boost::make_shared<FieldCache >( rfa, lb-6.0f, ub+6.0f, field_resl, "", false, oversample );
					field_by_atype[itype] = new FieldCache( rfa, lb-6.0f, ub+6.0f, field_resl, "", false, oversample );
					utility::io::ozstream out( cachefile , std::ios::binary );
					field_by_atype[itype]->save( out );
					out.close();
				}
				// if( opts.cache_mismatch_tolerance < 9e8 ){
				// 	double erf = static_cast<FieldCache&>(*field_by_atype[itype]).check_against_field( rfa, oversample, opts.cache_mismatch_tolerance );
				// 	if( erf > 0.0 ){
				// 		#ifdef USE_OPENMP
				// 		#pragma omp critical
				// 		#endif
				// 		{
				// 			cout << "FIELD MISMATCH ERROR ATYPE " << itype << " error_frac: " << erf << endl;
				// 			cout << "ERROR ON " << cachefile << endl;
				// 			utility_exit_with_message("field cache errors!");
				// 		}
				// 	}
				// }
			} catch( ... ) {
				#ifdef USE_OPENMP
				#pragma omp critical
				#endif
				exception = std::current_exception();
			}
		}
		if( exception ) std::rethrow_exception(exception);


		return cache_prefix;

}






std::string
get_rosetta_bounding_fields(
	std::vector<float> const & RESLS,
	std::string const & target_fname,
	core::pose::Pose const & target,
	utility::vector1<core::Size> const & target_res,
	RosettaFieldOptions const & opts,
	std::vector<              ::scheme::objective::voxel::VoxelArray<3,float> *   > & field_by_atype,
	std::vector< std::vector< ::scheme::objective::voxel::VoxelArray<3,float> *	> > & bounding_by_atype,
	bool verbose
){
	typedef ::scheme::objective::voxel::VoxelArray<3,float> VoxelArray;
	typedef ::scheme::objective::voxel::BoundingFieldCache3D<float> BoundingGrid;

	using ObjexxFCL::format::I;

	std::string cache_prefix = devel::scheme::get_rosetta_fields(
		target_fname,
		target,
		target_res,
		opts,
		field_by_atype,
		verbose
	);

	get_rosetta_bounding_fields_from_fba(
		RESLS,
		target_fname,
		target,
		target_res,
		opts,
		field_by_atype,
		bounding_by_atype,
		verbose,
		cache_prefix
	);

	return cache_prefix;

}

void
get_rosetta_bounding_fields_from_fba(
	std::vector<float> const & RESLS,
	std::string const & target_fname,
	core::pose::Pose const & target,
	utility::vector1<core::Size> const & target_res,
	RosettaFieldOptions const & opts,
	std::vector<              ::scheme::objective::voxel::VoxelArray<3,float> *   > & field_by_atype,
	std::vector< std::vector< ::scheme::objective::voxel::VoxelArray<3,float> *	> > & bounding_by_atype,
	bool verbose,
	std::string cache_prefix
){

	typedef ::scheme::objective::voxel::VoxelArray<3,float> VoxelArray;
	typedef ::scheme::objective::voxel::BoundingFieldCache3D<float> BoundingGrid;

	using ObjexxFCL::format::I;

	std::vector<std::pair<float,int> > jobs;
	for( int iresl = RESLS.size()-1; iresl >= 0; --iresl ){
		for( int itype = 1; itype <= N_ATYPE; ++itype ){
            if( opts.one_atype_only && itype != opts.one_atype_only ) continue;
			jobs.push_back( std::make_pair( iresl, itype ) );
		}
	}
	bounding_by_atype.resize( RESLS.size() );
	for(int i = 0; i < RESLS.size(); ++i) bounding_by_atype[i].resize(25,nullptr);
	std::exception_ptr exception = nullptr;
	#ifdef USE_OPENMP
	#pragma omp parallel for schedule(dynamic,1)
	#endif
	for( int ijob = 0; ijob < jobs.size(); ++ijob ){
		if(exception) continue;
		try {
			int const iresl = jobs[ijob].first;
			int const itype = jobs[ijob].second;
			float const bound = RESLS[iresl];
			float bresl = std::max<float>( bound/opts.max_bounding_ratio, opts.field_resl );
					// bresl = std::min( 0.25f, bresl );
			std::string cachefile = cache_prefix
				+"_bounding"+boost::lexical_cast<std::string>(bound)+"_"+boost::lexical_cast<std::string>(bresl)
				+"_atype" + boost::lexical_cast<std::string>(itype)
				 +".rf.gz";
			VoxelArray * gp;
			if( utility::file::file_exists(cachefile) ){
				if( opts.generate_only ) continue;
				if(verbose||itype==1){
					#ifdef USE_OPENMP
					#pragma omp critical
					#endif
					std::cout << "thread " << I(3,omp_thread_num_1()) << " init bounding field " << I(2,iresl) << " " << I(2,itype) << " CACHE AT " << cachefile << std::endl;
				}
					// gp = boost::make_shared<BoundingGrid>( *field_by_atype[itype], bound, bresl, "", true );
				gp = new BoundingGrid( *field_by_atype[itype], bound, bresl, "", true );
				utility::io::izstream in( cachefile, std::ios::binary );
				gp->load(in);
				in.close();
			} else {
				#ifdef USE_OPENMP
				#pragma omp critical
				#endif
				std::cout << "thread " << I(3,omp_thread_num_1()) << " init bounding field " << I(2,itype) << " CACHE TO " << cachefile << std::endl;
					// gp = boost::make_shared< BoundingGrid >( *field_by_atype[itype], bound, bresl, "", false );
				if( bound/bresl < 1.1 ){
					#ifdef USE_OPENMP
					#pragma omp critical
					#endif
					std::cout << "WARNING: bound/resl: " << bound << "/" << bresl << " too small, using unmodified source grid" << std::endl;
					gp = field_by_atype[itype];
				} else {
					gp = new BoundingGrid( *field_by_atype[itype], bound, bresl, "", false );
				}
				utility::io::ozstream out( cachefile , std::ios::binary );
				gp->save( out );
				out.close();
			}
			bounding_by_atype.at(iresl).at(itype) = gp;
		} catch( ... ) {
			#ifdef USE_OPENMP
			#pragma omp critical
			#endif
			exception = std::current_exception();
		}
	}
	if( exception ) std::rethrow_exception(exception);


}


}}

