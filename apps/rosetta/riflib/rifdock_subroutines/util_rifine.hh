// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:


#ifndef INCLUDED_riflib_rifdock_subroutines_util_hh
#define INCLUDED_riflib_rifdock_subroutines_util_hh

#include <ObjexxFCL/format.hh>

#include <riflib/types.hh>

#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>
#include<core/scoring/rms_util.hh>

#include <numeric/xyzTransform.hh>
#include <numeric/xyzVector.hh>
#include <utility/vector1.hh>
#include <scheme/kinematics/Director.hh>
#include <scheme/search/HackPack.hh>
#include <riflib/util.hh>

#include <rifdock_legacy.hh>
#include <riflib/rotamer_energy_tables.hh>
#include <riflib/RifFactory.hh>


#include <utility/file/file_sys_util.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>


using ::scheme::make_shared;
using ::scheme::shared_ptr;

typedef int32_t intRot;

namespace devel {
namespace scheme {

inline
Eigen::Vector3f
pose_center(
    core::pose::Pose const & pose,
    utility::vector1<core::Size> const & useres = utility::vector1<core::Size>()
){
    typedef numeric::xyzVector<core::Real> Vec;
    Vec cen(0,0,0);
    int count = 0;
    for( int ir = 1; ir <= pose.size(); ++ir ) {
        if( useres.size()==0 || std::find(useres.begin(),useres.end(),ir)!=useres.end() ){
            for( int ia = 1; ia <= pose.residue(ir).nheavyatoms(); ++ia ){
                cen += pose.xyz(core::id::AtomID(ia,ir));
                ++count;
            }
        // } else {
            // std::cout << "pose_center skip " << ir << std::endl;
        }
    }
    cen /= double(count);
    // ::scheme::util::SimpleArray<3,float> center;
    Eigen::Vector3f center;
    center[0] = cen[0];
    center[1] = cen[1];
    center[2] = cen[2];
    return center;
}

inline
void
get_rg_radius(
    core::pose::Pose const & pose,
    float & rg,
    float & radius,
    utility::vector1<core::Size> const & useres = utility::vector1<core::Size>(),
    bool allatom = false
){
    Eigen::Vector3f centmp = pose_center( pose, useres );
    numeric::xyzVector<double> cen;
    float maxdis = -9e9, avgdis2 = 0.0;
    for( int i = 0; i < 3; ++i ) cen[i] = centmp[i];
    for( int iri = 1; iri <= useres.size(); ++iri ){
        int ir = useres[iri];
        if( allatom ){
            for( int ia = 1; ia <= pose.residue(ir).nheavyatoms(); ++ia ){
                numeric::xyzVector<double> coord = pose.residue(ir).xyz(ia);
                avgdis2 += cen.distance_squared( coord );
                maxdis = std::max( maxdis, (float)cen.distance( coord ) );
            }
        } else {
            numeric::xyzVector<double> coord;
            if(      pose.residue(ir).has("CB") ) coord = pose.residue(ir).xyz("CB");
            else if( pose.residue(ir).has("CA") ) coord = pose.residue(ir).xyz("CA");
            else                                  coord = pose.residue(ir).nbr_atom_xyz();
            avgdis2 += cen.distance_squared( coord );
            maxdis = std::max( maxdis, (float)cen.distance( coord ) );
        }
    }
    avgdis2 /= useres.size();
    rg = sqrt( avgdis2 );
    radius = maxdis;
}


inline
void xform_pose( core::pose::Pose & pose, numeric::xyzTransform<float> s, core::Size sres=1, core::Size eres=0 ) {
  if(eres==0) eres = pose.size();
  for(core::Size ir = sres; ir <= eres; ++ir) {
    for(core::Size ia = 1; ia <= pose.residue(ir).natoms(); ++ia) {
      core::id::AtomID const aid(core::id::AtomID(ia,ir));
      pose.set_xyz( aid, s*pose.xyz(aid) );
    }
  }
}

    
void
add_pdbinfo_if_missing( core::pose::Pose & pose ) {
    if ( ! pose.pdb_info() ) {
        core::pose::PDBInfoOP pdb_info = make_shared<core::pose::PDBInfo>( pose );
        pose.pdb_info(pdb_info);
    }
}
    
    

bool parse_seeding_file(std::string fname, std::vector<devel::scheme::EigenXform> & seeding_positions, bool seeding_by_patchdock = true, float patchdock_min_sasa = -1000.0, int patchdock_top_ranks = 99999)
{
        
    // the seeding
    
    if (fname == "" ) {
        return false;
    }
    runtime_assert_msg(utility::file::file_exists( fname ), "seeding file does not exits: " + fname );
    
    std::ifstream in;
    in.open(fname, std::ios::in);
    // utility::io::izstream in(fname);
		std::string s;
		devel::scheme::EigenXform xform;
		if ( seeding_by_patchdock ){

				bool flag = false;
				while ( std::getline(in, s)){
						if (s.empty()) continue;

						utility::vector1<std::string> splt = utility::quoted_split( s );
						
						if (!flag && splt[1].find("#") == 0  && splt[3].find("score") == 0 ) {
								flag = true;
								continue;
						}

						if(!flag) continue;

						// remove bad patchdock seeding pos based on the sasa
						if ( utility::string2float(splt[7]) < patchdock_min_sasa ) continue;
						if ( utility::string2int  (splt[1]) > patchdock_top_ranks) continue;

						float cx = cos(utility::string2float(splt[25]));
						float cy = cos(utility::string2float(splt[26]));
						float cz = cos(utility::string2float(splt[27]));
						float sx = sin(utility::string2float(splt[25]));
						float sy = sin(utility::string2float(splt[26]));
						float sz = sin(utility::string2float(splt[27]));
						float tx = utility::string2float(splt[28]);
						float ty = utility::string2float(splt[29]);
						float tz = utility::string2float(splt[30]);
						
						xform.linear().row(0) = Eigen::Vector3f( cz * cy, -sy * sx * cz - sz * cx, -sy * cx *cz + sz *sx);
						xform.linear().row(1) = Eigen::Vector3f( sz * cy, -sy * sx * sz + cx * cz, -sy * cx * sz - sx * cz);
						xform.linear().row(2) = Eigen::Vector3f( sy, cy * sx, cy * cx);
						xform.translation() = Eigen::Vector3f( tx, ty, tz );

						seeding_positions.push_back(xform);
				}
		} else {
				// the format of the seeding file is Rosetta Xforms. the easy case
				while (std::getline(in, s)) {
						if (s.empty()) continue;
						
						utility::vector1<std::string> splt = utility::quoted_split( s );
						if (splt[1].find("#") == 0  || splt.size() == 0) {
								continue;
						}
        
						xform.linear().row(0) = Eigen::Vector3f(utility::string2float(splt[1]), utility::string2float(splt[2]), utility::string2float(splt[3]) );
						xform.linear().row(1) = Eigen::Vector3f(utility::string2float(splt[4]), utility::string2float(splt[5]), utility::string2float(splt[6]) );
						xform.linear().row(2) = Eigen::Vector3f(utility::string2float(splt[7]), utility::string2float(splt[8]), utility::string2float(splt[9]) );
						xform.translation() = Eigen::Vector3f(utility::string2float(splt[10]), utility::string2float(splt[11]), utility::string2float(splt[12]) );
        
						seeding_positions.push_back(xform);
				}
    }
        
    return true;
}
    

bool parse_exhausitive_searching_file(std::string fname, std::vector<std::pair< int64_t, devel::scheme::EigenXform > > & searching_positions, double maximum_ang = 999 )
{
        
    // the seeding
        
    if (fname == "" ) {
        return false;
    }
    runtime_assert_msg(utility::file::file_exists( fname ), "exhausitive searching file does not exits: " + fname );
    
    std::ifstream in;
    in.open(fname, std::ios::in);
    // utility::io::izstream in(fname);
    std::string s;
    devel::scheme::EigenXform xform;
    while (std::getline(in, s)) {
        if (s.empty()) continue;
            
        utility::vector1<std::string> splt = utility::string_split_simple(s, ' ');
        if (splt[1].find("#") == 0  || splt.size() == 0) {
            continue;
        }
        
        if ( std::abs( utility::string2float(splt[3]) ) > maximum_ang ) {
            continue;
        }
        
        xform.linear().row(0) = Eigen::Vector3f(utility::string2float(splt[4]), utility::string2float(splt[5]), utility::string2float(splt[6]) );
        xform.linear().row(1) = Eigen::Vector3f(utility::string2float(splt[7]), utility::string2float(splt[8]), utility::string2float(splt[9]) );
        xform.linear().row(2) = Eigen::Vector3f(utility::string2float(splt[10]), utility::string2float(splt[11]), utility::string2float(splt[12]) );
        xform.translation() = Eigen::Vector3f(utility::string2float(splt[13]), utility::string2float(splt[14]), utility::string2float(splt[15]) );
        
        searching_positions.push_back( std::pair< int64_t, devel::scheme::EigenXform >( utility::string2int(splt[2]), xform) );
    }
        
    return true;
}


template<class _DirectorBigIndex>
struct tmplRifDockResult {
    typedef _DirectorBigIndex DirectorBigIndex;
    typedef tmplRifDockResult<DirectorBigIndex> This;
    float dist0, packscore, nopackscore, rifscore, stericscore;
    uint64_t isamp;
    DirectorBigIndex scene_index;
    uint32_t prepack_rank;
    float cluster_score;
    bool operator< ( This const & o ) const { return packscore < o.packscore; }
    shared_ptr< std::vector< std::pair<intRot,intRot> > > rotamers_;
    core::pose::PoseOP pose_ = nullptr;
    size_t numrots() const { if(rotamers_==nullptr) return 0; return rotamers_->size(); }
    std::vector< std::pair<intRot,intRot> > const & rotamers() const { assert(rotamers_!=nullptr); return *rotamers_; }
};
    

#pragma pack (push, 4) // allows size to be 12 rather than 16
template<class _DirectorBigIndex>
struct tmplSearchPoint {
    typedef _DirectorBigIndex DirectorBigIndex;
    typedef tmplSearchPoint<DirectorBigIndex> This;
    float score;
    DirectorBigIndex index;
    tmplSearchPoint() : score(9e9) {}
    tmplSearchPoint(DirectorBigIndex i) : score(9e9), index(i) {}
    bool operator < (This const & o) const {
        return score < o.score;
    }
};
#pragma pack (pop)
    
template <class __Director>
using _SearchPoint = tmplSearchPoint<_DirectorBigIndex<__Director>>;
typedef _SearchPoint<DirectorBase> SearchPoint;
    
    
    
    
    
template<class _DirectorBigIndex>
struct tmplSearchPointWithRots {
    typedef _DirectorBigIndex DirectorBigIndex;
    typedef tmplSearchPointWithRots<DirectorBigIndex> This;
    float score;
    uint32_t prepack_rank;
    DirectorBigIndex index;
    shared_ptr< std::vector< std::pair<intRot,intRot> > > rotamers_;
    core::pose::PoseOP pose_ = nullptr;
    tmplSearchPointWithRots() : score(9e9), prepack_rank(0), rotamers_(nullptr) {}
    tmplSearchPointWithRots(DirectorBigIndex i, uint32_t orank) : score(9e9), prepack_rank(orank), index(i), rotamers_(nullptr) {}
    // ~SearchPointWithRots() { delete rotamers_; }
    void checkinit() { if( rotamers_==nullptr ) rotamers_ = make_shared< std::vector< std::pair<intRot,intRot> > > ();  }
    std::vector< std::pair<intRot,intRot> > & rotamers() { checkinit(); return *rotamers_; }
    std::vector< std::pair<intRot,intRot> > const & rotamers() const { runtime_assert(rotamers_!=nullptr); return *rotamers_; }
    size_t numrots() const { if(rotamers_==nullptr) return 0; return rotamers_->size(); }
    bool operator < (This const & o) const {
        return score < o.score;
    }
};

    
template <class __Director>
using _SearchPointWithRots = tmplSearchPointWithRots<_DirectorBigIndex<__Director>>;
typedef _SearchPointWithRots<DirectorBase> SearchPointWithRots;

template <class __Director>
using _RifDockResult = tmplRifDockResult<_DirectorBigIndex<__Director>>;
typedef _RifDockResult<DirectorBase> RifDockResult;
    


    // This function is just a simplified version of awful_compile_output_helper
    // designed to do redundancy filtering after hackpack
    // the only function is redundancy filtering, no score check, no pose check.
    // purely based on xform position.
    template<
    class EigenXform,
    class ScenePtr,
    class ObjectivePtr
    >
    void real_rmsd_xform(
                              int iresl,
                              std::vector< SearchPointWithRots > const & packed_results,
                              std::vector< ScenePtr > & scene_pt,
                              DirectorBase director,
                              float redundancy_filter_rg,
															core::pose::Pose const & pose
    )
    {
				const int64_t total_samples = packed_results.size();
				std::cout << std::endl << "Now calculating the real rmsd and the xform_magnitute (Very slow, I didn't optimize this. Only for testing. for " << total_samples << std::endl << std::endl;;
				const int64_t total_num_scores = total_samples * (total_samples - 1) / 2;
				const int64_t samples_per_dot = total_samples / 109>1? total_samples/109 : 1;
				typedef std::pair<float, float> RMSD_XFORM;
				std::vector<RMSD_XFORM> results(total_num_scores);

				#ifdef USE_OPENMP
        #pragma omp parallel for schedule(dynamic,16)
        #endif
				for( int64_t isamp = 0; isamp < total_samples; ++isamp ) {
						if ( isamp % samples_per_dot == 0 ) std::cout << "*" << std::flush;
						SearchPointWithRots const & isp = packed_results[isamp];
						ScenePtr scene_minimal( scene_pt[omp_get_thread_num()] );
						director->set_scene( isp.index, iresl, *scene_minimal );
						EigenXform xposition1 = scene_minimal->position(1);
						EigenXform xposition1inv = xposition1.inverse();
						core::pose::Pose ipose(pose);
						core::pose::Pose jpose(pose);
						xform_pose(ipose, eigen2xyz(xposition1));
						for( int64_t jsamp = isamp+1; jsamp<total_samples; ++jsamp) {
								SearchPointWithRots const & jsp = packed_results[jsamp];
								director->set_scene( jsp.index, iresl, *scene_minimal );
								EigenXform xposition2 = scene_minimal->position(1);
								EigenXform xposition2inv = xposition2.inverse();
								EigenXform const xdiff = xposition1inv * xposition2;
								float xform_rmsd = devel::scheme::xform_magnitude(xdiff, redundancy_filter_rg);
								xform_pose(jpose, eigen2xyz(xposition2));
								float real_rmsd  = core::scoring::all_atom_rmsd_nosuper(ipose, jpose);
								xform_pose(jpose, eigen2xyz(xposition2inv));
								int64_t index = total_num_scores - (total_samples - isamp) * (total_samples - isamp - 1) / 2 + (jsamp - isamp) - 1;
								results[index].first = real_rmsd;
								results[index].second = xform_rmsd;
						}
						
				}
				
				utility::io::ozstream results_file( "rmsd_xformmag.list" );
				for ( RMSD_XFORM const & r : results ) results_file << r.first << " " << r.second << std::endl;
				results_file.close();

				std::ostringstream oss;
				oss << "#isamp seeding_index nest_index xx xy xz yx yy yz zx zy zz tx ty tz\n";
				for ( int64_t isamp = 0; isamp < total_samples; ++ isamp ) {
						SearchPointWithRots const & isp = packed_results[isamp];
						ScenePtr scene_minimal( scene_pt.back() );
						director->set_scene( isp.index, iresl, *scene_minimal );
						EigenXform xposition1 = scene_minimal->position(1);
						// now dump the top100 index, xform, ... for debuging
						oss << isamp << " " << isp.index.seeding_index << " " << isp.index.nest_index << " "
								<< xposition1.linear().row(0) << " " << xposition1.linear().row(1) << " " << xposition1.linear().row(2) << " "
								<< xposition1.translation().x() << " " << xposition1.translation().y() << " " << xposition1.translation().z() << std::endl;
				}
				utility::io::ozstream info_file("results_info.list");
				info_file << oss.str();
				info_file.close();

				return;
    }
    
    
    // This function is just a simplified version of awful_compile_output_helper
    // designed to do redundancy filtering after hackpack
    // the only function is redundancy filtering, no score check, no pose check.
    // purely based on xform position.
    template<
    class EigenXform,
    class ScenePtr,
    class ObjectivePtr
    >
    void redundancy_filtering(
                              int64_t isamp,
                              int iresl,
                              std::vector< SearchPointWithRots > const & packed_results,
                              std::vector< ScenePtr > & scene_pt,
                              DirectorBase director,
                              float redundancy_filter_rg,
                              float redundancy_filter_mag,
                              std::vector< SearchPointWithRots > & selected_results,
                              std::vector< EigenXform > & selected_xforms,
                              #ifdef USE_OPENMP
                              omp_lock_t & dump_lock
                              #endif
    )
    {
        SearchPointWithRots const & sp = packed_results[isamp];
        ScenePtr scene_minimal( scene_pt[omp_get_thread_num()] );
        director->set_scene( sp.index, iresl, *scene_minimal );
        EigenXform xposition1 = scene_minimal->position(1);
        EigenXform xposition1inv = xposition1.inverse();
        
				bool is_redundant = false;
        BOOST_FOREACH( EigenXform const & xsel, selected_xforms ){
            EigenXform const xdiff = xposition1inv * xsel;
            if( devel::scheme::xform_magnitude_redundancy_check( xdiff, redundancy_filter_rg, redundancy_filter_mag) ) {
								is_redundant = true;
								break;
						}
        }
        if( ! is_redundant ) {
            // need to double check this because of multi-threading issues
            #ifdef USE_OPENMP
            omp_set_lock( &dump_lock );
            #endif
            {
                bool is_not_redundant = true;
                BOOST_FOREACH( EigenXform const & xsel, selected_xforms ){
                    EigenXform const xdiff = xposition1inv * xsel;
										if( devel::scheme::xform_magnitude_redundancy_check( xdiff, redundancy_filter_rg, redundancy_filter_mag) ) {
												is_not_redundant = false;
												break;
										}
                }
                if ( is_not_redundant ) {
                    if( redundancy_filter_mag > 0.0001 ) {
                        selected_xforms.push_back( xposition1 );
                    }
                    selected_results.push_back( sp );
                }
            }
            #ifdef USE_OPENMP
            omp_unset_lock( &dump_lock );
            #endif
        }
    }
    
    // how can I fix this??? make the whole prototype into a class maybe???
    // what does it do?
    //  set and rescore scene with nopackscore, record more score detail
    //  compute dist0
    //  select results with some redundancy filtering
    template<
    class EigenXform,
    // class Scene,
    class ScenePtr,
    class ObjectivePtr
    >
    void
    awful_compile_output_helper(
                                int64_t isamp,
                                int iresl,
                                std::vector< SearchPointWithRots > const & packed_results,
                                std::vector< ScenePtr > & scene_pt,
                                DirectorBase director,
                                float redundancy_filter_rg,
                                float redundancy_filter_mag,
                                Eigen::Vector3f scaffold_center,
                                std::vector< std::vector< RifDockResult > > & allresults_pt,
                                std::vector< RifDockResult >   & selected_results,
                                std::vector< std::pair<EigenXform, uint64_t> > & selected_xforms,
                                int n_pdb_out,
                                #ifdef USE_OPENMP
                                omp_lock_t & dump_lock,
                                #endif
                                ObjectivePtr objective,
                                int & nclose,
                                int nclosemax,
                                float nclosethresh,
                                EigenXform scaffold_perturb
                                ) {
        
        SearchPointWithRots const & sp = packed_results[isamp];
        if( sp.score >= 0.0f ) return;
				// if the pose is null, return. Longxing, 20190323
				if ( sp.pose_ == nullptr ) return;
        ScenePtr scene_minimal( scene_pt[omp_get_thread_num()] );
        director->set_scene( sp.index, iresl, *scene_minimal );
        std::vector<float> sc = objective->scores(*scene_minimal);
        float const nopackscore = sc[0]+sc[1]; //result.sum();
        float const rifscore = sc[0]; //result.template get<MyScoreBBActorRIF>();
        float const stericscore = sc[1]; //result.template get<MyClashScore>();
        
        // dist0 is only important to the nclose* options
        float dist0; {
            EigenXform x = scene_minimal->position(1);
            x = scaffold_perturb * x;
            x.translation() -= scaffold_center;
            dist0 = ::devel::scheme::xform_magnitude( x, redundancy_filter_rg );
        }
        
        RifDockResult r; // float dist0, packscore, nopackscore, rifscore, stericscore;
        r.isamp = isamp;
        r.prepack_rank = sp.prepack_rank;
        r.scene_index = sp.index;
        r.packscore = sp.score;
        r.nopackscore = nopackscore;
        r.rifscore = rifscore;
        r.stericscore = stericscore;
        r.dist0 = dist0;
        r.cluster_score = 0.0;
        r.pose_ = sp.pose_;
        allresults_pt.at( omp_get_thread_num() ).push_back( r ); // recorded w/o rotamers here
        
        bool force_selected = ( dist0 < nclosethresh && ++nclose < nclosemax ); // not thread-safe... is this important?
        
        if( selected_xforms.size() < n_pdb_out || force_selected ){
            
            EigenXform xposition1 = scene_minimal->position(1);
            EigenXform xposition1inv = xposition1.inverse();
            
            float mindiff_candidate = 9e9;
            int64_t i_closest_result;
            typedef std::pair< EigenXform, int64_t > XRpair;
            BOOST_FOREACH( XRpair const & xrp, selected_xforms ){
                EigenXform const & xsel = xrp.first;
                EigenXform const xdiff = xposition1inv * xsel;
                float diff = devel::scheme::xform_magnitude( xdiff, redundancy_filter_rg );
                if( diff < mindiff_candidate ){
                    mindiff_candidate = diff;
                    i_closest_result = xrp.second;
                }
                // todo: also compare AA composition of rotamers
            }
            
            if( mindiff_candidate < redundancy_filter_mag ){ // redundant result
                selected_results[i_closest_result].cluster_score += 1.0; //sp.score==0.0 ? nopackscore : sp.score;
            }
            
            if( mindiff_candidate > redundancy_filter_mag || force_selected ){
                
                #ifdef USE_OPENMP
                omp_set_lock( &dump_lock );
                #endif
                {
                    // std::cout << "checking again to add selected " << selected_xforms.size() << " " << omp_get_thread_num() << std::endl;
                    float mindiff_actual = 9e9;
                    BOOST_FOREACH( XRpair const & xrp, selected_xforms ){
                        EigenXform const & xsel = xrp.first;
                        EigenXform const xdiff = xposition1inv * xsel;
                        mindiff_actual = std::min( mindiff_actual, devel::scheme::xform_magnitude( xdiff, redundancy_filter_rg ) );
                    }
                    if( mindiff_actual > redundancy_filter_mag || force_selected ){
                        if( redundancy_filter_mag > 0.0001 ) {
                            selected_xforms.push_back( std::make_pair(xposition1, (int64_t)selected_results.size()) );
                        }
                        r.rotamers_ = sp.rotamers_;
                        selected_results.push_back( r ); // recorded with rotamers here
                    } else {
                        // std::cout << " second check failed" << std::endl;
                    }
                }
                #ifdef USE_OPENMP
                omp_unset_lock( &dump_lock );
                #endif
                
            } // end if( mindiff > redundancy_filter_mag ){
            
        } // end    if( selected_xforms.size() < n_pdb_out || force_selected )
        
    }
    
    
    // the reason that I didn't pass the reference of the option struct is that "too many different option structs"
    // To be safe.
    void
    output_results(
                   std::vector< RifDockResult > & selected_results,
                   std::vector< shared_ptr< RifBase> > & rif_ptrs,
                   DirectorBase director,
                   ScenePtr scene_minimal,
                   int iresl,
                   std::vector<int> & scaffres_l2g,
                   bool align_to_scaffold,
                   std::string scafftag,
                   std::string outdir,
                   std::string output_tag,
                   
                   utility::io::ozstream & dokout) {

				using ObjexxFCL::format::F;
				using ObjexxFCL::format::I;

        if( align_to_scaffold ) std::cout << "ALIGN TO SCAFFOLD" << std::endl;
        else                        std::cout << "ALIGN TO TARGET"   << std::endl;
        
        for( int i_selected_result=0; i_selected_result < selected_results.size(); ++i_selected_result ){
            RifDockResult const & selected_result = selected_results.at(i_selected_result);
            std::string pdboutfile = outdir + "/" + scafftag + "_" + devel::scheme::str(i_selected_result,9)+".pdb.gz";
            if( output_tag.size() ){
                pdboutfile = outdir + "/" + scafftag+"_" + output_tag + "_" + devel::scheme::str(i_selected_result,9)+".pdb.gz";
            }
            //std::string resfileoutfile = rdd.opt.outdir + "/" + scafftag+"_"+devel::scheme::str(i_selected_result,9)+".resfile";
            //std::string allrifrotsoutfile = rdd.opt.outdir + "/" + scafftag+"_allrifrots_"+devel::scheme::str(i_selected_result,9)+".pdb.gz";
            
            std::ostringstream oss;
            oss << "rif score: " << I(4,i_selected_result)
            << " rank "       << I(9,selected_result.isamp)
            << " dist0:    "  << F(7,2,selected_result.dist0)
            << " packscore: " << F(7,3,selected_result.packscore)
            // << " score: "     << F(7,3,selected_result.nopackscore)
            // << " rif: "       << F(7,3,selected_result.rifscore)
            << " steric: "    << F(7,3,selected_result.stericscore)
            << " cluster: "   << I(7,selected_result.cluster_score)
            << " " << pdboutfile
            << std::endl;
            std::cout << oss.str();
            dokout << oss.str(); dokout.flush();
            
						// be careful here. Usually, the pose has been minimized, so the positions are different from
						// the BBActors....
            core::pose::Pose pose_to_dump(*selected_result.pose_);

						director->set_scene( selected_result.scene_index, iresl, *scene_minimal );
            std::vector< std::pair< int, std::string > > brians_infolabels;
            for( int i_actor = 0; i_actor < scene_minimal-> template num_actors<BBActor>(1); ++i_actor ) {
                BBActor bba = scene_minimal->template get_actor<BBActor>(1,i_actor);
                int const ires = 1+scaffres_l2g.at( bba.index_ );
                
                std::vector< std::pair<float, int > > rotscores;
                rif_ptrs.back()->get_rotamers_for_xform( bba.position(), rotscores );
								typedef std::pair<float, int> PairFI;
                BOOST_FOREACH( PairFI const & p, rotscores) {
                    int const irot = p.second;
                    std::pair< int, int > sat1_sat2 = rif_ptrs.back()->get_sat1_sat2(bba.position(), irot);
                    if (sat1_sat2.first > -1) {
                        std::pair< int, std::string > brian_pair;
                        brian_pair.first = ires;
                        brian_pair.second = "HOT_IN:" + str(sat1_sat2.first);
                        brians_infolabels.push_back(brian_pair);
                    }
                }
                
            } // loop for all BBActor positions
            
            for ( auto p : brians_infolabels ) {
                pose_to_dump.pdb_info()->add_reslabel(p.first, p.second);
            }
            
            pose_to_dump.dump_pdb(pdboutfile);
        } // end of the for loop that iterating over all the results
    } // end of the output_results function
    
    } // namespace scheme
} // namespace devel



#endif
