// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:


#ifndef INCLUDED_riflib_rifdock_subroutines_util_hh
#define INCLUDED_riflib_rifdock_subroutines_util_hh


#include <riflib/types.hh>

#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>

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
    
    } // namespace scheme
} // namespace devel



#endif
