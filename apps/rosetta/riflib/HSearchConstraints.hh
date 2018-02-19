// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols

#ifndef INCLUDED_riflib_HSeachConstraints_hh
#define INCLUDED_riflib_HSeachConstraints_hh

#include <riflib/types.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/xyzStripeHashPose.hh>
#include <utility/file/file_sys_util.hh>


namespace devel {
namespace scheme {

// constrain codes, maybe I should move it to a separate file.

class CstBase
{
public:
    virtual void set_coordinates( core::pose::Pose const & target, core::pose::Pose const & scaffold ) = 0;
    virtual void init( core::pose::Pose const & target, core::pose::Pose const & scaffold, double resl ) = 0;
    virtual bool apply(devel::scheme::EigenXform const & trans) = 0;
    virtual void reset() = 0;
    virtual shared_ptr<CstBase> clone() = 0;
    virtual ~CstBase() {};
};


typedef shared_ptr<CstBase> CstBaseOP;

class AtomPairCst: public CstBase
{
private:
    typedef numeric::xyzVector<core::Real> Vec;
    int32_t target_atom_resnum;
    std::string target_atom_name;
    int32_t scaffold_atom_resnum;
    std::string scaffold_atom_name;
    double distance;
    bool close;
    Eigen::Vector3f target_atom_coor;
    Eigen::Vector3f scaffold_atom_coor;
    double resolution;
    bool coordinates_set;
public:
    AtomPairCst(){ }
    AtomPairCst( std::string tar_name, int32_t tar_resnum, std::string scaff_name, int32_t scaff_resnum, double dist):
    target_atom_name(tar_name),
    target_atom_resnum(tar_resnum),
    scaffold_atom_resnum(scaff_resnum),
    scaffold_atom_name(scaff_name),
    distance(dist),
    close(true),
    coordinates_set(false)
    {
        if (distance < 0 ) {
            distance = -distance;
            close = false;
        }
    }
    AtomPairCst( AtomPairCst const & src ) {
        *this = src;
    }
    AtomPairCst&
    operator=( AtomPairCst const & src ) {
        if ( this != & src ) {
            target_atom_resnum = src.target_atom_resnum;
            target_atom_name = src.target_atom_name;
            scaffold_atom_resnum = src.scaffold_atom_resnum;
            scaffold_atom_name = src.scaffold_atom_name;
            distance = src.distance;
            close = src.close;
            target_atom_coor = src.target_atom_coor;
            scaffold_atom_coor = src.scaffold_atom_coor;
            resolution = src.resolution;
        }
    }
    CstBaseOP clone() {
        return make_shared<AtomPairCst>(*this);
    }
    void reset() {
        coordinates_set = false;
    } 
    void set_coordinates(core::pose::Pose const & target, core::pose::Pose const & scaffold)
    {
        if (coordinates_set) return;

        target_atom_coor[0] = target.residue(target_atom_resnum).xyz(target_atom_name).x();
        target_atom_coor[1] = target.residue(target_atom_resnum).xyz(target_atom_name).y();
        target_atom_coor[2] = target.residue(target_atom_resnum).xyz(target_atom_name).z();
        
        scaffold_atom_coor[0] = scaffold.residue(scaffold_atom_resnum).xyz(scaffold_atom_name).x();
        scaffold_atom_coor[1] = scaffold.residue(scaffold_atom_resnum).xyz(scaffold_atom_name).y();
        scaffold_atom_coor[2] = scaffold.residue(scaffold_atom_resnum).xyz(scaffold_atom_name).z();
        
        coordinates_set = true;
        return;
    }
    void init(core::pose::Pose const & target, core::pose::Pose const & scaffold, double resl) { resolution = resl; }
    bool apply(devel::scheme::EigenXform const & trans)
    {
        Eigen::Vector3f v = trans * scaffold_atom_coor;
        double d = (v - target_atom_coor).norm();
        if ( close ){
            if ( (d - resolution) < distance ) {
                return true;
            }else{
                return false;
            }
        }else{
            if ( (d + resolution) > distance ) {
                return true;
            }else{
                return false;
            }
        }
    }
};

class AtomToPoseCst: public CstBase
{
private:
    typedef numeric::xyzVector<core::Real> Vec;
    int32_t atom_resnum;
    std::string atom_name;
    double distance;
    bool close;
    Eigen::Vector3f atom_coor;
    core::pose::xyzStripeHashPoseOP hashed_pose;
    bool atom_on_target;
    bool coordinates_set;
public:
    AtomToPoseCst(){}
    AtomToPoseCst(std::string name, int32_t resnum, double dist, bool a_on_t = true /* atom on targete or scaffold*/):
    atom_resnum(resnum),
    atom_name(name),
    distance(dist),
    atom_on_target(a_on_t),
    hashed_pose(nullptr),
    close(true),
    coordinates_set(false)
    {
        if (distance < 0 ) {
            distance = -distance;
            close = false;
        }
    }
    AtomToPoseCst( AtomToPoseCst const & src ) {
        *this = src;
    }
    AtomToPoseCst&
    operator=( AtomToPoseCst const & src ) {
        if ( this != & src ) {
            atom_resnum = src.atom_resnum;
            atom_name = src.atom_name;
            distance = src.distance;
            close = src.close;
            atom_coor = src.atom_coor;
            hashed_pose = src.hashed_pose;      // this might be a problem
            atom_on_target = src.atom_on_target;
        }
    }
    CstBaseOP clone() {
        return make_shared<AtomToPoseCst>(*this);
    }
    void reset() {
        hashed_pose = nullptr;
        coordinates_set = false;
    } 
    void set_coordinates(core::pose::Pose const & target, core::pose::Pose const & scaffold)
    {
        if (coordinates_set) return;

        core::pose::Pose p;
        if (atom_on_target) {
            p = target;
        }else{
            p = scaffold;
        }
        atom_coor[0] = p.residue(atom_resnum).xyz(atom_name).x();
        atom_coor[1] = p.residue(atom_resnum).xyz(atom_name).y();
        atom_coor[2] = p.residue(atom_resnum).xyz(atom_name).z();

        coordinates_set = true;
        return;
    }
    void init(core::pose::Pose const & target, core::pose::Pose const & scaffold, double resl)
    {
        core::pose::Pose p;
        if (atom_on_target) {
            p = scaffold;
        }else{
            p = target;
        }
        if (close){
            hashed_pose = make_shared<core::pose::xyzStripeHashPose>(p, core::pose::PoseCoordPickMode_BB, distance + resl);
        }else{
            if ( (distance - resl) < 0) {
                hashed_pose = nullptr;
            }else{
                hashed_pose = make_shared<core::pose::xyzStripeHashPose>(p, core::pose::PoseCoordPickMode_BB, distance - resl);
            }
        }
        return;
    }
    bool apply(devel::scheme::EigenXform const & trans)
    {
        if (!hashed_pose) return true;
        Eigen::Vector3f v;
        if (atom_on_target) {
            v = trans.inverse() * atom_coor;
        }else{
            v = trans * atom_coor;
        }
        bool clash = hashed_pose->clash(Vec(v[0], v[1], v[2]));
        if ( clash ^ close )
            return false;
        else
            return true;
    }
    
};


class RayToPoseCst: public CstBase
{
private:
    typedef numeric::xyzVector<core::Real> Vec;
    int32_t res1_num;
    std::string res1_atom_name;
    int32_t res2_num;
    std::string res2_atom_name;
    double length;
    double distance;
    bool close;
    std::vector<Eigen::Vector3f > ray_coors;
    core::pose::xyzStripeHashPoseOP hashed_pose;
    bool atom_on_target;
    bool coordinates_set;
public:
    RayToPoseCst(){ }
    RayToPoseCst(std::string r1_name, int32_t r1_num, std::string r2_name, int32_t r2_num, double len, double dist, bool a_on_t = true /*atom on target or scaffold*/):
    res1_num(r1_num),
    res1_atom_name(r1_name),
    res2_num(r2_num),
    res2_atom_name(r2_name),
    length(len),
    distance(dist),
    close(true),
    hashed_pose(nullptr),
    atom_on_target(a_on_t),
    coordinates_set(false)
    {
        if (distance < 0) {
            distance = -distance;
            close = false;
        }
        if (distance > length) {
            distance = length;
        }
    }
    RayToPoseCst( RayToPoseCst const & src ) {
        *this = src;
    }
    RayToPoseCst&
    operator=( RayToPoseCst const & src ) {
        if ( this != & src ) {
            res1_num = src.res1_num;
            res1_atom_name = src.res1_atom_name;
            res2_num = src.res2_num;
            res2_atom_name = src.res2_atom_name;
            length = src.length;
            distance = src.distance;
            close = src.close;
            ray_coors = src.ray_coors;
            hashed_pose = src.hashed_pose;      // this might be a problem
            atom_on_target = src.atom_on_target;
        }
    }
    CstBaseOP clone() {
        return make_shared<RayToPoseCst>(*this);
    }
    void reset() {
        ray_coors.resize(0);
        hashed_pose = nullptr;
        coordinates_set = false;
    } 
    void set_coordinates(core::pose::Pose const & target, core::pose::Pose const & scaffold)
    {   
        if (coordinates_set) return;

        core::pose::Pose p;
        if (atom_on_target) {
            p = target;
        }else{
            p = scaffold;
        }
        Eigen::Vector3f atom1_coor, atom2_coor;
        atom1_coor[0] = p.residue(res1_num).xyz(res1_atom_name).x();
        atom1_coor[1] = p.residue(res1_num).xyz(res1_atom_name).y();
        atom1_coor[2] = p.residue(res1_num).xyz(res1_atom_name).z();
        
        atom2_coor[0] = p.residue(res2_num).xyz(res2_atom_name).x();
        atom2_coor[1] = p.residue(res2_num).xyz(res2_atom_name).y();
        atom2_coor[2] = p.residue(res2_num).xyz(res2_atom_name).z();
        
        Eigen::Vector3f ray_direct = (atom1_coor - atom2_coor).normalized();
        
        double two_point_dist = distance > 3.0 ? distance : 3.0;
        int32_t points_num = int((length - distance) / two_point_dist);
        for(int32_t ii = 0; ii <= points_num; ++ii)
        {
            ray_coors.push_back(atom2_coor + ii * two_point_dist * ray_direct);
        }
        coordinates_set = true;
        return;
    }
    
    void init(core::pose::Pose const & target, core::pose::Pose const & scaffold, double resl)
    {
        core::pose::PoseCoordPickMode mode;
        core::pose::Pose p;
        if (atom_on_target) {
            mode = core::pose::PoseCoordPickMode_BB;
            p = scaffold;
        }else{
            //mode = core::pose::PoseCoordPickMode_ALL; ALL or BB ????????
            mode = core::pose::PoseCoordPickMode_BB;
            p = target;
        }
        if (close) {
            hashed_pose = make_shared<core::pose::xyzStripeHashPose>(p, mode, distance + resl);
        }else{
            if(distance - resl < 0){
                hashed_pose = nullptr;
            }else{
                hashed_pose = make_shared<core::pose::xyzStripeHashPose>(p, mode, distance - resl);
            }
        }
        return;
        
    }
    
    bool apply(devel::scheme::EigenXform const & trans)
    {
        bool clash = true;
        if (!hashed_pose || ray_coors.size() == 0) {
            return true;
        }
        Eigen::Vector3f v;
        if (close) {
            for (int32_t ii = 0; ii < ray_coors.size(); ++ii)
            {
                if (atom_on_target) {
                    v = trans.inverse() * ray_coors[ii];
                }else{
                    v = trans * ray_coors[ii];
                }
                clash = hashed_pose->clash(Vec(v[0], v[1], v[2]));
                if (clash) {
                    return true;
                }
            }
            return false;
        }else{
            for (int32_t ii=0; ii < ray_coors.size(); ++ii) {
                if (atom_on_target) {
                    v = trans.inverse() * ray_coors[ii];
                }else{
                    v = trans * ray_coors[ii];
                }
                clash = hashed_pose->clash(Vec(v[0], v[1], v[2]));
                if (clash) {
                    return false;
                }
            }
            return true;
        }
    }
};

typedef shared_ptr<AtomPairCst> AtomPairCstOP;
typedef shared_ptr<AtomToPoseCst> AtomToPoseCstOP;
typedef shared_ptr<RayToPoseCst> RayToPoseCstOP;

inline
bool parse_constrains_file(std::string fname, std::vector<CstBaseOP> & csts)
{
    if (fname == "" ) {
        return false;
    }
    runtime_assert_msg(utility::file::file_exists( fname ), "constrain file does not exits: " + fname );
        std::ifstream in;
        in.open(fname, std::ios::in);
    // utility::io::izstream in(fname);
    std::string s;
    while (std::getline(in, s)) {
                if (s.empty()) continue;
        utility::vector1<std::string> splt = utility::string_split_simple(s, ' ');
        if (splt[1].find("#") == 0  || splt.size() == 0) {
            continue;
        }
                std::cout << "Parsing constrians: " + s << std::endl; 
        if ( "AtomPair" == splt[1]) {
            csts.push_back(make_shared<AtomPairCst>(splt[2], utility::string2int(splt[3]), splt[4], utility::string2int(splt[5]), utility::string2float(splt[6]) ));
        }else if("AtomToScaffold" == splt[1]){
            csts.push_back(make_shared<AtomToPoseCst>(splt[2], utility::string2int(splt[3]), utility::string2float(splt[4]), true ));
        }else if("AtomToTarget" == splt[1]){
            csts.push_back(make_shared<AtomToPoseCst>(splt[2], utility::string2int(splt[3]), utility::string2float(splt[4]), false ));
        }else if("RayToScaffold" == splt[1]){
            csts.push_back(make_shared<RayToPoseCst>(splt[2], utility::string2int(splt[3]), splt[4], utility::string2int(splt[5]), utility::string2float(splt[6]), utility::string2float(splt[7]), true ));
        }else if("RayToTarget" == splt[1]){
            csts.push_back(make_shared<RayToPoseCst>(splt[2], utility::string2int(splt[3]), splt[4], utility::string2int(splt[5]), utility::string2float(splt[6]), utility::string2float(splt[7]), false ));
        }else{
            return false;
        }
    }
    return true;
}


}}

#endif