// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols

#ifndef INCLUDED_riflib_scaffold_NineAManager_hh
#define INCLUDED_riflib_scaffold_NineAManager_hh

#include <riflib/types.hh>
#include <riflib/rifdock_typedefs.hh>
#include <scheme/scaffold/ScaffoldProviderBase.hh>
#include <riflib/scaffold/ScaffoldDataCache.hh>
#include <riflib/scaffold/nineA_util.hh>
#include <riflib/HSearchConstraints.hh>
#include <riflib/scaffold/MorphingScaffoldProvider.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyz.functions.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/PDBInfo.hh>

#include <ObjexxFCL/format.hh>

#include <string>
#include <vector>
#include <boost/any.hpp>
#include <boost/format.hpp>
#include <unordered_map>





namespace devel {
namespace scheme {



struct NineAMember {
    ParametricSceneConformationCOP conformation;
    uint64_t cdindex;
    uint64_t clust;
    ::scheme::scaffold::TreeRelation tree_relation;
};


namespace DanielAtomIndex {
    const uint64_t N = 3; 
    const uint64_t CA = 0; 
    const uint64_t C = 1; 
    const uint64_t O = 2; 
}

enum IntFields {
    Clust = 0,
    NumCon,
    PrevAssig,
    last_int_field
};



enum FloatFields {
    AvgRad = 0,
    MaxRad,
    X_1,
    X_2,
    X_3,
    X_4,
    X_5,
    X_6,
    X_7,
    X_8,
    X_9,
    X_10,
    X_11,
    X_12,
    X_13,
    X_14,
    X_15,
    X_16,
    X_17,
    X_18,
    X_19,
    X_20,
    X_21,
    X_22,
    X_23,
    X_24,
    X_25,
    X_26,
    X_27,
    X_28,
    X_29,
    X_30,
    X_31,
    X_32,
    X_33,
    X_34,
    X_35,
    X_36,
    Y_1,
    Y_2,
    Y_3,
    Y_4,
    Y_5,
    Y_6,
    Y_7,
    Y_8,
    Y_9,
    Y_10,
    Y_11,
    Y_12,
    Y_13,
    Y_14,
    Y_15,
    Y_16,
    Y_17,
    Y_18,
    Y_19,
    Y_20,
    Y_21,
    Y_22,
    Y_23,
    Y_24,
    Y_25,
    Y_26,
    Y_27,
    Y_28,
    Y_29,
    Y_30,
    Y_31,
    Y_32,
    Y_33,
    Y_34,
    Y_35,
    Y_36,
    Z_1,
    Z_2,
    Z_3,
    Z_4,
    Z_5,
    Z_6,
    Z_7,
    Z_8,
    Z_9,
    Z_10,
    Z_11,
    Z_12,
    Z_13,
    Z_14,
    Z_15,
    Z_16,
    Z_17,
    Z_18,
    Z_19,
    Z_20,
    Z_21,
    Z_22,
    Z_23,
    Z_24,
    Z_25,
    Z_26,
    Z_27,
    Z_28,
    Z_29,
    Z_30,
    Z_31,
    Z_32,
    Z_33,
    Z_34,
    Z_35,
    Z_36,
    Phi_1,
    Psi_1,
    Ome_1,
    Phi_2,
    Psi_2,
    Ome_2,
    Phi_3,
    Psi_3,
    Ome_3,
    Phi_4,
    Psi_4,
    Ome_4,
    Phi_5,
    Psi_5,
    Ome_5,
    Phi_6,
    Psi_6,
    Ome_6,
    Phi_7,
    Psi_7,
    Ome_7,
    Phi_8,
    Psi_8,
    Ome_8,
    Phi_9,
    Psi_9,
    Ome_9,
    last_float_field
};



struct StatRow {
    std::vector<uint32_t> int_fields;
    std::vector<float> float_fields;
    std::string filename;

    StatRow() {
        int_fields.resize(last_int_field);
        float_fields.resize(last_float_field);
    }
};

typedef utility::vector1<StatRow> StatTable;
typedef std::shared_ptr<StatTable> StatTableOP;
typedef std::shared_ptr<StatTable const> StatTableCOP;

typedef utility::vector1<std::vector<uint64_t>> ChildrenByClust;
typedef std::shared_ptr<ChildrenByClust> ChildrenByClustOP;
typedef std::shared_ptr<ChildrenByClust const> ChildrenByClustCOP;

const std::vector<std::string> CLUSTER_DATA_NAMES {
    "_res5.22",
    "_res5.26",
    "_res5.31",
    "_res5.24",
    "_res5.32",
    "_res2.23",
    "_res1.90",
    "_res1.86",
    "_res1.65",
    "_res1.58",
    "final_res1.41"
};

const uint64_t NUM_CLUSTERS = CLUSTER_DATA_NAMES.size();


struct NineAManager : public utility::pointer::ReferenceCount {

    static
    shared_ptr<NineAManager>
    get_instance(
        shared_ptr< RotamerIndex > rot_index_p_in,
        RifDockOpt const & opt_in ) {

        if ( ! single_instance_ ) {
            single_instance_ = make_shared<NineAManager>( rot_index_p_in, opt_in );
        }

        return single_instance_;
    }

// Don't call this it's supposed to be private
    NineAManager(
        shared_ptr< RotamerIndex > rot_index_p_in,
        RifDockOpt const & opt_in ) :

        rot_index_p( rot_index_p_in ),
        opt( opt_in ) {
            tables_.resize( NUM_CLUSTERS );
        }



    NineAMember
    get_nineA_member( uint64_t cdindex, uint64_t clust ) {
        NineAMember nine;
        nine.cdindex = cdindex;
        nine.clust = clust;

        core::pose::PoseOP pose_p = make_fragment( cdindex, clust );
        std::string scafftag = pose_p->pdb_info()->name();
        utility::vector1<core::Size> scaffold_res;
        EigenXform scaffold_perturb = EigenXform::Identity();
        std::vector<CstBaseOP> csts;

        get_default_scaffold_res(*pose_p, scaffold_res);

        ScaffoldDataCacheOP temp_data_cache = make_shared<ScaffoldDataCache>(
            *pose_p,
            scaffold_res,
            scafftag,
            scaffold_perturb,
            rot_index_p,
            opt,
            csts);

        nine.conformation = make_conformation_from_data_cache(temp_data_cache, false);

        return nine;
    }


    uint64_t
    size( uint64_t cdindex ) {
        return load_stat_table( cdindex )->size();
    }



private:
    core::pose::PoseOP
    make_fragment( uint64_t cdindex, uint64_t clust ) {
        using numeric::xyzVector;
        using core::Real;

        StatTableCOP stats = load_stat_table( cdindex );
        StatRow const & sr = stats->at( clust );

        core::pose::PoseOP pose_p = make_shared<core::pose::Pose>();
        core::pose::Pose & pose = *pose_p;
        core::pose::make_pose_from_sequence(pose, "AAAAAAAAA", "fa_standard", false);
        add_pdbinfo_if_missing( pose );

        // Set the phis and psis
        for ( uint64_t i = 0; i < 9; i++ ) {
            pose.set_phi( i + 1, sr.float_fields[static_cast<int>(Phi_1) + i*3]);
            pose.set_psi( i + 1, sr.float_fields[static_cast<int>(Psi_1) + i*3]);
            pose.set_omega( i + 1, sr.float_fields[static_cast<int>(Ome_1) + i*3]);
        }

        // Set the name
        std::string name = boost::str( boost::format( "da9_%02i_%07i" ) % cdindex % clust );
        pose.pdb_info()->name( name );


        // align the pose to the xyz data
        std::vector<xyzVector<Real>> atom_vec = get_fragment_xyzVectors( sr );
        xyzVector<Real> my_first_vec = pose.residue(1).xyz("CA") - pose.residue(1).xyz("N");
        xyzVector<Real> my_last_vec = pose.residue(9).xyz("CA") - pose.residue(9).xyz("N");

        xyzVector<Real> real_first_vec = atom_vec[(1-1)*4+DanielAtomIndex::CA] - atom_vec[(1-1)*4+DanielAtomIndex::N];
        xyzVector<Real> real_last_vec =  atom_vec[(9-1)*4+DanielAtomIndex::CA] - atom_vec[(9-1)*4+DanielAtomIndex::N];

        numeric::xyzMatrix<Real> rot = numeric::alignVectorSets(my_first_vec, my_last_vec, real_first_vec, real_last_vec);
        xyzVector<Real> z( 0, 0, 0 );
        pose.apply_transform_Rx_plus_v(rot, z);

        xyzVector<Real> trans = atom_vec[(1-1)*4+DanielAtomIndex::CA] - pose.residue(1).xyz("CA");
        numeric::xyzMatrix<Real> no_rot = numeric::z_rotation_matrix(0);

        pose.apply_transform_Rx_plus_v(no_rot, trans);

        Real msd = 0;
        for ( uint64_t i = 1; i <= 9; i ++ ) {
            msd += atom_vec[(i-1)*4+DanielAtomIndex::CA].distance_squared(pose.residue(i).xyz("CA"));
        }

        if ( msd > 0.3f*0.3f ) {
            std::cerr << "Warning: fragment " << name << " rmsd " << std::sqrt(msd) << std::endl;
        }

        return pose_p;

    }





    std::vector<numeric::xyzVector<core::Real>>
    get_fragment_xyzVectors( StatRow const & sr ) {
        std::vector<numeric::xyzVector<core::Real>> to_ret;
        for ( uint64_t i = 0; i < 36; i ++ ) {
            to_ret.push_back( numeric::xyzVector<core::Real>( 
                sr.float_fields[static_cast<int>(X_1) + i],
                sr.float_fields[static_cast<int>(Y_1) + i],
                sr.float_fields[static_cast<int>(Z_1) + i]
            ));
        }
        return to_ret;
    }


    ChildrenByClustCOP
    load_children_by_clust( uint64_t parent_cdindex, uint64_t child_cdindex ) {
        using ObjexxFCL::format::I;

        std::string key = I(3, parent_cdindex) + "_" + I(3, child_cdindex );

        if ( children_by_clusts_.count( key ) == 0 ) {

            StatTableCOP stats = load_stat_table( child_cdindex );

            ChildrenByClustOP by_clust = make_shared<ChildrenByClust>();


            for ( uint64_t clust_in_child = 1; clust_in_child <= stats->size(); clust_in_child++ ) {
                StatRow const & sr = (*stats)[clust_in_child];
                uint64_t clust_in_parent = sr.int_fields[PrevAssig];

                if (clust_in_parent >= by_clust->size()) {
                    by_clust->resize(clust_in_parent - by_clust->size() + 1);
                }
                by_clust->at( clust_in_parent ).push_back( clust_in_child );
            }

            children_by_clusts_[key] = by_clust;

        }

        return children_by_clusts_.at(key);
    }

    std::vector<uint64_t>
    get_fragment_children_clusts( uint64_t parent_clust, uint64_t parent_cdindex, uint64_t child_cdindex ) {
        ChildrenByClustCOP by_clust = load_children_by_clust( parent_cdindex, child_cdindex );
        return by_clust->at(parent_clust);
    }


    void
    read_int( StatRow & sr, uint64_t & pos, utility::vector1< std::string > & sp, IntFields intf ) {
        sr.int_fields[intf] = stoi(sp[pos++]); 
    }

    void
    read_float( StatRow & sr, uint64_t & pos, utility::vector1< std::string > & sp, FloatFields floatf ) {
        sr.float_fields[floatf] = stof(sp[pos++]); 
    }


    StatTableCOP
    load_stat_table( uint64_t cdindex ) {
        using ObjexxFCL::format::I;

        if ( ! tables_.at(cdindex) ) {
            


            std::string filename = opt.nineA_cluster_path 
                                    + "/kcenters_stats_al1.dat" 
                                    + CLUSTER_DATA_NAMES.at( cdindex );

            std::ifstream f( filename );

            if ( !f ) {
                utility_exit_with_message("nineA_cluster_path file not found: " + filename);
            }

            std::string line;
            std::getline( f, line );

            utility::vector1< std::string > string_split = utility::string_split( line, '\t' );

            if ( string_split.size() != 395 ) {
                utility_exit_with_message("nineA_cluster_path file has wrong number of columns:  " 
                    + I(3,string_split.size()) + " != 395 " + filename);
            }


            StatTableOP table = make_shared<StatTable>();

            while ( std::getline( f, line ) ) {
                utility::vector1< std::string > sp = utility::string_split( line, '\t' );
                if ( sp.size() <= 1 ) {
                    break;
                }
                if ( sp.size() != 395 ) {
                    utility_exit_with_message("nineA_cluster_path file has wrong number of columns inside:  " 
                    + I(3,sp.size()) + " != 395 " + filename);
                }

                StatRow sr;
                uint64_t pos = 1;

                pos++;
                read_int(sr,pos,sp, Clust);
                read_int(sr,pos,sp, NumCon);
                read_float(sr,pos,sp, AvgRad);
                read_float(sr,pos,sp, MaxRad);
                for ( int i = 0; i < 36*3; i++ ) {
                    read_float(sr,pos,sp, static_cast<FloatFields>(static_cast<int>(X_1) + i));
                }
                for ( int i = 0; i < 3*9; i++ ) {
                    read_float(sr,pos,sp, static_cast<FloatFields>(static_cast<int>(Phi_1) + i));
                }
                pos += 393 - 141;
                sr.filename = sp[pos++];
                pos++;
                read_int(sr,pos,sp, PrevAssig);

                table->push_back( sr );
            }

            f.close();

            tables_[cdindex] = table;
        }

        return tables_[cdindex];

    }

    static shared_ptr<NineAManager> single_instance_;

    std::vector< StatTableOP > tables_;
    std::unordered_map<std::string, ChildrenByClustOP> children_by_clusts_;

    shared_ptr< RotamerIndex > rot_index_p;
    RifDockOpt const & opt;

};



}}

#endif
