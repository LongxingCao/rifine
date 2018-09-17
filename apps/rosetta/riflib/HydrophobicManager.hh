// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols

#ifndef INCLUDED_riflib_HydrophobicManager_hh
#define INCLUDED_riflib_HydrophobicManager_hh

#include <riflib/types.hh>
//#include <riflib/rifdock_typedefs.hh>

#include "scheme/util/SimpleArray.hh"
#include <riflib/RotamerGenerator.hh>

#include <core/pose/Pose.hh>

#include <string>
#include <vector>

// on my laptop, I must add this to use the unordered_map, how about on the digs?
#include<tr1/unordered_map>

#include <boost/any.hpp>
#include <boost/format.hpp>

#include <riflib/util.hh>

#include <ObjexxFCL/format.hh>


namespace devel {
namespace scheme {


// This class is pretty sketch
// We're trying to figure out whether the scaffold interacts with each target hydrophobic residue with more than
//   0.5 kcal fa_sol, fa_atr, fa_rep
//
// The plan:
// Make a 0.5 A grid
// Give each hydrophobic a Hyd number
// Start filling in the grid for each hydrophobic carbon
// From 3.0 - 4.5: fill the grid twice with the Hyd // not doing this for now
// From 4.5 - 6.0: fill the grid once with the Hyd

// Then after packing:
// Scan through all rifres carbon atoms and start tallying the Hyds
// Any Hyd that gets more than 2 is tallied
// This count is the hydrophobic_residue_contacts



struct HydrophobicManager {

    static constexpr float DDG_PER_COUNT = -0.28f;

    typedef uint16_t Hyd;
    static int const CACHE_MAX_HYD = std::numeric_limits<Hyd>::max();

    typedef ::scheme::util::SimpleArray<3,size_t> Indices;
    typedef Eigen::Vector3f Bounds;
    Bounds lb_,ub_,cs_;

    const float max_interaction_range_ = 7.0;   // It's actually 6.0, but just to be safe
    size_t max_hyds_;

    const std::set<char> hydrophobic_name1s_ {'A', 'C', 'F', 'I', 'L', 'M', 'P', 'T', 'V', 'W', 'Y'}; 
    std::vector<core::Size> hydrophobic_res_;

    // Cache these for speed
    std::vector<int> rif_atype_map_;
    std::vector<bool> rif_hydrophobic_map_;
    /*
    std::vector<bool> rif_pi_map_;
    */


    std::vector<Hyd> voxel_map_;

    // This only actually needs 1. But include a few extra in case it accidentally gets written to
    std::vector<Hyd> OOB_LIST = {CACHE_MAX_HYD, CACHE_MAX_HYD, CACHE_MAX_HYD, CACHE_MAX_HYD, CACHE_MAX_HYD, CACHE_MAX_HYD};

    Indices shape_;

    /*
    Bounds cat_lb_,cat_ub_,cat_cs_;
    std::vector<core::Size> cation_res_;
    std::vector<Hyd> cation_voxel_map_;
    Indices cation_shape_;
    size_t cation_max_hyds_;
    */


    shared_ptr< RotamerIndex > rot_index_p;

    float hydrophobic_ddg_cut_ = 99999;
    int hydrophobic_residue_contacts_cut_ = -99999;

    float one_hydrophobic_better_than_ = 0;
    float two_hydrophobics_better_than_ = 0;
    float three_hydrophobics_better_than_ = 0;
    float hydrophobic_ddg_per_atom_cut_ = 0;
    bool doing_better_than_ = false;

    /*
    int num_cation_pi_ = 0;
    */


    HydrophobicManager(
        core::pose::Pose const & target,
        utility::vector1<core::Size> const & target_res,
        shared_ptr< RotamerIndex > rot_index_p_in
    ) : rot_index_p( rot_index_p_in )
    {

        rif_atype_map_ = get_rif_atype_map();

        // move the get_rif_hydrophobic_map here, actually this is corresponding to the util.cc file
        rif_hydrophobic_map_.resize( rif_atype_map_.size(), false );
        rif_hydrophobic_map_[ 1] = false; //"CNH2"
        rif_hydrophobic_map_[ 2] = false; //"COO" 
        rif_hydrophobic_map_[ 3] =  true; //"CH1" 
        rif_hydrophobic_map_[ 4] =  true; //"CH2" 
        rif_hydrophobic_map_[ 5] =  true; //"CH3" 
        rif_hydrophobic_map_[ 6] =  true; //"aroC"
        rif_hydrophobic_map_[ 7] = false; //"Ntrp"
        rif_hydrophobic_map_[ 8] = false; //"Nhis"
        rif_hydrophobic_map_[ 9] = false; //"NH2O"
        rif_hydrophobic_map_[10] = false; //"Nlys"
        rif_hydrophobic_map_[11] = false; //"Narg"
        rif_hydrophobic_map_[12] = false; //"Npro"
        rif_hydrophobic_map_[13] = false; //"OH"  
        rif_hydrophobic_map_[14] = false; //"ONH2"
        rif_hydrophobic_map_[15] = false; //"OOC" 
        rif_hydrophobic_map_[16] = false; //"Oaro"
        rif_hydrophobic_map_[17] =  true; //"S"   
        rif_hydrophobic_map_[18] = false; //"Nbb" 
        rif_hydrophobic_map_[19] =  true; //"CAbb"
        rif_hydrophobic_map_[20] = false; //"CObb"
        rif_hydrophobic_map_[21] = false; //"OCbb"

        /*
        rif_pi_map_.resize( rif_atype_map_.size(), false );
        hydrophobic_map[ 1] = false; //"CNH2"
        hydrophobic_map[ 2] =  true; //"COO" 
        hydrophobic_map[ 3] = false; //"CH1" 
        hydrophobic_map[ 4] = false; //"CH2" 
        hydrophobic_map[ 5] = false; //"CH3" 
        hydrophobic_map[ 6] =  true; //"aroC"
        hydrophobic_map[ 7] =  true; //"Ntrp"
        hydrophobic_map[ 8] =  true; //"Nhis"
        hydrophobic_map[ 9] = false; //"NH2O"
        hydrophobic_map[10] = false; //"Nlys"
        hydrophobic_map[11] = false; //"Narg"
        hydrophobic_map[12] = false; //"Npro"
        hydrophobic_map[13] = false; //"OH"  
        hydrophobic_map[14] = false; //"ONH2"
        hydrophobic_map[15] = false; //"OOC" 
        hydrophobic_map[16] = false; //"Oaro"
        hydrophobic_map[17] = false; //"S"   
        hydrophobic_map[18] = false; //"Nbb" 
        hydrophobic_map[19] = false; //"CAbb"
        hydrophobic_map[20] = false; //"CObb"
        hydrophobic_map[21] = false; //"OCbb"
        */

        identify_hydrophobic_residues( target, target_res );
        prepare_bounds( target );
        std::vector<std::vector<Hyd>> early_map = first_pass_fill( target );
        create_and_fill_voxel_map( early_map );

        /*
        identify_cation_residues( target, target_res );
        cation_prepare_bounds( target );
        std::vector<std::set<Hyd>> cation_early_map = cation_first_pass_fill( target );
        cation_create_and_fill_voxel_map( cation_early_map );
        */



    }

    void
    setup_criteria( float hydrophobic_ddg, int hydrophobic_residue_contacts, float one, float two, float three, float ddg_per_atom_cut ) {
        hydrophobic_ddg_cut_ = hydrophobic_ddg;
        hydrophobic_residue_contacts_cut_ = hydrophobic_residue_contacts;
        one_hydrophobic_better_than_ = one;
        two_hydrophobics_better_than_ = two;
        three_hydrophobics_better_than_ = three;
        hydrophobic_ddg_per_atom_cut_ = ddg_per_atom_cut;

        doing_better_than_ = ( one < 0 || two < 0 || three < 0 );
    }

    /*
    void
    set_num_cation_pi( int num_pi ) {
        num_cation_pi_ = num_pi;
    }
    */


    bool
    is_rosetta_atom_hydrophobic( core::conformation::Residue const & res, core::Size atno ) {
        int iatype = rif_atype_map_.at(res.atom_type_index(atno));
        if (iatype >= rif_hydrophobic_map_.size()) {
            return false;
        }
        return rif_hydrophobic_map_.at(iatype);
    }

    void
    identify_hydrophobic_residues( 
        core::pose::Pose const & target,
        utility::vector1<core::Size> const & target_res
    ) {
        hydrophobic_res_.clear();

        for ( core::Size seqpos : target_res ) {
            if ( hydrophobic_name1s_.count( target.residue(seqpos).name1() ) != 0 ) {
                hydrophobic_res_.push_back(seqpos);
            }
        }
    }

    void
    prepare_bounds( core::pose::Pose const & target ) {

        Eigen::Vector3f lbs( 9e9, 9e9, 9e9 );
        Eigen::Vector3f ubs( -9e9, -9e9, -9e9 );

        for ( core::Size seqpos : hydrophobic_res_ ) {
            core::conformation::Residue const & res = target.residue(seqpos);
            for ( core::Size atno = 4; atno < res.nheavyatoms(); atno++ ) {  // Don't include backbone
                if ( ! is_rosetta_atom_hydrophobic( res, atno ) ) continue;

                numeric::xyzVector<core::Real> xyz = res.xyz(atno);
                for ( int i = 0; i < 3; i++ ) {
                    lbs[i] = std::min<float>( lbs[i], xyz[i] );
                    ubs[i] = std::max<float>( ubs[i], xyz[i] );
                }
            }
            
        }

        for ( int i = 0; i < 3; i++ ) {
            lbs[i] -= max_interaction_range_;
            ubs[i] += max_interaction_range_;
        }

        lb_ = lbs;
        ub_ = ubs;
        cs_ = Eigen::Vector3f( 0.5, 0.5, 0.5 );

        Indices extents = floats_to_index( ub_ );
        shape_ = extents + Indices(1);

    }

    std::vector<std::vector<Hyd>>
    first_pass_fill( core::pose::Pose const & target ) {

        std::vector<std::vector<Hyd>> early_map;

        size_t elements = shape_[0] * shape_[1] * shape_[2];
        // std::cout << shape_ << std::endl;
        // std::cout << floats_to_index( ub_ ) << std::endl;
        // std::cout << elements << " " << index_to_map_index( floats_to_index( ub_ ) ) << std::endl;
        // this is because the indexing is from 0.
        runtime_assert( elements - 1 == index_to_map_index( floats_to_index( ub_ ) ) );

        early_map.resize(elements);

        for ( Hyd ihyd = 0; ihyd < hydrophobic_res_.size(); ihyd++ ) {
            core::conformation::Residue const & res = target.residue(hydrophobic_res_[ihyd]);

            for ( core::Size atno = 4; atno < res.nheavyatoms(); atno++ ) {  // Don't include backbone
                if ( ! is_rosetta_atom_hydrophobic( res, atno ) ) continue;

                numeric::xyzVector<core::Real> _xyz = res.xyz( atno );
                Eigen::Vector3f xyz; xyz[0] = _xyz[0]; xyz[1] = _xyz[1]; xyz[2] = _xyz[2];

                Eigen::Vector3f lbs( xyz[0] - max_interaction_range_, xyz[1] - max_interaction_range_, xyz[2] - max_interaction_range_ );
                Eigen::Vector3f ubs( xyz[0] + max_interaction_range_, xyz[1] + max_interaction_range_, xyz[2] + max_interaction_range_ );

                const float low_rad_sq = 3.0f*3.0f;
                const float med_rad_sq = 4.5f*4.5f;
                const float long_rad_sq = 4.8f*4.8f;

                const float step = cs_[0];

                Eigen::Vector3f worker;

                for ( float x = lbs[0] - step/2; x < ubs[0] + step; x += step ) {
                    if ( x < lb_[0] || x > ub_[0] ) continue;
                    worker[0] = x;

                    for ( float y = lbs[1] - step/2; y < ubs[1] + step; y += step ) {
                        if ( y < lb_[1] || y > ub_[1] ) continue;
                        worker[1] = y;

                        for ( float z = lbs[2] - step/2; z < ubs[2] + step; z += step ) {
                            if ( z < lb_[2] || z > ub_[2] ) continue;
                            worker[2] = z;

                            const float squared_dist = ( xyz - worker ).squaredNorm();

                            if ( squared_dist < low_rad_sq ) {
                                // do nothing
                            } else if ( squared_dist < med_rad_sq ) {
                                size_t offset = index_to_map_index( floats_to_index( worker ) );
                                early_map.at(offset).push_back( ihyd );
                                // here only fill it once
                                // early_map.at(offset).push_back( ihyd );
                            } else if ( squared_dist < long_rad_sq ) {
                                size_t offset = index_to_map_index( floats_to_index( worker ) );
                                early_map.at(offset).push_back( ihyd );
                            } else {
                                // do nothing
                            }
                        }
                    }
                }
            }
        }
        return early_map;
    }


    void
    create_and_fill_voxel_map( std::vector<std::vector<Hyd>> const & early_map ) {

        max_hyds_ = 0;

        for ( std::vector<Hyd> these_hyds : early_map ) {
            max_hyds_ = std::max<size_t>( max_hyds_, these_hyds.size() );
        }

        std::cout << "Max hydrophobics at one voxel: " << max_hyds_ << std::endl;

        // max_hyds_ is 1 too big. This means that the receiver can always loop until they hit a CACHE_MAX_HYD
        max_hyds_ += 1;

        voxel_map_.resize( max_hyds_ * shape_[0] * shape_[1] * shape_[2], CACHE_MAX_HYD );

        // We already asserted the two maps have the same size, so just array fill
        for ( size_t imap = 0; imap < early_map.size(); imap++ ) {

            std::vector<Hyd> these_hyds = early_map[imap];
            bool wrote_a_null = false;
            for ( size_t idx = 0; idx < max_hyds_; idx++ ) {
                size_t offset = imap * max_hyds_ + idx;
                if ( idx < these_hyds.size() ) {
                    voxel_map_.at(offset) = these_hyds.at(idx);
                } else {
                    voxel_map_.at(offset) = CACHE_MAX_HYD;
                    wrote_a_null = true;
                }
            }
            runtime_assert( wrote_a_null );     // Redundancy. Make sure each list is terminated
        }
    }

///////////////////////////// CATION PI ////////////////////////////
    /*
    void
    identify_cation_residues( 
        core::pose::Pose const & target,
        utility::vector1<core::Size> const & target_res
    ) {
        cation_res_.clear();

        for ( core::Size seqpos : target_res ) {
            if ( target.residue(seqpos).name1() == 'R'  ) {
                cation_res_.push_back(seqpos);
            }
        }
    }
    

    void
    cation_prepare_bounds( core::pose::Pose const & target ) {

        Eigen::Vector3f lbs( 9e9, 9e9, 9e9 );
        Eigen::Vector3f ubs( -9e9, -9e9, -9e9 );

        for ( core::Size seqpos : cation_res_ ) {
            core::conformation::Residue const & res = target.residue(seqpos);
            for ( core::Size atno = 9; atno <= 12; atno++ ) {  // The arg head

                numeric::xyzVector<core::Real> xyz = res.xyz(atno);
                for ( int i = 0; i < 3; i++ ) {
                    lbs[i] = std::min<float>( lbs[i], xyz[i] );
                    ubs[i] = std::max<float>( ubs[i], xyz[i] );
                }
            }
            
        }

        for ( int i = 0; i < 3; i++ ) {
            lbs[i] -= max_interaction_range_;
            ubs[i] += max_interaction_range_;
        }

        cat_lb_ = lbs;
        cat_ub_ = ubs;
        cat_cs_ = Eigen::Vector3f( 0.5, 0.5, 0.5 );

        Indices extents = cation_floats_to_index( cat_ub_ );
        cation_shape_ = extents + Indices(1);

    }

    std::vector<std::set<Hyd>>
    cation_first_pass_fill( core::pose::Pose const & target ) {

        std::vector<std::set<Hyd>> early_map;

        size_t elements = cation_shape_[0] * cation_shape_[1] * cation_shape_[2];
        // std::cout << shape_ << std::endl;
        // std::cout << floats_to_index( ub_ ) << std::endl;
        // std::cout << elements << " " << index_to_map_index( floats_to_index( ub_ ) ) << std::endl;
        runtime_assert( elements - 1 == cation_index_to_map_index( cation_floats_to_index( cat_ub_ ) ) );

        early_map.resize(elements);

        for ( Hyd ihyd = 0; ihyd < cation_res_.size(); ihyd++ ) {
            core::conformation::Residue const & res = target.residue(cation_res_[ihyd]);

            numeric::xyzVector<core::Real> _xyz;

            _xyz = res.xyz("NE");
            Eigen::Vector3f ne; ne[0] = _xyz[0]; ne[1] = _xyz[1]; ne[2] = _xyz[2];

            _xyz = res.xyz("CZ");
            Eigen::Vector3f cz; cz[0] = _xyz[0]; cz[1] = _xyz[1]; cz[2] = _xyz[2];

            _xyz = res.xyz("NH1");
            Eigen::Vector3f nh1; nh1[0] = _xyz[0]; nh1[1] = _xyz[1]; nh1[2] = _xyz[2];

            _xyz = res.xyz("NH2");
            Eigen::Vector3f nh2; nh2[0] = _xyz[0]; nh2[1] = _xyz[1]; nh2[2] = _xyz[2];

            Eigen::Vector3f normal = ( nh2 - cz ).cross( nh1 - cz );
            normal /= normal.norm();

            //https://stackoverflow.com/questions/47932955/how-to-check-if-a-3d-point-is-inside-a-cylinder

            // Define a cylinder with radius 2.2 A from the center of ARG
            // then on either side of the cylinder, you have to be between 2.5 and 4.5 A away
            // Actually do it as 2 cylinders

            float radius = 2.2;

            Eigen::Vector3f cy1_p1 = cz + normal * 2.5f;
            Eigen::Vector3f cy1_p2 = cz + normal * 4.5f;
            Eigen::Vector3f cy2_p1 = cz + normal * -2.5f;
            Eigen::Vector3f cy2_p2 = cz + normal * -4.5f;

            Eigen::Vector3f cy1_p1_min_p2 = cy1_p1 - cy1_p2;
            Eigen::Vector3f cy2_p1_min_p2 = cy2_p1 - cy2_p2;

            float cy1_p1_min_p2_norm = cy1_p1_min_p2.norm();
            float cy2_p1_min_p2_norm = cy2_p1_min_p2.norm();


            Eigen::Vector3f lbs( cz[0] - max_interaction_range_, cz[1] - max_interaction_range_, cz[2] - max_interaction_range_ );
            Eigen::Vector3f ubs( cz[0] + max_interaction_range_, cz[1] + max_interaction_range_, cz[2] + max_interaction_range_ );

            const float step = cs_[0];

            Eigen::Vector3f worker;
                                                                    // way over-sample
            for ( float x = lbs[0] - step/2; x < ubs[0] + step; x += step/4.0 ) {
                if ( x < lb_[0] || x > ub_[0] ) continue;
                worker[0] = x;

                for ( float y = lbs[1] - step/2; y < ubs[1] + step; y += step/4.0 ) {
                    if ( y < lb_[1] || y > ub_[1] ) continue;
                    worker[1] = y;

                    for ( float z = lbs[2] - step/2; z < ubs[2] + step; z += step/4.0 ) {
                        if ( z < lb_[2] || z > ub_[2] ) continue;
                        worker[2] = z;

                        // Cylinder 1

                        // Test between planes
                        if ( (worker - cy1_p1).dot( -cy1_p1_min_p2) >= 0 ) {
                        if ( (worker - cy1_p2).dot(  cy1_p1_min_p2) >= 0 ) {
                        if ( 
                            ( ( worker - cy1_p1).cross( -cy1_p1_min_p2 ) ).norm() 
                                                    /
                                            cy1_p1_min_p2_norm                       <= radius
                                                    ) {

                            size_t offset = cation_index_to_map_index( cation_floats_to_index( worker ) );
                            early_map.at(offset).insert( ihyd );

                        }
                        }
                        }
                        // Cylinder 2

                        // Test between planes
                        if ( (worker - cy2_p1).dot( -cy2_p1_min_p2) >= 0 ) {
                        if ( (worker - cy2_p2).dot(  cy2_p1_min_p2) >= 0 ) {
                        // Test inside curved space
                        if ( 
                            ( ( worker - cy2_p1).cross( -cy2_p1_min_p2 ) ).norm() 
                                                    /
                                            cy2_p1_min_p2_norm                       <= radius
                                                    ) {

                            size_t offset = cation_index_to_map_index( cation_floats_to_index( worker ) );
                            early_map.at(offset).insert( ihyd );

                        }
                        }
                        }

                    }
                }
            }
            
        }
        return early_map;
    }



    void
    cation_create_and_fill_voxel_map( std::vector<std::set<Hyd>> const & early_map ) {

        cation_max_hyds_ = 0;

        for ( std::set<Hyd> these_hyds : early_map ) {
            cation_max_hyds_ = std::max<size_t>( cation_max_hyds_, these_hyds.size() );
        }

        std::cout << "Max cations at one voxel: " << cation_max_hyds_ << std::endl;

        // max_hyds_ is 1 too big. This means that the receiver can always loop until they hit a CACHE_MAX_HYD
        cation_max_hyds_ += 1;

        cation_voxel_map_.resize( cation_max_hyds_ * cation_shape_[0] * cation_shape_[1] * cation_shape_[2], CACHE_MAX_HYD );

        // We already asserted the two maps have the same size, so just array fill
        for ( size_t imap = 0; imap < early_map.size(); imap++ ) {

            std::set<Hyd> these_hyds_set = early_map[imap];
            std::vector<Hyd> these_hyds;
            for ( Hyd hyd : these_hyds_set ) {
                these_hyds.push_back(hyd);
            }


            bool wrote_a_null = false;
            for ( size_t idx = 0; idx < cation_max_hyds_; idx++ ) {
                size_t offset = imap * cation_max_hyds_ + idx;
                if ( idx < these_hyds.size() ) {
                    cation_voxel_map_.at(offset) = these_hyds.at(idx);
                } else {
                    cation_voxel_map_.at(offset) = CACHE_MAX_HYD;
                    wrote_a_null = true;
                }
            }
            runtime_assert( wrote_a_null );     // Redundancy. Make sure each list is terminated
        }
    }


/////////////////////////////////////////////////////////////////////////

*/

    bool
    apply( 
        std::vector<std::pair<int, EigenXform> > const & irot_and_bbpos,
        std::vector<int> & hyd_counts,
        std::vector<int> & per_irot_counts,
        bool & pass_hydrophobic_residue_contacts,
        bool & pass_hydrophobic_ddg,
        bool & pass_better_than
        /*
        bool & pass_cation_pi
        */
    ) const {
        hyd_counts.clear();
        hyd_counts.resize( hydrophobic_res_.size(), 0 );

        per_irot_counts.clear();
        per_irot_counts.resize( irot_and_bbpos.size(), 0 );

        std::vector<int> has_ok_atoms( irot_and_bbpos.size(), 0 );

        typedef typename RotamerIndex::Atom Atom;

        for ( int ipair = 0; ipair < irot_and_bbpos.size(); ipair++ ) {
            std::pair<int, EigenXform> const & pair = irot_and_bbpos[ipair];
            int irot = pair.first;
            EigenXform const & bbpos = pair.second;


            std::tr1::unordered_map<Hyd, int> with_whom;

            int ok_atoms = 0;
            int this_irot_count = 0;
            for( int iatom = 3; iatom < rot_index_p->nheavyatoms(irot); ++iatom )
            {
                Atom const & atom = rot_index_p->rotamer(irot).atoms_.at(iatom);

                if ( ! rif_hydrophobic_map_.at(atom.type()) ) continue;
                ok_atoms ++;
                typename Atom::Position pos = bbpos * atom.position();

                std::vector<Hyd>::const_iterator hyds_iter = this->at( pos );


                Hyd this_hyd = 0;
                while ( (this_hyd = *(hyds_iter++)) != CACHE_MAX_HYD ) {
                    hyd_counts[this_hyd] ++;
                    this_irot_count++;
                    if (with_whom.count(this_hyd) == 0) {
                        with_whom[this_hyd] = 1;
                    } else {
                        with_whom[this_hyd] += 1;
                    }
                }
            }
            has_ok_atoms[ipair] = ok_atoms;

            int partners = 0;
            for ( auto pair : with_whom ) {
                if ( pair.second >= 2) {
                    partners += 1;
                }
            }

            per_irot_counts[ipair] = ( partners >= 2 ? this_irot_count : 0 );
        }

        int total_sum = 0;
        int hydrophobic_residue_contacts = 0;
        for ( int count : hyd_counts ) {
            if (count >= 2) {
                hydrophobic_residue_contacts++;
            }
            total_sum += count;
        }

// HYDROPHOBIC DDG
        float hydrophobic_ddg = DDG_PER_COUNT * total_sum; // This is how it was parameterized, but holy crap this is a hack.

        pass_hydrophobic_ddg = (hydrophobic_ddg <= hydrophobic_ddg_cut_ );

        pass_hydrophobic_residue_contacts = (hydrophobic_residue_contacts >= hydrophobic_residue_contacts_cut_);
// HYDROPHOBICS BETTER THAN

        // ohhhhhhhhhhhhhhhh..
        // my god, what's going on here?????????????????????
        pass_better_than = true;
        if ( doing_better_than_ ) {
            int one_better_than = 0;
            int two_better_than = 0;
            int three_better_than = 0;
            for ( int ipair = 0; ipair < per_irot_counts.size(); ipair ++) {
                float ddg_hyd = per_irot_counts[ipair] * DDG_PER_COUNT;
                float ddg_ratio = ddg_hyd / has_ok_atoms[ipair];
                if (ddg_ratio > hydrophobic_ddg_per_atom_cut_) continue;
                if ( ddg_hyd < one_hydrophobic_better_than_ ) one_better_than++;
                if ( ddg_hyd < two_hydrophobics_better_than_ ) two_better_than++;
                if ( ddg_hyd < three_hydrophobics_better_than_ ) three_better_than++;
            }
            if ( one_hydrophobic_better_than_ < 0 ) pass_better_than &= (one_better_than >= 1);
            if ( two_hydrophobics_better_than_ < 0 ) pass_better_than &= (two_better_than >= 2);
            if ( three_hydrophobics_better_than_ < 0 ) pass_better_than &= (three_better_than >= 3);
        }


        // CATION PI
        /*
        pass_cation_pi = true;
        if ( num_cation_pi_ > 0 ) {
            int found_cation_pi = 0;
            for ( int ipair = 0; ipair < irot_and_bbpos.size(); ipair++ ) {
                std::pair<intRot, EigenXform> const & pair = irot_and_bbpos[ipair];
                int irot = pair.first;
                EigenXform const & bbpos = pair.second;

                std::map<Hyd, int> this_irot_counts;
                int this_irot_count = 0;
                for( int iatom = 3; iatom < rot_index_p->nheavyatoms(irot); ++iatom )
                {
                    Atom const & atom = rot_index_p->rotamer(irot).atoms_.at(iatom);

                    if ( ! rif_pi_map_.at(atom.type()) ) continue;

                    typename Atom::Position pos = bbpos * atom.position();

                    std::vector<Hyd>::const_iterator hyds_iter = this->cat_at( pos );

                    Hyd this_hyd = 0;
                    while ( (this_hyd = *(hyds_iter++)) != CACHE_MAX_HYD ) {
                        if ( this_irot_counts.count( this_hyd ) == 0) {
                            this_irot_counts[ this_hyd ] = 1;
                        } else {
                            this_irot_counts[ this_hyd ] += 1;
                        }
                    }
                }
                for ( std::pair<Hyd, int> pair : this_irot_counts ) {
                    if ( pair.second >= 5 ) {       // 5 for histidine
                        found_cation_pi += 1;
                    }
                }
            }
            pass_cation_pi = found_cation_pi >= num_cation_pi_;
        }
        */


        return pass_hydrophobic_ddg && pass_hydrophobic_residue_contacts && pass_better_than;
    }

    void
    print_hydrophobic_counts( 
        core::pose::Pose const & target,
        std::vector<int> const & hyd_counts,
        std::vector<std::pair<int, EigenXform>> const & irot_and_bbpos,
        std::vector<int> const & seqposs,
        std::vector<int> const & per_irot_counts,
        int scaff_size,
        std::ostream & ostream
    ) const {
        using ObjexxFCL::format::F;
        ostream << "Hydrophobic counts:" << std::endl;

        runtime_assert( hyd_counts.size() == hydrophobic_res_.size() );

        for ( int i = 0; i < hyd_counts.size(); i++ ) {
            int count = hyd_counts[i];
            if ( count >= 2 ) {
                core::Size seqpos = hydrophobic_res_[i];
                core::conformation::Residue const & res = target.residue(seqpos);
                ostream << " " << seqpos << "/" << seqpos + scaff_size << " " << res.name3() << ": " << count << std::endl;
            }
        }

        if ( doing_better_than_ ) {
            ostream << "Hotspot ddGs:" << std::endl;
            runtime_assert( irot_and_bbpos.size() == seqposs.size() );
            runtime_assert( per_irot_counts.size() == seqposs.size() );

            for ( int i = 0; i < irot_and_bbpos.size(); i++ ) {
                int irot = irot_and_bbpos[i].first;
                ostream << " " << seqposs[i] << " " << rot_index_p->resname(irot) << ": " 
                          << F(5,1,per_irot_counts[i] * DDG_PER_COUNT) << std::endl;
            }
        }

    }




    template<class Floats> Indices floats_to_index(Floats const & f) const {
        Indices ind;
        for(int i = 0; i < 3; ++i){
            float tmp = ((f[i]-lb_[i])/cs_[i]);
            ind[i] = tmp;
        }
        return ind;
    }

    size_t index_to_map_index( Indices const & ind ) const {

        size_t accum = ind[0];
        accum = accum * shape_[1] + ind[1];
        accum = accum * shape_[2] + ind[2];

        return accum;

    }

    // flattening the 3D voxel
    size_t index_to_offset( Indices const & ind ) const {

        return index_to_map_index( ind ) * max_hyds_;

    }



    /*
    template<class Floats> Indices cation_floats_to_index(Floats const & f) const {
        Indices ind;
        for(int i = 0; i < 3; ++i){
            float tmp = ((f[i]-cat_lb_[i])/cat_cs_[i]);
            ind[i] = tmp;
        }
        return ind;
    }

    size_t cation_index_to_map_index( Indices const & ind ) const {

        size_t accum = ind[0];
        accum = accum * cation_shape_[1] + ind[1];
        accum = accum * cation_shape_[2] + ind[2];

        return accum;

    }

    size_t cation_index_to_offset( Indices const & ind ) const {

        return cation_index_to_map_index( ind ) * cation_max_hyds_;

    }
    */
    // template<class Floats>
    // std::vector<Hyd>::const_iterator
    // operator[](Floats const & floats) const { 
    //     Indices ind = floats_to_index(floats);
    //     return voxel_map_.cbegin() + index_to_offset(ind);
    // }

    // template<class Floats>
    // std::vector<Hyd>::iterator
    // operator[](Floats const & floats){
    //     Indices ind = floats_to_index(floats);
    //     return voxel_map_.begin() + index_to_offset(ind);
    // }

    std::vector<Hyd>::const_iterator 
    at( float f, float g, float h ) const {
        Indices idx = floats_to_index( Bounds( f, g, h ) );
        if( idx[0] < shape_[0] && idx[1] < shape_[1] && idx[2] < shape_[2] )
            return voxel_map_.cbegin() + index_to_offset(idx);
        else return OOB_LIST.cbegin();
    }

    template<class V>
    std::vector<Hyd>::const_iterator
    at( V const & v ) const {
        Indices idx = floats_to_index( Bounds( v[0], v[1], v[2] ) );
        if( idx[0] < shape_[0] && idx[1] < shape_[1] && idx[2] < shape_[2] )
            return voxel_map_.cbegin() + index_to_offset(idx);
        else return OOB_LIST.cbegin();
    }

    /*
    std::vector<Hyd>::const_iterator 
    cat_at( float f, float g, float h ) const {
        Indices idx = cation_floats_to_index( Bounds( f, g, h ) );
        if( idx[0] < cation_shape_[0] && idx[1] < cation_shape_[1] && idx[2] < cation_shape_[2] )
            return cation_voxel_map_.cbegin() + cation_index_to_offset(idx);
        else return OOB_LIST.cbegin();
    }

    template<class V>
    std::vector<Hyd>::const_iterator
    cat_at( V const & v ) const {
        Indices idx = cation_floats_to_index( Bounds( v[0], v[1], v[2] ) );
        if( idx[0] < cation_shape_[0] && idx[1] < cation_shape_[1] && idx[2] < cation_shape_[2] )
            return cation_voxel_map_.cbegin() + cation_index_to_offset(idx);
        else return OOB_LIST.cbegin();
    }
    */



};




}}

#endif
