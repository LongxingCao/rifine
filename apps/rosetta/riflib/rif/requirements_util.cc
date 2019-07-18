// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include "requirements_util.hh"



#include <utility/file/file_sys_util.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

#include <utility/string_util.hh>

namespace devel {
namespace scheme {
namespace rif {
    
    // reture the donor definitions from a tuning file.
    std::vector< DonorDefinition > get_donor_definitions( std::string tuning_file )
    {
        std::vector< DonorDefinition > donors;
        
        if ( tuning_file == "" )
        {
            return donors;
        }
        runtime_assert_msg(utility::file::file_exists( tuning_file ), "tunning file does not exits: " + tuning_file );
        std::ifstream in;
        std::string s;
        in.open( tuning_file , std::ios::in );
        std::vector<std::string> lines;
        bool flag = false;
        while ( std::getline(in, s) ){
            if (s.empty() || s.find("#") == 0) continue;
            if (s.find("DONOR_DEFINITION") != std::string::npos && s.find("END_DONOR_DEFINITION") == std::string::npos ) { flag = true; continue; }
            else if (s.find("END_DONOR_DEFINITION") != std::string::npos ) { flag = false; break; }
            
            if ( flag )
            {
                DonorDefinition donor_temp;
                utility::vector1<std::string> splt = utility::quoted_split( s );
                runtime_assert_msg(splt.size() >=1, "something is wrong with the donor definition block, please check the tuning file." );
                donor_temp.res_num = utility::string2int( splt[1] );
                donor_temp.allowed_donor_res.clear();
                for(int ii = 2; ii <= splt.size(); ++ii )
                {
                    donor_temp.allowed_donor_res.push_back( splt[ii] );
                }
                donors.push_back(donor_temp);
            }
        }
        return donors;
    }
    
    // reture the hbond definitions from a tuning file.
    std::vector< AcceptorDefinition > get_acceptor_definitions( std::string tuning_file )
    {
        std::vector< AcceptorDefinition > acceptors;
        
        if ( tuning_file == "" )
        {
            return acceptors;
        }
        runtime_assert_msg(utility::file::file_exists( tuning_file ), "tunning file does not exits: " + tuning_file );
        std::ifstream in;
        std::string s;
        in.open( tuning_file , std::ios::in );
        std::vector<std::string> lines;
        bool flag = false;
        while ( std::getline(in, s) ){
            if (s.empty() || s.find("#") == 0) continue;
            if (s.find("ACCEPTOR_DEFINITION") != std::string::npos && s.find("END_ACCEPTOR_DEFINITION") == std::string::npos ) { flag = true; continue; }
            else if (s.find("END_ACCEPTOR_DEFINITION") != std::string::npos ) { flag = false; break; }
            
            if ( flag )
            {
                AcceptorDefinition acceptor_temp;
                utility::vector1<std::string> splt = utility::quoted_split( s );
                runtime_assert_msg(splt.size() >=1, "something is wrong with the acceptor definition block, please check the tuning file." );
                acceptor_temp.res_num = utility::string2int( splt[1] );
                acceptor_temp.allowed_acceptor_res.clear();
                for(int ii = 2; ii <= splt.size(); ++ii )
                {
                    acceptor_temp.allowed_acceptor_res.push_back( splt[ii] );
                }
                acceptors.push_back(acceptor_temp);
            }
        }
        return acceptors;
    }


    
    // reture the hbond definitions from a tuning file.
    std::vector< HBondDefinition > get_hbond_definitions( std::string tuning_file )
    {
        std::vector< HBondDefinition > hbs;
        
        if ( tuning_file == "" )
        {
            return hbs;
        }
        runtime_assert_msg(utility::file::file_exists( tuning_file ), "tunning file does not exits: " + tuning_file );
        std::ifstream in;
        std::string s;
        in.open( tuning_file , std::ios::in );
        std::vector<std::string> lines;
        bool flag = false;
        while ( std::getline(in, s) ){
            if (s.empty() || s.find("#") == 0) continue;
            if (s.find("HBOND_DEFINITION") != std::string::npos && s.find("END_HBOND_DEFINITION") == std::string::npos ) { flag = true; continue; }
            else if (s.find("END_HBOND_DEFINITION") != std::string::npos ) { flag = false; break; }
            
            if ( flag )
            {
                HBondDefinition hb_temp;
                utility::vector1<std::string> splt = utility::quoted_split( s );
                runtime_assert_msg(splt.size() >=2, "something is wrong with the hbond definition block, please check the tuning file." );
                hb_temp.atom_name = splt[1];
                hb_temp.res_num = utility::string2int( splt[2] );
                hb_temp.allowed_rot_names.clear();
                for(int ii = 3; ii <= splt.size(); ++ii )
                {
                    hb_temp.allowed_rot_names.push_back( splt[ii] );
                }
                hbs.push_back(hb_temp);
            }
        }
        return hbs;
    }
    
    std::vector< BidentateDefinition > get_bidentate_definitions( std::string tuning_file )
    {
        std::vector< BidentateDefinition > bdhbs;
        BidentateDefinition bdhb_temp;
        
        if ( tuning_file == "" ) {
            return bdhbs;
        }
        runtime_assert_msg(utility::file::file_exists( tuning_file ), "tunning file does not exits: " + tuning_file );
        std::ifstream in;
        std::string s;
        in.open( tuning_file , std::ios::in );
        std::vector<std::string> lines;
        bool flag = false;
        while ( std::getline(in, s) ){
            if (s.empty() || s.find("#") == 0) continue;
            if (s.find("BIDENTATE_DEFINITION") != std::string::npos && s.find("END_BIDENTATE_DEFINITION") == std::string::npos ) { flag = true; continue; }
            else if (s.find("END_BIDENTATE_DEFINITION") != std::string::npos ) { flag = false; break; }
            
            if ( flag )
            {
                utility::vector1<std::string> splt = utility::quoted_split( s );
                runtime_assert_msg(splt.size() == 4, "something is wrong with the bidentate hydrogen bonds definition block, please check the tuning file." );
                bdhb_temp.atom1_name = splt[1];
                bdhb_temp.res1_num = utility::string2int( splt[2] );
                bdhb_temp.atom2_name = splt[3];
                bdhb_temp.res2_num = utility::string2int( splt[4] );
                bdhbs.push_back(bdhb_temp);
            }
        }
        return bdhbs;
    }

    bool check_requirement_definition_exists( std::string tuning_file )
    {
        // TODO: there are lots of duplicate codes here, it is not smart, but I need to make it smarter.
        if ( tuning_file == "" ) {
            return false;
        }
        runtime_assert_msg(utility::file::file_exists( tuning_file ), "tunning file does not exits: " + tuning_file );
        std::ifstream in;
        std::string s;
        in.open( tuning_file , std::ios::in );
        std::vector<std::string> lines;
        
        bool flag = false;
        bool any_definition = false;
        while ( std::getline(in, s) ){
            if (s.empty() || s.find("#") == 0) continue;
            if (s.find("REQUIREMENT_DEFINITION") != std::string::npos && s.find("END_REQUIREMENT_DEFINITION") == std::string::npos ) { flag = true; continue; }
            else if (s.find("END_REQUIREMENT_DEFINITION") != std::string::npos ) { flag = false; break; }
            
            // if there are some lines of definitions, return true.
				    if ( flag ) any_definition = true;
        }
        return any_definition;
    }
    
    std::vector< HbondRequirement > get_hbond_requirement_definitions( std::string tuning_file )
    {
        std::vector< HbondRequirement > hbond_reqs;
        
        if ( tuning_file == "" ) {
            return hbond_reqs;
        }
        runtime_assert_msg(utility::file::file_exists( tuning_file ), "tunning file does not exits: " + tuning_file );
        std::ifstream in;
        std::string s;
        in.open( tuning_file , std::ios::in );
        std::vector<std::string> lines;
        bool flag = false;
        while ( std::getline(in, s) ){
            if (s.empty() || s.find("#") == 0) continue;
            if (s.find("REQUIREMENT_DEFINITION") != std::string::npos && s.find("END_REQUIREMENT_DEFINITION") == std::string::npos ) { flag = true; continue; }
            else if (s.find("END_REQUIREMENT_DEFINITION") != std::string::npos ) { flag = false; break; }
            
            if ( flag )
            {
                HbondRequirement req_temp;
                utility::vector1<std::string> splt = utility::quoted_split( s );
                
                //std::cout << req_temp.req_num << std::endl;
                
                if ( splt.size() == 4 && splt[2] == "HBOND" ) {
                    runtime_assert_msg( splt.size() == 4, "something is wrong with the HBOND requirement definition!" );
                    runtime_assert_msg( utility::string2int( splt[1] ) >=0, "The requirement number must be a positive integer!" );
                    req_temp.req_num = utility::string2int(splt[1]);
                    req_temp.atom_name = splt[3];
                    req_temp.res_num = utility::string2int( splt[4] );
                    hbond_reqs.push_back(req_temp);
                }
            }
        }
        return hbond_reqs;
    }
    
    std::vector< BidentateRequirement > get_bidentate_requirement_definitions( std::string tuning_file )
    {
        std::vector< BidentateRequirement > bidentate_reqs;
        
        if ( tuning_file == "" ) {
            return bidentate_reqs;
        }
        runtime_assert_msg(utility::file::file_exists( tuning_file ), "tunning file does not exits: " + tuning_file );
        std::ifstream in;
        std::string s;
        in.open( tuning_file , std::ios::in );
        std::vector<std::string> lines;
        bool flag = false;
        while ( std::getline(in, s) ){
            if (s.empty() || s.find("#") == 0) continue;
            if (s.find("REQUIREMENT_DEFINITION") != std::string::npos && s.find("END_REQUIREMENT_DEFINITION") == std::string::npos ) { flag = true; continue; }
            else if (s.find("END_REQUIREMENT_DEFINITION") != std::string::npos ) { flag = false; break; }
            
            if ( flag )
            {
                BidentateRequirement req_temp;
                utility::vector1<std::string> splt = utility::quoted_split( s );
                
                //std::cout << req_temp.req_num << std::endl;
                
                if ( splt.size() == 6 && splt[2] == "BIDENTATE" ) {
                    runtime_assert_msg( splt.size() == 6, "something is wrong with the BIDENTATE requirement definition!" );
                    runtime_assert_msg( utility::string2int( splt[1] ) >=0, "The requirement number must be a positive integer!" );
                    req_temp.req_num = utility::string2int(splt[1]);
                    req_temp.atom1_name = splt[3];
                    req_temp.res1_num = utility::string2int( splt[4] );
                    req_temp.atom2_name = splt[5];
                    req_temp.res2_num = utility::string2int( splt[6] );
                    bidentate_reqs.push_back(req_temp);
                }
            }
        }
        return bidentate_reqs;
    }

    std::vector< TridentateRequirement > get_tridentate_requirement_definitions( std::string tuning_file )
    {
        std::vector< TridentateRequirement > tridentate_reqs;
        
        if ( tuning_file == "" ) {
            return tridentate_reqs;
        }
        runtime_assert_msg(utility::file::file_exists( tuning_file ), "tunning file does not exits: " + tuning_file );
        std::ifstream in;
        std::string s;
        in.open( tuning_file , std::ios::in );
        std::vector<std::string> lines;
        bool flag = false;
        while ( std::getline(in, s) ){
            if (s.empty() || s.find("#") == 0) continue;
            if (s.find("REQUIREMENT_DEFINITION") != std::string::npos && s.find("END_REQUIREMENT_DEFINITION") == std::string::npos ) { flag = true; continue; }
            else if (s.find("END_REQUIREMENT_DEFINITION") != std::string::npos ) { flag = false; break; }
            
            if ( flag )
            {
                TridentateRequirement req_temp;
                utility::vector1<std::string> splt = utility::quoted_split( s );
                
                //std::cout << req_temp.req_num << std::endl;
                
                if ( splt.size() == 8 && splt[2] == "TRIDENTATE" ) {
                    runtime_assert_msg( splt.size() == 8, "something is wrong with the BIDENTATE requirement definition!" );
                    runtime_assert_msg( utility::string2int( splt[1] ) >=0, "The requirement number must be a positive integer!" );
                    req_temp.req_num = utility::string2int(splt[1]);
                    req_temp.atom1_name = splt[3];
                    req_temp.res1_num = utility::string2int( splt[4] );
                    req_temp.atom2_name = splt[5];
                    req_temp.res2_num = utility::string2int( splt[6] );
                    req_temp.atom3_name = splt[7];
                    req_temp.res3_num = utility::string2int( splt[8] );
                    tridentate_reqs.push_back(req_temp);
                }
            }
        }
        return tridentate_reqs;
    }
    
    std::vector< HotspotRequirement > get_hotspot_requirement_definitions( std::string tuning_file )
    {
        std::vector< HotspotRequirement > hotspot_reqs;
        
        if ( tuning_file == "" ) {
            return hotspot_reqs;
        }
        runtime_assert_msg(utility::file::file_exists( tuning_file ), "tunning file does not exits: " + tuning_file );
        std::ifstream in;
        std::string s;
        in.open( tuning_file , std::ios::in );
        std::vector<std::string> lines;
        bool flag = false;
        while ( std::getline(in, s) ){
            if (s.empty() || s.find("#") == 0) continue;
            if (s.find("REQUIREMENT_DEFINITION") != std::string::npos && s.find("END_REQUIREMENT_DEFINITION") == std::string::npos ) { flag = true; continue; }
            else if (s.find("END_REQUIREMENT_DEFINITION") != std::string::npos ) { flag = false; break; }
            
            if ( flag )
            {
                HotspotRequirement req_temp;
                utility::vector1<std::string> splt = utility::quoted_split( s );
                
                //std::cout << req_temp.req_num << std::endl;
                
                if ( splt.size() == 3 && splt[2] == "HOTSPOT" ) {
                    runtime_assert_msg( splt.size() == 3, "something is wrong with the HOTSPOT requirement definition!" );
                    runtime_assert_msg( utility::string2int( splt[1] ) >=0, "The requirement number must be a positive integer!" );
                    req_temp.req_num = utility::string2int(splt[1]);
                    req_temp.hotspot_num = utility::string2int(splt[3]);
                    hotspot_reqs.push_back(req_temp);
                }
            }
        }
        return hotspot_reqs;
    }
    
    std::vector< ApoRequirement > get_Apo_requirement_definitions( std::string tuning_file )
    {
        std::vector< ApoRequirement > apo_reqs;
        
        if ( tuning_file == "" ) {
            return apo_reqs;
        }
        runtime_assert_msg(utility::file::file_exists( tuning_file ), "tunning file does not exits: " + tuning_file );
        std::ifstream in;
        std::string s;
        in.open( tuning_file , std::ios::in );
        std::vector<std::string> lines;
        
        // define here, because of the logic ...
        // how to parse the apolar definition
        ApoRequirement req_temp;
        bool flag = false;
        bool found_one = false;
        while ( std::getline(in, s) ){
            if (s.empty() || s.find("#") == 0) continue;
            if (s.find("REQUIREMENT_DEFINITION") != std::string::npos && s.find("END_REQUIREMENT_DEFINITION") == std::string::npos ) { flag = true; continue; }
            else if (s.find("END_REQUIREMENT_DEFINITION") != std::string::npos ) { flag = false; break; }
            
            if ( flag && !found_one )
            {
                
                utility::vector1<std::string> splt = utility::quoted_split( s );
                if ( splt.size() >=2 && splt[2] == "APOLAR" ) {
                    runtime_assert_msg( utility::string2int( splt[1] ) >=0, "The requirement number must be a positive integer!" );
										req_temp.req_num = utility::string2int( splt[1] );
                    for (int ii = 3; ii <= splt.size(); ++ii) {
                        req_temp.allowed_rot_names.push_back( splt[ii] );
                    }
                    found_one = true;
                    continue;
                }
            } else if ( flag && found_one ) {
                utility::vector1<std::string> splt = utility::quoted_split( s );
                if ( splt.size() == 3 ) {
                    ApoReqTerm term_temp;
                    term_temp.atom_name = splt[1];
                    term_temp.res_num = utility::string2int( splt[2] );
                    term_temp.distance = utility::string2float( splt[3] );
										req_temp.terms.push_back( term_temp );
                } else if ( splt.size() == 1 ) {
                    runtime_assert_msg( splt[1] == "END_APOLAR", "something is wrong with the APOLAR requirement definition!" );
                    runtime_assert_msg( req_temp.terms.size() != 0, "You must assign some definitions for the apolar req, talk with longxing");
                    apo_reqs.push_back(req_temp);
										// reinitialize the req_temp
										req_temp.req_num = -1;
										req_temp.allowed_rot_names.clear();
										req_temp.terms.clear();
                    found_one = false;
                } else {
                    utility_exit_with_message("something is wrong with the APOLAR requirement definition!" );
                }
            } else {
                // not in the
            }
        }
        return apo_reqs;
    }
    
    std::vector< PiPiStackingRequirement > get_pipi_stacking_requirement_definitions( std::string tuning_file )
    {
        std::vector< PiPiStackingRequirement > pipi_reqs;
        
        if ( tuning_file == "" ) {
            return pipi_reqs;
        }
        runtime_assert_msg(utility::file::file_exists( tuning_file ), "tunning file does not exits: " + tuning_file );
        std::ifstream in;
        std::string s;
        in.open( tuning_file , std::ios::in );
        std::vector<std::string> lines;
        
        std::vector< std::string > residues_allowed_for_pipi_stacking;
        residues_allowed_for_pipi_stacking.push_back( "TRP" );
        residues_allowed_for_pipi_stacking.push_back( "PHE" );
        residues_allowed_for_pipi_stacking.push_back( "TYR" );
        //residues_allowed_for_pipi_stacking.push_back( "HIS" );
        
        // define here, because of the logic ...
        // how to parse the apolar definition
        PiPiStackingRequirement req_temp;
        bool flag = false;
        bool found_one = false;
        while ( std::getline(in, s) ){
            if (s.empty() || s.find("#") == 0) continue;
            if (s.find("REQUIREMENT_DEFINITION") != std::string::npos && s.find("END_REQUIREMENT_DEFINITION") == std::string::npos ) { flag = true; continue; }
            else if (s.find("END_REQUIREMENT_DEFINITION") != std::string::npos ) { flag = false; break; }
            
            if ( flag && !found_one )
            {
                
                utility::vector1<std::string> splt = utility::quoted_split( s );
                if ( splt.size() >=2 && splt[2] == "PIPISTACKING" ) {
                    runtime_assert_msg( utility::string2int( splt[1] ) >=0, "The requirement number must be a positive integer!" );
                    req_temp.req_num = utility::string2int( splt[1] );
                    if ( splt.size() == 2 ) {
                        req_temp.allowed_rot_names = residues_allowed_for_pipi_stacking;
                    } else {
                        for (int ii = 3; ii <= splt.size(); ++ii) {
                            req_temp.allowed_rot_names.push_back( splt[ii] );
                            runtime_assert_msg( std::find( residues_allowed_for_pipi_stacking.begin(), residues_allowed_for_pipi_stacking.end(), splt[ii] ) != residues_allowed_for_pipi_stacking.end(), "Can this residue form pipi stacking interactin?? How?? Let me know." );
                        }
                    }
                    found_one = true;
                    continue;
                }
            } else if ( flag && found_one ) {
                utility::vector1<std::string> splt = utility::quoted_split( s );
                if ( splt.size() == 2 ) {
                    PiPiStackingTerm term_temp;
                    term_temp.atom_name = splt[1];
                    term_temp.res_num = utility::string2int( splt[2] );
                    req_temp.terms.push_back( term_temp );
                } else if ( splt.size() == 1 ) {
                    runtime_assert_msg( splt[1] == "END_PIPISTACKING", "something is wrong with the PIPI stacking requirement definition!" );
                    runtime_assert_msg( req_temp.terms.size() >= 3, "You must assign at least three points to define the ring, talk with longxing to fix the code if you have better ideas");
                    pipi_reqs.push_back(req_temp);
                    // reinitialize the req_temp
                    req_temp.req_num = -1;
                    req_temp.allowed_rot_names.clear();
                    req_temp.terms.clear();
                    found_one = false;
                } else {
                    utility_exit_with_message("something is wrong with the APOLAR requirement definition!" );
                }
            } else {
                // not in the
            }
        }
        return pipi_reqs;
    }

    
    std::vector< CationPiRequirement > get_cationpi_requirement_definitions( std::string tuning_file )
    {
        std::vector< CationPiRequirement > cationpi_reqs;
        
        if ( tuning_file == "" ) {
            return cationpi_reqs;
        }
        
        std::vector< std::string > residues_allowed_for_cationpi;
        residues_allowed_for_cationpi.push_back( "TRP" );
        residues_allowed_for_cationpi.push_back( "PHE" );
        residues_allowed_for_cationpi.push_back( "TYR" );
        residues_allowed_for_cationpi.push_back( "HIS" );
        
        runtime_assert_msg(utility::file::file_exists( tuning_file ), "tunning file does not exits: " + tuning_file );
        std::ifstream in;
        std::string s;
        in.open( tuning_file , std::ios::in );
        std::vector<std::string> lines;
        bool flag = false;
        while ( std::getline(in, s) ){
            if (s.empty() || s.find("#") == 0) continue;
            if (s.find("REQUIREMENT_DEFINITION") != std::string::npos && s.find("END_REQUIREMENT_DEFINITION") == std::string::npos ) { flag = true; continue; }
            else if (s.find("END_REQUIREMENT_DEFINITION") != std::string::npos ) { flag = false; break; }
            
            if ( flag )
            {
                CationPiRequirement req_temp;
                utility::vector1<std::string> splt = utility::quoted_split( s );
                
                //std::cout << req_temp.req_num << std::endl;
                
                if ( splt.size() >= 2 &&  splt[2] == "CATIONPI" ) {
                    runtime_assert_msg( splt.size() >= 3, "something is wrong with the CATIONPI requirement definition, talk with longxing about this. You can specify the allow residues to form a cationpi, such as TRP, PHE, HIS!" );
                    runtime_assert_msg( utility::string2int( splt[1] ) >=0, "The requirement number must be a positive integer!" );
                    req_temp.req_num = utility::string2int(splt[1]);
                    req_temp.res_num = utility::string2int(splt[3]);
                    if ( splt.size() == 3 ) {
                        req_temp.allowed_rot_names = residues_allowed_for_cationpi;
                    } else {
                        for(int ii = 4; ii <= splt.size(); ++ii )
                        {
                            runtime_assert_msg( std::find( residues_allowed_for_cationpi.begin(), residues_allowed_for_cationpi.end(), splt[ii] ) != residues_allowed_for_cationpi.end(), "can this residue form cationpi interactin?? we should talk about this." );
                            req_temp.allowed_rot_names.push_back( splt[ii] );
                        }
                    }
                    cationpi_reqs.push_back(req_temp);
                }
            }
        }
        return cationpi_reqs;
    }
    
}
}
}
