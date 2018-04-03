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
    
    std::vector< RequirementDefinition > get_requirement_definitions( std::string tuning_file )
    {
        std::vector< RequirementDefinition > reqs;
        
        if ( tuning_file == "" ) {
            return reqs;
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
                RequirementDefinition req_temp;
                utility::vector1<std::string> splt = utility::quoted_split( s );
                runtime_assert_msg( utility::string2int( splt[1] ) >=0, "The requirement number must be a positive integer!" );
                req_temp.req_num = utility::string2int( splt[1] );
                req_temp.require = splt[2];
                for (int ii = 3; ii <= splt.size(); ++ii) {
                    req_temp.definition.push_back( splt[ii] );
                }
                
                
                //std::cout << req_temp.req_num << std::endl;
                
                
                if ( req_temp.require == "HBOND" ) {
                    runtime_assert_msg( req_temp.definition.size() == 2, "something is wrong with the HBOND requirement definition!" );
                } else if ( req_temp.require == "BIDENTATE" ) {
                    runtime_assert_msg( req_temp.definition.size() == 4, "something is wrong with the BIDENTATE requirement definition!" );
                } else if ( req_temp.require == "HOTSPOT" ) {
                    runtime_assert_msg( req_temp.definition.size() == 1, "something is wrong with the HOTSPOT requirement definition!" );
                } else {
                    utility_exit_with_message("Unknown requirement definition: " + std::string(req_temp.require) );
                }
                
                reqs.push_back(req_temp);
            }
        }
        return reqs;
    }
    
    
    
}
}
}