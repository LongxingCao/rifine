// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols

#ifndef INCLUDED_riflib_scaffold_morph_util_hh
#define INCLUDED_riflib_scaffold_morph_util_hh

#include <riflib/types.hh>
#include <core/types.hh>
#include <rif_dock_test.hh>
#include <riflib/rifdock_typedefs.hh>
#include <scheme/scaffold/ScaffoldProviderBase.hh>

#include <utility/string_util.hh>
#include <utility/file/file_sys_util.hh>


namespace devel {
namespace scheme {



struct MorphAction {
    core::Size first_seq_pos;
    core::Size last_seq_pos;
    core::Size frag_length;
    std::string frag_id;
};




struct MorphRule {
    uint64_t low_cut_site;
    uint64_t high_cut_site;
    uint64_t max_deletion;
    uint64_t max_insertion;
    uint64_t max_fragments;
    float fragment_cluster_tolerance;
    float fragment_max_rmsd;
};



typedef std::vector<MorphRule> MorphRules;


struct MorphMember {
    ParametricSceneConformationCOP conformation;
    MorphRules morph_rules;
    std::vector<MorphRule> morph_history;
    ::scheme::scaffold::TreeRelation tree_relation;
};




inline
MorphRule
morph_rule_from_options(RifDockOpt const & opt) {
    MorphRule r;

    r.low_cut_site = opt.low_cut_site;
    r.high_cut_site = opt.high_cut_site;
    r.max_deletion = opt.max_deletion;
    r.max_insertion = opt.max_insertion;
    r.max_fragments = opt.max_fragments;
    r.fragment_cluster_tolerance = opt.fragment_cluster_tolerance;
    r.fragment_max_rmsd = opt.fragment_max_rmsd;

    return r;
}



// S entry
// low_cut high_cut max_deletion max_insertion num_structures cluster_tolerance max_rmsd

inline
bool
parse_morph_rules_file(std::string fname, MorphRules & rules, RifDockOpt const & opt) {

    runtime_assert_msg(utility::file::file_exists( fname ), "morph_rules file does not exist: " + fname );
    std::ifstream in;
    in.open(fname, std::ios::in);
    // utility::io::izstream in(fname);
    std::string s;
    while (std::getline(in, s)) {
        std::string save_s = s;
        if (s.empty()) continue;

        s = utility::replace_in( s, ":", " ");
        s = utility::replace_in( s, "#", " #");
        utility::vector1<std::string> comment_splt = utility::string_split_simple(s, ' ');
        utility::vector1<std::string> splt;
        for ( std::string item : comment_splt ) {
            item = utility::strip(item, " \t\n");
            if (item.empty()) continue;
            if (item[0] == '#') break;
            splt.push_back(item);
        }
        if (splt.size() == 0) continue;

        std::cout << "Parsing morph rules: " + save_s << std::endl;
        // std::cout << splt << std::endl; 
        if ( "S" == splt[1]) {
            MorphRule r;
            if (splt.size() < 5) {
                std::cout << "Bad line format: " << save_s << std::endl;
                return false;
            }
            r.low_cut_site = utility::string2int(splt[2]);
            r.high_cut_site = utility::string2int(splt[3]);
            r.max_deletion = utility::string2int(splt[4]);
            r.max_insertion = utility::string2int(splt[5]);
            r.max_fragments = splt.size() >= 6 ? utility::string2int(splt[6]) : opt.max_fragments;
            r.fragment_cluster_tolerance = splt.size() >= 7 ? utility::string2float(splt[7]) : opt.fragment_cluster_tolerance;
            r.fragment_max_rmsd = splt.size() >= 8 ? utility::string2float(splt[8]) : opt.fragment_max_rmsd;
            rules.push_back(r);

        }else{
            std::cout << "Error parsing line: " << save_s << std::endl;
            return false;
        }
    }

    if (rules.size() == 0) {
        std::cout << "Error, no morph rules in: " << fname << std::endl;
        return false;
    }
    
    return true;
}



}}

#endif
