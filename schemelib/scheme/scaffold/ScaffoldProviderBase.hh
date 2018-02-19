// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols

#ifndef INCLUDED_scaffold_ScaffoldProviderBase_hh
#define INCLUDED_scaffold_ScaffoldProviderBase_hh

#include <scheme/types.hh>

#include <string>
#include <vector>
#include <limits>
#include <boost/any.hpp>
#include <boost/functional/hash.hpp>

#include <core/pose/Pose.hh>


namespace scheme {
namespace scaffold {

const uint64_t BOGUS_INDEX = std::numeric_limits<uint64_t>::max();

template<
    class _Scaffold,
    class _ScaffoldIndex,
    class _ScaffoldIndexLimits
>
struct ScaffoldProviderBase {
    typedef _Scaffold Scaffold;
    typedef _ScaffoldIndex ScaffoldIndex;
    typedef _ScaffoldIndexLimits ScaffoldIndexLimits;
    typedef typename Scaffold::CacheData ScaffoldCacheData;
    ScaffoldProviderBase() {}


    virtual shared_ptr<Scaffold const> get_scaffold(ScaffoldIndex i) = 0;

    virtual ScaffoldIndexLimits get_scaffold_index_limits() const = 0;

    virtual shared_ptr<ScaffoldCacheData> get_data_cache_slow(ScaffoldIndex i) = 0;

    virtual void set_fa_mode( bool is_fa ) = 0;

    virtual ScaffoldIndex get_representative_scaffold_index() = 0;

    virtual void setup_twobody_tables( ScaffoldIndex i ) = 0;

    virtual void modify_pose_for_output( ScaffoldIndex i, core::pose::Pose & pose ) {}

};


// Key Assumptions of this class:


struct TreeIndex {
    uint64_t depth;
    uint64_t member;
    TreeIndex() : depth(BOGUS_INDEX), member(BOGUS_INDEX) {};
    TreeIndex(uint64_t _depth, uint64_t _member) : depth(_depth), member(_member) {};

    bool operator==(const TreeIndex &other) const { 
        return (depth == other.depth
            && member == other.member);
    }

};

inline
std::ostream & operator << ( std::ostream & out, TreeIndex const & ti ){
    out << "TreeIndex: " << ti.depth << " " << ti.member;
    return out;
}

struct TreeRelation {
    uint64_t depth;
    uint64_t parent_member;
    uint64_t first_child;
    uint64_t last_child;
};

typedef std::vector<uint64_t> TreeLimits;


template<
    class _Scaffold
>
struct TreeScaffoldProvider :
    public ScaffoldProviderBase<_Scaffold, TreeIndex, TreeLimits> {
    typedef _Scaffold Scaffold;


    virtual shared_ptr<Scaffold const> get_scaffold(TreeIndex i) = 0;

    // virtual void fill_children(TreeIndex i) = 0;

    virtual TreeLimits get_scaffold_index_limits() const = 0;


protected:

    bool
    is_valid_index(TreeIndex i) const {
        TreeLimits limits = get_scaffold_index_limits();

        if ( i.depth >= limits.size() )
            return false;

        if ( i.member >= limits[i.depth] )
            return false;

        return true;
    }

};



inline
uint64_t scaffold_index_default_value(uint64_t) {
    return BOGUS_INDEX;
}

inline
TreeIndex scaffold_index_default_value(TreeIndex) {
    return TreeIndex();
}




}}





namespace std {

    template <>
    struct hash<scheme::scaffold::TreeIndex>
    {
        std::size_t operator()(const scheme::scaffold::TreeIndex& ti) const {
            using std::size_t;
            using boost::hash;
            using boost::hash_combine;

            std::size_t seed = 0;

            boost::hash<int> hasher;
            hash_combine(seed, hasher(ti.depth));
            hash_combine(seed, hasher(ti.member));

            return seed;
        }
    };

}



#endif
