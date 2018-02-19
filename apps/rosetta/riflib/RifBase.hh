// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols

#ifndef INCLUDED_riflib_RifBase_hh
#define INCLUDED_riflib_RifBase_hh

#include <riflib/types.hh>

#include <string>
#include <vector>
#include <boost/any.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <riflib/RotamerGenerator.hh>

#include <core/conformation/Residue.hh>



namespace scheme { namespace search { struct HackPackOpts; }}

namespace devel {
namespace scheme {

template<class Index>
struct KeyIterHelperBase {
    virtual Index get_key() const = 0;
    virtual void next() = 0;
    virtual bool equal(KeyIterHelperBase const & that) const = 0;
};

struct RifBaseKeyIter : public boost::iterator_facade<
    RifBaseKeyIter, uint64_t, boost::forward_traversal_tag, uint64_t
>{
    RifBaseKeyIter() : irange_(nullptr) {}
    RifBaseKeyIter(std::shared_ptr<KeyIterHelperBase<uint64_t>> irange) : irange_(irange) {}
 private:
    friend class boost::iterator_core_access;
    std::shared_ptr<KeyIterHelperBase<uint64_t>> irange_;
    void increment() { irange_->next(); }
    bool equal(RifBaseKeyIter const& that) const { return this->irange_->equal(*that.irange_); }
    uint64_t dereference() const { return irange_->get_key(); }
};

struct RifBaseKeyRange {
    typedef RifBaseKeyIter const_iterator;
    const_iterator begin_, end_;
    RifBaseKeyRange(const_iterator begin, const_iterator end) : begin_(begin), end_(end) {}
    const_iterator begin() const { return begin_; }
    const_iterator end() const { return end_; }
};

// template<class _Index=uint64_t>
struct RifBase
{
    typedef uint64_t Key;
	RifBase() : type_("") {}
	RifBase( std::string type ) : type_(type){}
	std::string type() const { return type_; }
	std::string type_;

    virtual Key get_bin_key(EigenXform const & x) const = 0;
    virtual EigenXform get_bin_center( Key const & k) const = 0;
    virtual void get_rotamers_for_key( Key const & k, std::vector< std::pair< float, int > > & rotscores ) const = 0;
    virtual void get_rotamers_for_xform( EigenXform const & x, std::vector< std::pair< float, int > > & rotscores ) const = 0;
	virtual void get_rotamer_ids_in_use( std::vector<bool> & using_rot ) const = 0;
	virtual void print( std::ostream & out ) const = 0;
    // virtual void super_print( std::ostream & out, shared_ptr< RotamerIndex > rot_index_p ) const = 0;
	virtual void collision_analysis( std::ostream & out ) const = 0;
	virtual std::string value_name() const = 0;

    //Brian
    virtual std::pair< int, int > get_sat1_sat2( EigenXform const & x, int roti ) const {
        return std::pair< int, int >( -1, -1);
    }

	virtual size_t size() const = 0;
	virtual float  load_factor() const = 0;
	virtual float  cart_resl() const = 0;
	virtual float  ang_resl() const = 0;
	virtual size_t mem_use() const = 0;
	virtual size_t sizeof_value_type() const = 0;
	virtual bool  has_sat_data_slots() const = 0;

	template< class XMap > bool get_xmap_ptr( shared_ptr<XMap> & xmap_ptr );
	template< class XMap > bool get_xmap_const_ptr( shared_ptr<XMap const> & xmap_ptr ) const;
	template< class XMap > bool set_xmap_ptr( shared_ptr<XMap> const & xmap_ptr );
	virtual bool set_xmap_ptr( boost::any * any_p ) = 0;
	virtual bool get_xmap_ptr( boost::any * any_p ) = 0;
	virtual bool get_xmap_const_ptr( boost::any * any_p ) const = 0;

	virtual bool load( std::istream & in , std::string & description ) = 0;
	virtual bool save( std::ostream & out, std::string & description ) = 0;

	virtual void finalize_rif() = 0;

    virtual RifBaseKeyRange key_range() const = 0;

    virtual bool dump_rotamers_near_res( core::conformation::Residue const & res, std::string const & file_name, 
                                           float dump_dist, float dump_frac, shared_ptr<RotamerIndex> rot_index_p ) const = 0;

    // for convenience
    std::vector< std::pair< float, int > >
    get_rotamers_for_key(Key k) const {
        std::vector< std::pair< float, int > > rotscores;
        get_rotamers_for_key(k, rotscores);
        return rotscores;
    }
    std::vector< std::pair< float, int > >
    get_rotamers_for_xform(EigenXform const & x) const {
        std::vector< std::pair< float, int > > rotscores;
        get_rotamers_for_xform(x, rotscores);
        return rotscores;
    }

};

typedef shared_ptr<RifBase> RifPtr;
typedef shared_ptr<RifBase const> RifConstPtr;

template< class XMap >
bool RifBase::get_xmap_ptr( shared_ptr<XMap> & xmap_ptr )
{
	boost::any any = static_cast< shared_ptr<XMap> const>( xmap_ptr );
	if( get_xmap_ptr( &any ) ){
		xmap_ptr = boost::any_cast< shared_ptr<XMap> const>( any );
		return true;
	}
	return false;
}
template< class XMap >
bool RifBase::get_xmap_const_ptr( shared_ptr<XMap const> & xmap_ptr ) const
{
	boost::any any = static_cast< shared_ptr<XMap const> const>( xmap_ptr );
	if( get_xmap_const_ptr( &any) ){
		xmap_ptr = boost::any_cast< shared_ptr<XMap const> const>( any );
		return true;
	}
	return false;
}
template< class XMap >
bool RifBase::set_xmap_ptr( shared_ptr<XMap> const & xmap_ptr )
{
	boost::any any = static_cast< shared_ptr<XMap> >( xmap_ptr );
	return set_xmap_ptr( &any );
}


}}

#endif
