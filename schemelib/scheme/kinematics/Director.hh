#ifndef INCLUDED_kinematics_Director_HH
#define INCLUDED_kinematics_Director_HH


#include <scheme/kinematics/Director.meta.hh>

#include "scheme/types.hh"
#include "scheme/kinematics/Scene.hh"
#include "scheme/nest/MultiNest.hh"

#include <boost/any.hpp>

#include <vector>


/*
	NOTES

	the boost any stuff does have a performance impact...

	AS IS: set_scene rate: 1.26338e+06 / sec
	without new/del       ~1.75000e+06 / sec

	directly setting scene from NEST without boost::any or virtual func calls or multinest:
	set_scene rate: 1.86553e+07 / sec [       OK ] Director.basic_test (37 ms)


*/

namespace scheme {
namespace kinematics {


template<
	class _Position,
	class _BigIndex, /*= uint64_t, */ // These are defined in Director.meta.hh
	class _Index /* = uint64_t */
>
struct Director
{
	typedef _Position Position;
	typedef _BigIndex BigIndex;
	typedef _Index Index;
	typedef SceneBase<Position,Index> Scene;

	virtual
	bool
	set_scene(
		BigIndex const & i,
		int resl,
		Scene & scene
	) const = 0;

	virtual BigIndex size(int resl, BigIndex sizes) const = 0;

};

inline
void
set_nest_size(uint64_t set, uint64_t & sizes) {
	sizes = set;
}

template<class RifDockIndex>
void
set_nest_size(uint64_t set, RifDockIndex & sizes) {
	sizes.nest_index = set;
}

inline
uint64_t
get_nest_index(uint64_t index) {
	return index;
}

template<class RifDockIndex>
uint64_t
get_nest_index(RifDockIndex index) {
	return index.nest_index;
}


// This class is dumb, if used in a CompositeDirector, it must be first
template< class _Nest, class _BigIndex >
struct NestDirector
 :  public Director<
		typename _Nest::Value,
		_BigIndex
		// typename _Nest::Index // brian - This is related to nbodies. Should not be same as _Nest::Index
	> {
	typedef _Nest Nest;
	typedef typename Nest::Value Position;
	typedef _BigIndex BigIndex;
	typedef typename Nest::Index Index;
	typedef SceneBase<Position,Index> Scene;

	typedef typename Nest::Index NestIndex;

	Index ibody_;
	Nest nest_;

	NestDirector() : ibody_(0),nest_() {}
	NestDirector( Index ibody ) : ibody_(ibody),nest_() {}
	template<class A>
	NestDirector( A const & a, Index ibody ) : ibody_(ibody),nest_(a) {}
	template<class A, class B>
	NestDirector( A const & a, B const & b, Index ibody ) : ibody_(ibody),nest_(a,b) {}
	template<class A, class B, class C>
	NestDirector( A const & a, B const & b, C const & c, Index ibody ) : ibody_(ibody),nest_(a,b,c) {}
	template<class A, class B, class C, class D>
	NestDirector( A const & a, B const & b, C const & c, D const & d, Index ibody ) : ibody_(ibody),nest_(a,b,c,d) {}

	Nest const & nest() const { return nest_; }

	virtual
	bool
	set_scene(
		BigIndex const & i,
		int resl,
		Scene & scene
	) const {
		Position p;

		NestIndex ni = get_nest_index(i);
		bool success = nest_.get_state( ni, resl, p );
		if( !success ) return false;
		scene.set_position( ibody_, p );
		return true;
	}
        

	virtual BigIndex size(int resl, BigIndex sizes) const {
		set_nest_size(nest_.size(resl), sizes);
		return sizes;
	}


};

template< class Nest, class BigIndex >
std::ostream & operator << ( std::ostream & out, NestDirector<Nest,BigIndex> const & d ){
	out << "NestDirector " << d.nest();
	return out;
}




template<
	class _Position,
	class _RifDockIndex
>
struct CompositeDirector
 : 	public Director<
 		_Position,
 		_RifDockIndex
 	> {

 	typedef Director< _Position, _RifDockIndex> Base;
	typedef typename Base::Position Position;
	typedef typename Base::BigIndex BigIndex;
	typedef typename Base::Index Index;
	typedef typename Base::Scene Scene;

	std::vector<shared_ptr<Base>> directors_;

	CompositeDirector( std::vector<shared_ptr<Base>> director_list ) :
		directors_( director_list ) {}

	virtual
	bool
	set_scene(
		BigIndex const & i,
		int resl,
		Scene & scene
	) const override {

		for ( shared_ptr<Base> director : directors_ ) {
			bool success = director->set_scene( i, resl, scene );
			if ( !success ) return false;
		}

		return true;
	}

	virtual BigIndex size(int resl, BigIndex sizes) const override {
		for ( shared_ptr<Base> director : directors_ ) {
			sizes = director->size(resl, sizes);
		}
		return sizes;
	}



};
    

    // This must be used as part of a CompositeDirector
    template<
    class _Position,
    class _SeedingPositions,
    class _RifDockIndex
    >
    struct SeedingDirector
    :  public Director<
    _Position,
    _RifDockIndex
    > {
        typedef Director< _Position, _RifDockIndex> Base;
        typedef _SeedingPositions SeedingPositions;
        typedef shared_ptr<SeedingPositions> SeedingPositionsOP;
        typedef typename Base::Position Position;
        typedef typename Base::BigIndex BigIndex;
        typedef typename Base::Index Index;
        typedef typename Base::Scene Scene;
        
        
        Index ibody_;
        SeedingPositionsOP seeding_positions_;
        double maximum_allowed_rotation_ang_;
        
        SeedingDirector(SeedingPositionsOP sp) : ibody_(1),seeding_positions_(sp), maximum_allowed_rotation_ang_ ( -1.0 ) {}
        // I thinks this ho
        SeedingDirector(SeedingPositionsOP sp, Index ibody, double max_ang_degree ) : ibody_(ibody),seeding_positions_(sp), maximum_allowed_rotation_ang_ (max_ang_degree / 180 * 3.14)  {}
        
        
        virtual
        bool
        set_scene(
                  BigIndex const & i,
                  int resl,
                  Scene & scene
                  ) const override {
            runtime_assert( seeding_positions_ );
            
            uint64_t si = i.seeding_index;
            
            Position ori_pos = scene.position(1);
            
            if (  maximum_allowed_rotation_ang_ > 0 && std::abs (Eigen::AngleAxisf( ori_pos.rotation() ).angle() ) > maximum_allowed_rotation_ang_ ) {
                return false;
            }
            
            Position new_pos(Position::Identity());
            new_pos.rotate( ori_pos.rotation() * seeding_positions_->at(si).rotation() );
            
            new_pos.translation() = ori_pos.translation() + seeding_positions_->at(si).translation();
            
            
            scene.set_position(1, new_pos);
            
            
            return true;
        }
        
        // change the BigIndex Struct.
        virtual BigIndex size(int resl, BigIndex sizes) const override {
            sizes.seeding_index = seeding_positions_->size();
            return sizes;
        }
        
        
        
    };
    
    template< class Position, class SeedingPositions, class RifDockIndex >
    std::ostream & operator << ( std::ostream & out, SeedingDirector<Position,SeedingPositions,RifDockIndex> const & d ){
        out << "ScaffoldNestDirector ";
        return out;
    }



// This must be used as part of a CompositeDirector
template< 
	class _Position,
	class _ScaffoldProvider,
	class _RifDockIndex
>
struct ScaffoldDirector
 :  public Director<
		_Position,
		_RifDockIndex
	> {
	typedef Director< _Position, _RifDockIndex> Base;
	typedef _ScaffoldProvider ScaffoldProvider;
	typedef shared_ptr<ScaffoldProvider> ScaffoldProviderOP;
	typedef typename Base::Position Position;
	typedef typename Base::BigIndex BigIndex;
	typedef typename Base::Index Index;
	typedef typename Base::Scene Scene;

	typedef typename ScaffoldProvider::ScaffoldIndex ScaffoldIndex;

	Index ibody_;
	ScaffoldProviderOP scaffold_provider_;

	ScaffoldDirector(ScaffoldProviderOP sp) : ibody_(0),scaffold_provider_(sp) {}
	ScaffoldDirector(ScaffoldProviderOP sp, Index ibody ) : ibody_(ibody),scaffold_provider_(sp) {}


	virtual
	bool
	set_scene(
		BigIndex const & i,
		int resl,
		Scene & scene
	) const override {
		runtime_assert( scaffold_provider_ );

		ScaffoldIndex si = i.scaffold_index;

		scene.replace_body( ibody_, scaffold_provider_->get_scaffold( si ) );

		return true;
	}

	// This doesn't actually work, the format of ScaffoldProvider.size() can't be returned in this format
	virtual BigIndex size(int resl, BigIndex sizes) const override {
		sizes.scaffold_index = scaffold::scaffold_index_default_value(ScaffoldIndex());
		return sizes;
	}



};

template< class Position, class ScaffoldProvider, class RifDockIndex >
std::ostream & operator << ( std::ostream & out, ScaffoldDirector<Position,ScaffoldProvider,RifDockIndex> const & d ){
	out << "ScaffoldNestDirector ";
	return out;
}


inline
uint64_t
director_index_default_value() {
	return 0;
}



class SymElem {};

template<
	class _Position,
	class _Index = uint64_t
>
struct SceneTree {
	typedef SceneTree<_Position,_Index> THIS;
	typedef _Position Position;
	typedef _Index Index;
	typedef SceneBase<Position,Index> Scene;
	typedef nest::NestBase<Index> Nest;
	typedef shared_ptr<Nest> NestP;
	typedef shared_ptr<SceneTree> SceneTreeP;

	Position edge_xform_;
	std::vector<NestP> position_nests_;
	std::vector<NestP> sym_elem_nests_;
	std::vector<NestP> conformation_nests_;

	std::vector<Index> bodies_;
	std::vector<SymElem> sym_elems_;

	std::vector<SceneTreeP> children_;


	// mutable Position tmp_;

	SceneTree( Position const & edge_xform = Position::Identity() ) : edge_xform_( edge_xform ) {}

	void add_body(Index ibody) { bodies_.push_back(ibody); }
	void add_position_nest( NestP nestp ) { position_nests_.push_back(nestp); }
	void add_child( SceneTreeP child ) { children_.push_back(child); }

	void
	get_nests(
		std::vector<NestP> & nests
	) const {
		nests.insert( nests.end(),     position_nests_.begin(),     position_nests_.end() );
		nests.insert( nests.end(),     sym_elem_nests_.begin(),     sym_elem_nests_.end() );
		nests.insert( nests.end(), conformation_nests_.begin(), conformation_nests_.end() );
		BOOST_FOREACH( SceneTreeP child, children_ ){
			child->get_nests( nests );
		}
	}

	void
	get_dofs( std::vector<boost::any> & dofs ) const {
		for( int ipos = 0; ipos < position_nests_.size(); ++ipos){
			// std::cerr << "create temporary Position" << std::endl;
			dofs.push_back( new Position );
			// Position *tmp = &tmp_;
			// dofs.push_back( tmp );
		}
		if(     sym_elem_nests_.size() ) throw std::logic_error("sym_elem_nests not implemented");
		if( conformation_nests_.size() ) throw std::logic_error("conformation nests not implemented");
		BOOST_FOREACH( SceneTreeP child, children_ ){
			child->get_dofs( dofs );
		}

	}
	void
	clear_empty_dofs( std::vector<boost::any>::iterator idofs ) const {
		for( int ipos = 0; ipos < position_nests_.size(); ++ipos){
			// std::cerr << "delete temporary Position" << std::endl;
			delete boost::any_cast<Position*>(*idofs);
			++idofs;
		}
		if(     sym_elem_nests_.size() ) throw std::logic_error("sym_elem_nests not implemented");
		if( conformation_nests_.size() ) throw std::logic_error("conformation nests not implemented");
		BOOST_FOREACH( SceneTreeP child, children_ ){
			child->clear_empty_dofs( idofs );
		}

	}

	void
	set_scene(
		Position const & parent_position,
		std::vector<boost::any>::const_iterator idof,
		Scene & scene
	) const {
		Position position( edge_xform_ * parent_position );
		for( int ipos = 0; ipos < position_nests_.size(); ++ipos){
			// std::cerr << "set_scene get nest xform" << std::endl;
			Position & nest_xform( *boost::any_cast<Position*>(*idof) );
			// std::cerr << "set_scene get nest xform DONE" << std::endl;
			position = nest_xform * position;
			++idof;
		}
		if(     sym_elem_nests_.size() ) throw std::logic_error("sym_elem_nests not implemented");
		if( conformation_nests_.size() ) throw std::logic_error("conformation nests not implemented");
		BOOST_FOREACH( Index ibody, bodies_ ){
			scene.set_position(ibody,position);
		}
		BOOST_FOREACH( SceneTreeP child, children_ ){
			child->set_scene( parent_position, idof, scene );
		}
	}

};

	// MultiNest(Nests const & nests)
	// void add_nest( NestP nestp )

template<
	class _Position,
	class _BigIndex = uint64_t,
	class _Index = uint64_t
>
struct TreeDirector : public Director<_Position,_BigIndex,_Index> {

	typedef _BigIndex BigIndex;
	typedef _Index Index;
	typedef _Position Position;
	// typedef SceneTree<Position,Index> SceneTree;
	typedef shared_ptr<SceneTree<Position,Index> > SceneTreeP;
	typedef SceneBase<Position,Index> Scene;
	typedef nest::NestBase<Index> Nest;
	typedef shared_ptr<Nest> NestP;
	typedef nest::MultiNest<Index,BigIndex> MultiNest;

	SceneTreeP root_;
	MultiNest multinest_;

	TreeDirector( SceneTreeP root ) : root_(root) { init(); }

	void init(){
		std::vector<NestP> nests;
		root_->get_nests( nests );
		multinest_.init( nests );
	}

	virtual
	bool
	set_scene(
		BigIndex const & i,
		int resl,
		Scene & scene
	) const {
		std::vector<boost::any> tmp_dofs;
		// std::cerr << "get_dofs" << std::endl;
		root_->get_dofs( tmp_dofs );
		assert( tmp_dofs.size() == multinest_.size() );
		// std::cerr << "get_states" << std::endl;
		bool success = multinest_.get_states( i, resl, tmp_dofs );
		if( !success ) return false;
		// std::cerr << *boost::any_cast<numeric::X1dim*>(tmp_dofs.front()) << std::endl;
		assert( root_ != 0 );
		// std::cerr << "set_scene" << std::endl;
		root_->set_scene( Position::Identity(), tmp_dofs.begin(), scene );
		// std::cerr << "clear_empty_dofs" << std::endl;
		root_->clear_empty_dofs( tmp_dofs.begin() );
		return true;
	}

	virtual void size(int resl, BigIndex & sizes) const { 
		sizes = multinest_.size(resl);
		return sizes;
	}


private:
	TreeDirector() {}
};


}
}

#endif
