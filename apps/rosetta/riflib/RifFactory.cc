// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols

#include <riflib/RifFactory.hh>

#include <riflib/util.hh>
#include <riflib/rifdock_typedefs.hh>

#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/file/file_sys_util.hh>
#include <ObjexxFCL/format.hh>
#include <boost/format.hpp>

#include <riflib/rif/RifAccumulators.hh>


#include <scheme/objective/hash/XformMap.hh>
#include <scheme/objective/storage/RotamerScores.hh>

#include <scheme/actor/Atom.hh>
#include <scheme/actor/BackboneActor.hh>
#include <scheme/actor/VoxelActor.hh>
#include <scheme/kinematics/Scene.hh>
#include <scheme/nest/pmap/OriTransMap.hh>
#include <scheme/objective/ObjectiveFunction.hh>
#include <scheme/search/HackPack.hh>

#include <riflib/scaffold/ScaffoldDataCache.hh>
#include <complex>


namespace devel {
namespace scheme {


template<class C> void call_sort_rotamers( C & v ){ v.second.sort_rotamers(); }
template<class C> void assert_is_sorted  ( C & v ){ runtime_assert( v.second.is_sorted() ); }

template<class XmapIter>
struct XmapKeyIterHelper : public KeyIterHelperBase<typename RifBase::Key> {
    typedef typename RifBase::Key Key;
    XmapKeyIterHelper(XmapIter iter) : iter_(iter) {}
    Key get_key() const override { return iter_->first; }
    void next() override { ++iter_; }
    bool equal(KeyIterHelperBase const & that) const override {
        XmapKeyIterHelper const & concrete = static_cast<XmapKeyIterHelper const &>(that);
        return iter_ == concrete.iter_;
    }
    XmapIter iter_;
};

template< class XMap >
class RifWrapper : public RifBase {

	shared_ptr<XMap> xmap_ptr_;

public:

private:
	RifWrapper() : RifBase() {}
public:
	RifWrapper( shared_ptr<XMap> xmap_ptr, std::string type ) : RifBase(type), xmap_ptr_(xmap_ptr) {}

	virtual bool load( std::istream & in , std::string & description )
	{
		size_t s;
		in.read((char*)&s,sizeof(size_t));
		char buf[9999];
		for(int i = 0; i < 9999; ++i) buf[i] = 0;
		in.read(buf,s);
		std::string type_in(buf);
		runtime_assert_msg( type_in == type_, "mismatched rif_types, expected: '" + type_ + "' , got: '" + type_in + "'" );
		return xmap_ptr_->load( in , description );
	}
	virtual bool save( std::ostream & out, std::string & description ) {
		size_t s = type_.size();
		out.write((char*)&s,sizeof(size_t));
		out.write(type_.c_str(),s);
		return xmap_ptr_->save( out, description );
	}

	virtual bool get_xmap_ptr( boost::any * any_p )	{
		bool is_compatible_type =     boost::any_cast< shared_ptr<XMap> const>( any_p );
		if( is_compatible_type ) *any_p = static_cast< shared_ptr<XMap> const>( xmap_ptr_ );
		return is_compatible_type;
	}
	virtual bool get_xmap_const_ptr( boost::any * any_p ) const	{
		bool is_compatible_type =     boost::any_cast< shared_ptr<XMap const> const>( any_p );
		if( is_compatible_type ) *any_p = static_cast< shared_ptr<XMap const> const>( xmap_ptr_ );
		return is_compatible_type;
	}
	virtual bool set_xmap_ptr( boost::any * any_p )	{
		shared_ptr<XMap> const * tmp = boost::any_cast< shared_ptr<XMap> >( any_p );
		if( !tmp ) return false;
		xmap_ptr_ = *tmp;
		return true;
	}


    //Brian
    virtual std::pair< int, int > get_sat1_sat2( EigenXform const & x, int roti ) const override {
    	std::pair< int, int > sat1_sat2( -1, -1);

    	Key k = get_bin_key(x);


        auto const & rs = (*xmap_ptr_)[k];
        int const Nrots = XMap::Value::N;
        for( int i = 0; i < Nrots; ++i ){
            if( rs.empty(i) ) break;
            if (rs.rotamer(i) != roti) continue; 
            std::vector<int> sats;
        	rs.rotamer_sat_groups( i, sats );
        	if ( sats.size() == 0) {
        		return sat1_sat2;
        	}
        	sat1_sat2.first = sats[0];
        	if ( sats.size() == 1 ) {
        		return sat1_sat2;
        	}
        	sat1_sat2.second = sats[1];
        	return sat1_sat2;
        }

        return sat1_sat2;
    }

	size_t size() const override { return xmap_ptr_->size(); }
	float load_factor() const override { return xmap_ptr_->map_.size()*1.f/xmap_ptr_->map_.bucket_count(); }
	size_t mem_use()    const override { return xmap_ptr_->mem_use(); }
	float cart_resl()   const override { return xmap_ptr_->cart_resl_; }
	float ang_resl()    const override { return xmap_ptr_->ang_resl_; }
	size_t sizeof_value_type() const override { return sizeof(typename XMap::Map::value_type); }
	bool  has_sat_data_slots() const override { return XMap::Value::RotScore::UseSat; }

	// will resize to accomodate highest number rotamer
	void get_rotamer_ids_in_use( std::vector<bool> & using_rot ) const override
	{
		typedef typename XMap::Map::value_type MapPair;
		typedef typename MapPair::second_type RotScores;
		for( auto const & v : xmap_ptr_->map_ ){
			RotScores const & xmrot( v.second );
			for( int i = 0; i < RotScores::N; ++i ){
				if( xmrot.empty(i) ) break;
				if( xmrot.rotamer(i) >= using_rot.size() ) using_rot.resize( xmrot.rotamer(i)+1 , false );
				using_rot[ xmrot.rotamer(i) ] = true;
			}
		}

	}

    Key get_bin_key( EigenXform const & x) const override {
        return xmap_ptr_->get_key(x);
    }

    EigenXform get_bin_center( Key const & k) const override {
        return xmap_ptr_->get_center(k);
    }

    void get_rotamers_for_key( Key const & k, std::vector< std::pair< float, int > > & rotscores ) const override{
        typename XMap::Value const & rs = (*xmap_ptr_)[k];
        int const Nrots = XMap::Value::N;
        for( int i = 0; i < Nrots; ++i ){
            if( rs.empty(i) ) break;
            rotscores.push_back( std::make_pair<float,int>( rs.score(i), rs.rotamer(i) ) );
        }
    }

	void
	get_rotamers_for_xform(
		EigenXform const & x,
		std::vector< std::pair< float, int > > & rotscores
	) const override {
        Key const & k = xmap_ptr_->get_key(x);
        get_rotamers_for_key(k, rotscores);
	}

	void finalize_rif() override {
		// sort the rotamers in each cell so best scoring is first
		__gnu_parallel::for_each( xmap_ptr_->map_.begin(), xmap_ptr_->map_.end(), call_sort_rotamers<typename XMap::Map::value_type> );
	}

	// void super_print( std::ostream & out, shared_ptr< RotamerIndex > rot_index_p ) const override { xmap_ptr_->super_print( out, rot_index_p  ); }
	void print( std::ostream & out ) const override { out << (*xmap_ptr_) << std::endl; }
	std::string value_name() const override { return XMap::Value::name(); }

	void collision_analysis( std::ostream & out ) const override {
		using namespace ObjexxFCL::format;
		typedef typename XMap::Value XMapVal;
		int64_t rif_num_collisions  [ XMapVal::N ];
		double  rif_avg_scores      [ XMapVal::N ];
		int64_t rif_avg_scores_count[ XMapVal::N ];
		for( int i = 0; i < XMapVal::N; ++i ){ rif_num_collisions[i]=0; rif_avg_scores[i]=0; rif_avg_scores_count[i]=0; }
		for( auto const & v : xmap_ptr_->map_ ){
			for( int i = 0; i < XMapVal::N; ++i ){
				bool not_empty = !v.second.rotscores_[i].empty();
				if( not_empty ){
					rif_num_collisions[i] += 1;
					rif_avg_scores[i] += v.second.rotscores_[i].score();
					// std::out << v.second.rotscores_[i].score() << std::endl; // WHY SOME WAY TOO LOW?????? fixed.
					rif_avg_scores_count[i]++;
				}
			}
		}
		for( int i = 0; i < XMapVal::N; ++i ) rif_avg_scores[i] /= rif_avg_scores_count[i];

		// out << "======================================================================" << std::endl;
		out << "======================= frac collisions in map =======================" << std::endl;
		out << "======================= cart: "<< xmap_ptr_->cart_resl_ <<
						", ori: " <<  xmap_ptr_->ang_resl_ <<
						" =======================" << std::endl;
		out << "======================================================================" << std::endl;
		float Ecollision = 0.0;
		for( int i = 0; i < XMapVal::N; ++i ){
			float colfrac = rif_num_collisions[i]*1.0/xmap_ptr_->map_.size();
			out << "   Nrots " << I(3,i+1) << " " << F(7,5,colfrac) << " " << F(7,3,rif_avg_scores[i]) << " " << rif_avg_scores_count[i] << std::endl;
			if( i > 0 ){
				float pcolfrac = rif_num_collisions[i-1]*1.0/xmap_ptr_->map_.size();
				Ecollision += i * (pcolfrac-colfrac);
			}
		}
		Ecollision += rif_num_collisions[XMapVal::N-1]*1.0/xmap_ptr_->map_.size() * XMapVal::N;
		out << "E(collisions) = " << Ecollision << std::endl;
		out << "======================================================================" << std::endl;

	}

    RifBaseKeyRange key_range() const override {
        auto b = std::make_shared<XmapKeyIterHelper<typename XMap::Map::const_iterator>>(
            ((typename XMap::Map const &)xmap_ptr_->map_).begin()  );
        auto e = std::make_shared<XmapKeyIterHelper<typename XMap::Map::const_iterator>>(
            ((typename XMap::Map const &)xmap_ptr_->map_).end()  );
        return RifBaseKeyRange(RifBaseKeyIter(b), RifBaseKeyIter(e));
    }


    // This looks for rifgen rotamers that have their N, CA, CB, and last atom within dump_dist of the residue given
    bool dump_rotamers_near_res( core::conformation::Residue const & res, std::string const & file_name, 
                                        float dump_dist, float dump_frac, shared_ptr<RotamerIndex> rot_index_p ) const override {

        std::string name3 = res.name3();


        std::pair<int, std::string> cb = ( name3 == "GLY" ? std::pair<int, std::string> {1, "CA"} : std::pair<int, std::string>{3, "CB"});
        // rif-number, rosetta number
        std::vector<std::pair<int, std::string>> align_pairs { 
            {0, "N"},
            {1, "CA"},
            cb
        };

        // prepare listed atoms
        std::vector<Eigen::Vector3f> scaff_atoms;
        for ( int i = 0; i < align_pairs.size(); i++ ) {
            std::pair<int, std::string> pair = align_pairs[i];
            numeric::xyzVector<core::Real> xyz = res.xyz(pair.second);
            Eigen::Vector3f vec;
            vec[0] = xyz.x();
            vec[1] = xyz.y();
            vec[2] = xyz.z();
            scaff_atoms.push_back(vec);
        }
        // prepare last atom
        numeric::xyzVector<core::Real> xyz = res.xyz(res.natoms());
        Eigen::Vector3f vec;
        vec[0] = xyz.x();
        vec[1] = xyz.y();
        vec[2] = xyz.z();
        scaff_atoms.push_back(vec);

    	std::cout << "Looking for rotamers within " << dump_dist << "A of in input pdb and of aa " << name3 << std::endl;
    	const RifBase * base = this;
		shared_ptr<XMap const> from;
		base->get_xmap_const_ptr( from );


		float coarse_dist_sq = (dump_dist + 5) * (dump_dist + 5);
		float dump_dist_sq = dump_dist * dump_dist;

		// transform and irot
		std::vector<std::pair<EigenXform, std::pair<int, float>>> to_dump;
		to_dump.reserve( from->map_.size() );


		std::pair<int,int> index_bounds = rot_index_p->index_bounds( name3 );


		for( auto const & v : from->map_ ){
			EigenXform x = from->hasher_.get_center( v.first );

			float dist_sq = (x.translation() - scaff_atoms[2]).squaredNorm();

            if (dist_sq > coarse_dist_sq) continue;

            typename XMap::Value const & rotscores = from->operator[]( x );
            static int const Nrots = XMap::Value::N;
            for( int i_rs = 0; i_rs < Nrots; ++i_rs ){
                if( rotscores.empty(i_rs) ) {
                    break;
                }

                int irot = rotscores.rotamer(i_rs);
                if (irot < index_bounds.first || irot >= index_bounds.second ) continue;

                // check listed atoms
                bool all_good = true;
                for ( int i = 0; i < align_pairs.size(); i++ ) {
                    int rif_atom_no = align_pairs[i].first;
                    Eigen::Vector3f scaff_vec = scaff_atoms[i];


                    SchemeAtom atom = rot_index_p->rotamers_.at( irot ).atoms_[rif_atom_no];
                    Eigen::Vector3f atom_vec = x * atom.position();
                    dist_sq = (atom_vec - scaff_vec).squaredNorm();

                    if (dist_sq > dump_dist_sq) {
                        all_good = false;
                        break;
                    }
                }

                if ( ! all_good ) continue;

                // check last atom
                Eigen::Vector3f scaff_vec = scaff_atoms.back();
                SchemeAtom atom = rot_index_p->rotamers_.at( irot ).atoms_.back();
                Eigen::Vector3f atom_vec = x * atom.position();
                dist_sq = (atom_vec - scaff_vec).squaredNorm();

                if (dist_sq > dump_dist_sq) continue;
                







				float score = rotscores.score(i_rs);

				to_dump.push_back(std::pair<EigenXform, std::pair<int, float>>(x, std::pair<int, float>(irot, score)));
			
			
			}
		}


        if (to_dump.size() == 0) {
            std::cout << "No rotamers found!!!!" << std::endl;
            return false;
        }

		uint64_t num_dump = to_dump.size() * dump_frac;
		uint64_t dump_every = to_dump.size() / num_dump;


		std::cout << "Found " << to_dump.size() << " rotamers. Dumping " << num_dump << " to " << file_name << " ..." << std::endl;


		utility::io::ozstream out( file_name );
		uint64_t dumped = 0;
		for ( uint64_t i = 0; i < to_dump.size(); i ++ ) {
			if ( i % dump_every != 0 ) {
				continue;
			}
			EigenXform x = to_dump[i].first;
			auto inner_pair = to_dump[i].second;
			int irot = inner_pair.first;
			float score = inner_pair.second;


			out << std::string("MODEL") << " " << boost::str(boost::format("%.3f")%score) << std::endl;

            BOOST_FOREACH( SchemeAtom a, rot_index_p->rotamers_.at( irot ).atoms_ ){
                a.set_position( x * a.position() ); 
                a.nonconst_data().resnum = dumped;
                a.nonconst_data().chain = 'A';
                ::scheme::actor::write_pdb( out, a, nullptr );
            }

            dumped ++;

			out << std::string("ENDMDL") << std::endl;

		}

		out.close();


        return true;

    }




};







std::string get_rif_type_from_file( std::string fname )
{
	runtime_assert( utility::file::file_exists(fname) );
	utility::io::izstream in( fname );
	runtime_assert( in.good() );
	size_t s;
	in.read((char*)&s,sizeof(size_t));
	char buf[9999];
	for(int i = 0; i < 9999; ++i) buf[i] = 0;
	in.read(buf,s);
	return std::string(buf);
}





//////////////////////////////////// ScoreBBActorVsRIF ////////////////////////////////////////////

	typedef int32_t intRot;

	// should move to libraries somewhere
	// empty class, serves only to position RIF in absolute space
	// struct RIFAnchor {
	// 	RIFAnchor() {}
	// };
	std::ostream & operator<<( std::ostream & out, RIFAnchor const& va ){
		return out << "RIFAnchor";
	}
	template< class MetaData >
	void write_pdb( std::ostream & , RIFAnchor const &, MetaData const & ){}

	struct ScoreBBActorvsRIFResult {
		float val_;
		std::vector< std::pair<intRot,intRot> > rotamers_;
		ScoreBBActorvsRIFResult() : val_(0) {}
		ScoreBBActorvsRIFResult( float f ) : val_(f) {}
		operator float() const { return val_; }
		void operator=( float f ) { val_ = f; }
		void operator+=( float f ) { val_ += f; }
		bool operator<( ScoreBBActorvsRIFResult const & other ) const { return val_ < other.val_; }
	};
	struct ScoreBBActorvsRIFScratch {
		shared_ptr< ::scheme::search::HackPack> hackpack_;
		std::vector<bool> is_satisfied_;
		std::vector<bool> has_rifrot_;
        std::vector<std::vector<float> > const * rotamer_energies_1b_ = nullptr;
        std::vector< std::pair<int,int> > const * scaffold_rotamers_ = nullptr;
		// sat group vector goes here
		//std::vector<float> is_satisfied_score_;
	};

	template< class BBActor, class RIF, class VoxelArrayPtr >
	struct ScoreBBActorVsRIF
	{
		typedef boost::mpl::true_ HasPre;
		typedef boost::mpl::true_ HasPost;
		typedef ScoreBBActorvsRIFScratch Scratch;
		typedef ScoreBBActorvsRIFResult Result;
		typedef std::pair<RIFAnchor,BBActor> Interaction;
		bool packing_ = false;
		::scheme::search::HackPackOpts packopts_;
		int n_sat_groups_ = 0, require_satisfaction_ = 0, require_n_rifres_ = 0;
		std::vector< shared_ptr< ::scheme::search::HackPack> > packperthread_;
	private:
		shared_ptr<RIF const> rif_ = nullptr;
	public:
		VoxelArrayPtr target_proximity_test_grid_ = nullptr;
		devel::scheme::ScoreRotamerVsTarget<
				VoxelArrayPtr, ::scheme::chemical::HBondRay, ::devel::scheme::RotamerIndex
			> rot_tgt_scorer_;
		std::vector<int> always_available_rotamers_;

		ScoreBBActorVsRIF() {}

		 void clear() {
		 	packperthread_.clear();
		 }

		static std::string name(){ return "ScoreBBActorVsRIF"; }

		void set_rif( shared_ptr< ::devel::scheme::RifBase const> rif_ptr ){
			rif_ptr->get_xmap_const_ptr( rif_ );
		}

		void init_for_packing(
			// ::scheme::objective::storage::TwoBodyTable<float> const & twob,
			shared_ptr< ::devel::scheme::RotamerIndex > rot_index_p,
			std::vector<VoxelArrayPtr> const & target_field_by_atype,
			std::vector< ::scheme::chemical::HBondRay > const & target_donors,
			std::vector< ::scheme::chemical::HBondRay > const & target_acceptors,
			::scheme::search::HackPackOpts const & hackpackopts
		){
			packing_ = true;
			packperthread_.clear();
			for( int i  = 0; i < ::devel::scheme::omp_max_threads_1(); ++i ){
				shared_ptr< ::scheme::search::HackPack> tmp = make_shared< ::scheme::search::HackPack>(hackpackopts,rot_index_p->ala_rot(),i);
				packperthread_.push_back( tmp );
			}
			rot_tgt_scorer_.rot_index_p_ = rot_index_p;
			rot_tgt_scorer_.target_field_by_atype_ = target_field_by_atype;
			rot_tgt_scorer_.target_donors_ = target_donors;
			rot_tgt_scorer_.target_acceptors_ = target_acceptors;
			rot_tgt_scorer_.hbond_weight_ = hackpackopts.hbond_weight;
			rot_tgt_scorer_.upweight_iface_ = hackpackopts.upweight_iface;
			rot_tgt_scorer_.upweight_multi_hbond_ = hackpackopts.upweight_multi_hbond;
			packopts_ = hackpackopts;
			always_available_rotamers_.clear();
			for( int irot = 0; irot < rot_index_p->n_primary_rotamers(); ++irot ){
				switch( packopts_.always_available_rotamers_level ){
					case 2:
						always_available_rotamers_.push_back(irot);
						break;
					case 1:
						if( rot_index_p->resname(irot) == "VAL" ) always_available_rotamers_.push_back(irot);
						if( rot_index_p->resname(irot) == "ILE" ) always_available_rotamers_.push_back(irot);
						if( rot_index_p->resname(irot) == "LEU" ) always_available_rotamers_.push_back(irot);
						if( rot_index_p->resname(irot) == "MET" ) always_available_rotamers_.push_back(irot);
						if( rot_index_p->resname(irot) == "PHE" ) always_available_rotamers_.push_back(irot);
						break;
					case 0:
						break;
					default:
						utility_exit_with_message("unknown always_available_rotamers level");
				}
			}
		}

		template<class Scene, class Config>
		void pre( Scene const & scene, Result & result, Scratch & scratch, Config const & config ) const
		{

			// Added by brian ////////////////////////
			ScaffoldDataCacheOP data_cache = scene.conformation_ptr(1)->cache_data_;
            runtime_assert( data_cache );
			scratch.rotamer_energies_1b_ = data_cache->local_onebody_p.get();
			//////////////////////////////////////////

			runtime_assert( rif_ );
			runtime_assert( scratch.rotamer_energies_1b_ );
			if( n_sat_groups_ > 0 ){
				scratch.is_satisfied_.resize(n_sat_groups_,false); // = new bool[n_sat_groups_];
				for( int i = 0; i < n_sat_groups_; ++i ) scratch.is_satisfied_[i] = false;
				//scratch.is_satisfied_score_.resize(n_sat_groups_,0.0);
				//for( int i = 0; i < n_sat_groups_; ++i ) scratch.is_satisfied_score_[i] = 0;
			}
			scratch.has_rifrot_.resize(scratch.rotamer_energies_1b_->size(), false);
			for ( int i = 0; i < scratch.has_rifrot_.size(); i++ ) scratch.has_rifrot_[i] = false;
			if( !packing_ ) return;

			// Added by brian ////////////////////////
			scratch.scaffold_rotamers_ = data_cache->local_rotamers_p.get();
			//////////////////////////////////////////

			runtime_assert( scratch.scaffold_rotamers_ );
			runtime_assert( rot_tgt_scorer_.rot_index_p_ );
			runtime_assert( rot_tgt_scorer_.target_field_by_atype_.size() == 22 );
			scratch.hackpack_ = packperthread_.at( ::devel::scheme::omp_thread_num() );
			scratch.hackpack_->reinitialize( data_cache->local_twobody_p );
		}

		template<class Config>
		Result operator()( RIFAnchor const &, BBActor const & bb, Scratch & scratch, Config const& c ) const
		{
			if( target_proximity_test_grid_ && target_proximity_test_grid_->at( bb.position().translation() ) == 0.0 ){
				return 0.0;
			}

			typename RIF::Value const & rotscores = rif_->operator[]( bb.position() );
			static int const Nrots = RIF::Value::N;
			int const ires = bb.index_;
			float bestsc = 0.0;
			for( int i_rs = 0; i_rs < Nrots; ++i_rs ){
				if( rotscores.empty(i_rs) ) {
					break;
				}
				int irot = rotscores.rotamer(i_rs);
				float const rot1be = (*scratch.rotamer_energies_1b_).at(ires).at(irot);
				float score_rot_v_target = rotscores.score(i_rs);

				bool rotamer_satisfies = rotscores.do_i_satisfy_anything(i_rs);


				if( packing_ && packopts_.packing_use_rif_rotamers ){

					if( rot1be <= packopts_.rotamer_onebody_inclusion_threshold || rotamer_satisfies){
						
						float const recalc_rot_v_tgt = rot_tgt_scorer_.score_rotamer_v_target( irot, bb.position(), 10.0, 4 );
						score_rot_v_target = recalc_rot_v_tgt;
				
						if (( score_rot_v_target + rot1be < packopts_.rotamer_inclusion_threshold &&
						      score_rot_v_target          < packopts_.rotamer_inclusion_threshold ) || rotamer_satisfies){

							float sat_bonus = 0;
							if (rotamer_satisfies) {
								sat_bonus = packopts_.user_rotamer_bonus_per_chi * rot_tgt_scorer_.rot_index_p_->nchi(irot) +
								            packopts_.user_rotamer_bonus_constant;
								// std::cout << "ires " << ires << " cdirot " << irot << std::endl;
								// std::cout << "Sat bonus: " << sat_bonus << " Score: " << score_rot_v_target + rot1be << std::endl;
							}
							scratch.hackpack_->add_tmp_rot( ires, irot, score_rot_v_target + rot1be + sat_bonus );
						}
					}
					if( packopts_.use_extra_rotamers ){
						auto child_rots = rot_tgt_scorer_.rot_index_p_->child_map_.at(irot);
						for( int crot = child_rots.first; crot < child_rots.second; ++crot ){
							float const crot1be = (*scratch.rotamer_energies_1b_).at(ires).at(crot);
							if( crot1be > packopts_.rotamer_onebody_inclusion_threshold ) continue;
							float const recalc_crot_v_tgt = rot_tgt_scorer_.score_rotamer_v_target( crot, bb.position(), 10.0, 4 );
							if( recalc_crot_v_tgt + rot1be < packopts_.rotamer_inclusion_threshold &&
							    recalc_crot_v_tgt          < packopts_.rotamer_inclusion_threshold ){
								scratch.hackpack_->add_tmp_rot( ires, crot, recalc_crot_v_tgt + crot1be );
							}
						}
					}
				}

				float const score_rot_tot = score_rot_v_target + rot1be;
				//TaYi change the score_rot_tot cutoff for mark_sat_groups to include not so good rotamers 
				if( n_sat_groups_ > 0 && score_rot_tot < 5.0 ){
					rotscores.mark_sat_groups( i_rs, scratch.is_satisfied_ );
				}
                if ( score_rot_tot < 0.0 ) {
                    // std::cout << "Adding " << ires << std::endl;
                    scratch.has_rifrot_[ires] = true;
                }

				bestsc = std::min( score_rot_tot , bestsc );
				//}
			}
			// // add native scaffold rotamers TODO: this is bugged somehow?
			if( packing_ && packopts_.add_native_scaffold_rots_when_packing ){
				for( int irot = scratch.scaffold_rotamers_->at(ires).first; irot < scratch.scaffold_rotamers_->at(ires).second; ++irot ){
					if( scratch.hackpack_->using_rotamer( ires, irot ) ){
						float const rot1be = (*scratch.rotamer_energies_1b_).at(ires).at(irot);
						float const recalc_rot_v_tgt = rot_tgt_scorer_.score_rotamer_v_target( irot, bb.position(), 10.0, 4 );
						float const rot_tot_1b = recalc_rot_v_tgt + rot1be;
						if( rot_tot_1b < -1.0 && recalc_rot_v_tgt < -0.1 ){ // what's this logic??
							scratch.hackpack_->add_tmp_rot( ires, irot, rot_tot_1b );
						}
					}
				}
			}
			if( packing_ ){
				for( int irot : always_available_rotamers_ ){
					float const irot1be = (*scratch.rotamer_energies_1b_).at(ires).at(irot);
					if( irot1be > packopts_.rotamer_onebody_inclusion_threshold ) continue;
					int sat1=-1, sat2=-1;
					float const recalc_rot_v_tgt = rot_tgt_scorer_.score_rotamer_v_target( irot, bb.position(), 10.0, 4 );
					if( recalc_rot_v_tgt + irot1be < packopts_.rotamer_inclusion_threshold &&
					    recalc_rot_v_tgt           < packopts_.rotamer_inclusion_threshold ){
						scratch.hackpack_->add_tmp_rot( ires, irot, recalc_rot_v_tgt + irot1be );
					}
				}
			}

			bestsc = std::min( bestsc, rotscores.score(Nrots-1) ); // in case not all rots stored...
			return bestsc;
		}

		template<class Scene, class Config>
		void post( Scene const & scene, Result & result, Scratch & scratch, Config const & config ) const
		{
			if( packing_ ){

				::scheme::search::HackPack & packer( *scratch.hackpack_ );
				result.val_ = packer.pack( result.rotamers_ );
				if( n_sat_groups_ > 0 ) for( int i = 0; i < n_sat_groups_; ++i ) scratch.is_satisfied_[i] = false;
				//std::vector< std::pair<intRot,intRot> > selected_rotamers;
				for( int i = 0; i < result.rotamers_.size(); ++i ){
					BBActor const & bb = scene.template get_actor<BBActor>( 1, result.rotamers_[i].first );
					int sat1=-1, sat2=-1;
					float const recalc_rot_v_tgt = rot_tgt_scorer_.score_rotamer_v_target_sat(
									result.rotamers_[i].second, bb.position(), sat1, sat2, 10.0, 4 );
					// todo: should do extra selection here?
					// if( recalc_rot_v_tgt < -1.0 ){
					// 	selected_rotamers.push_back( result.rotamers_[i] );
					// }
					if( n_sat_groups_ > 0 ){
						if( sat1 >= 0 ) scratch.is_satisfied_[ sat1 ] = true;
						if( sat2 >= 0 ) scratch.is_satisfied_[ sat2 ] = true;
					}
				}
				// result.rotamers_ = selected_rotamers;

			}

			if( n_sat_groups_ > 0 && !packing_ ){
				int nsat = 0;
				
				for( int i = 0; i < n_sat_groups_; ++i ){
					
					nsat += scratch.is_satisfied_[i];
					//result.val_ += scratch.is_satisfied_score_[i];
				}
				// if (nsat >= 4 ){
				// 	#pragma omp critical
				// 	{
				// 	std::cout << config << "     ";
					
				// 	for (int i = 0; i < 10; ++i){
				// 		std::cout << " "<< scratch.is_satisfied_[i];

				// 	}
				// 	std::cout << " " << std::endl;
				// 	}
				// }
				// std::cout << "here: " << nsat << std::endl;
				// runtime_assert( 0 <= nsat && nsat <= n_sat_groups_ );

				if( nsat - require_satisfaction_ < 0 ){
					result.val_ = 9e9;
				} //else {
					//result.val_ += -4.0f * nsat;
				//}

				// delete scratch.is_satisfied_;
				// scratch.is_satisfied_.clear();
			}

			// this is yolo code by Brian. I have no idea if th is will always work
			if ( require_n_rifres_ > 0 ) {
				if ( packing_ ) {
					std::map<int, bool> used_positions;
					for( int i = 0; i < result.rotamers_.size(); ++i ){
						int position = result.rotamers_[i].first;
						if ( used_positions.count(position) == 0 ) {
							used_positions[position] = true;
						}
					}
					// std::cout << result.rotamers_.size() << " ";
					if (used_positions.size() < require_n_rifres_ ) {
						result.val_ = 9e9;
					}
				} else {
					int count = 0;
					for ( int i = 0; i < scratch.has_rifrot_.size(); i++ ) {
						if ( scratch.has_rifrot_[i] ) {
							count ++;
						}
					}
					// std::cout << "Found " << count << std::endl;
					if (count < require_n_rifres_ ) {
						result.val_ = 9e9;
					}
				}
			}

				// #ifdef USE_OPENMP
				// #pragma omp critical
				// #endif
				// std::cout << result.rotamers_.size() << " " << result.val_ << std::endl;
			// if( result.rotamers_.size() == 0 ){
			// 	#ifdef USE_OPENMP
			// 	#pragma omp critical
			// 	#endif
			// 	std::cout << "no rotamers!" << std::endl;
			// }

			// std::cout << "bound score: " << result.val_ << ", packscore: " << packscore << std::endl;
			// for( int i = 0; i < result.rotamers_.size(); ++i ){
			// 	std::cout << "res: " << result.rotamers_[i].first << ", rotamer: " << result.rotamers_[i].second << std::endl;
			// }

		}
	};
	template< class B, class X, class V >
	std::ostream & operator<<( std::ostream & out, ScoreBBActorVsRIF<B,X,V> const& si ){ return out << si.name(); }


template< class XMap >
struct RifFactoryImpl :
	public RifFactory,
	enable_shared_from_this< RifFactoryImpl<XMap> >
 {


	typedef ScoreBBActorVsRIF<
			BBActor,
			// debugXMap,
			XMap, // no worky worky
			VoxelArrayPtr
		> MyScoreBBActorRIF;

	typedef ::scheme::objective::ObjectiveFunction<
			boost::mpl::vector<
				MyScoreBBActorRIF,
				MyClashScore
			>,
			int // Config type, just resl
		> MyRIFObjective;

	typedef ::scheme::objective::integration::SceneObjectiveParametric<
			ParametricScene,
			MyRIFObjective,
			MyScoreBBActorRIF
		> MySceneObjectiveRIF;

	RifFactoryImpl( RifFactoryConfig const & config ) : RifFactory(config) {}

	virtual RifPtr
	create_rif( float cart_resl=0, float ang_resl=0, float cart_bound=0 ) const
	{
	   if( cart_resl != 0 && ang_resl != 0 && cart_bound != 0 ){
	       return make_shared<RifWrapper<XMap> >( make_shared<XMap>( cart_resl, ang_resl, cart_bound ), this->config().rif_type );
	   } else if( cart_resl == 0 && ang_resl == 0 && cart_bound == 0 ){
	       return make_shared<RifWrapper<XMap> >( make_shared<XMap>(), this->config().rif_type );
	   } else {
	   		utility_exit_with_message("some XformMap constructor values specified, others not!");
	   }
	}

	virtual RifPtr
	create_rif_from_rif( RifConstPtr refrif, float cart_resl, float ang_resl, float cart_bound ) const {
		runtime_assert( this->config().rif_type == refrif->type() );
		RifPtr rif = create_rif( cart_resl, ang_resl, cart_bound );

		shared_ptr<XMap> to;
		rif->get_xmap_ptr( to );
		shared_ptr<XMap const> from;
		refrif->get_xmap_const_ptr( from );

		// std::cout << "create rif progress "; std::cout.flush();


		// old
		int progress0 = 0;
		for( auto const & v : from->map_ ){
			// if( ++progress0 % std::max((size_t)1,(from->size()/100)) == 0 ){
				// std::cout << '*'; std::cout.flush();
			// }
			EigenXform x = from->hasher_.get_center( v.first );

			uint64_t k = to->hasher_.get_key(x);
			typename XMap::Map::iterator iter = to->map_.find(k);
			if( iter == to->map_.end() ){
				to->map_.insert( std::make_pair(k,v.second) );
			} else {
				iter->second.merge( v.second );

			}
		}
		// // std::cout << std::endl;

		// new
		// int progress0 = 0;
		// for( auto const & v : from->map_ ){
		// 	// if( ++progress0 % std::max((size_t)1,(from->size()/100)) == 0 ){
		// 		// std::cout << '*'; std::cout.flush();
		// 	// }
		// 	EigenXform x = from->hasher_.get_center( v.first );
		// 	std::vector<uint64_t> keys = to->hasher_.get_key_and_nbrs(x);
		// 	for ( uint64_t const & k : keys ) {
		// 		typename XMap::Map::iterator iter = to->map_.find(k);
		// 		if( iter == to->map_.end() ){
		// 			to->map_.insert( std::make_pair(k,v.second) );
		// 		} else {
		// 			iter->second.merge( v.second );
		// 		}
		// 	}
		// }
		// std::cout << std::endl;


		return rif;
	}

	virtual	shared_ptr<rif::RifAccumulator>
	create_rif_accumulator( float cart_resl, float ang_resl, float cart_bound, size_t scratchM ) const {
		return make_shared< rif::RIFAccumulatorMapThreaded<XMap> >(
			this->shared_from_this(),
			cart_resl, ang_resl, cart_bound,
			scratchM
		);
	}

	virtual RifPtr
	create_rif_from_file( std::string const & fname, std::string & description ) const
	{
		RifPtr rif = this->create_rif();
		runtime_assert( rif );
		if( ! utility::file::file_exists(fname) ){
			utility_exit_with_message("create_rif_from_file missing file: " + fname );
		}
		utility::io::izstream in( fname );
		if( !in.good() ) return nullptr;
		bool success = rif->load( in, description );
		in.close();
		if( success ) return rif;
		else return nullptr;
	}

	virtual ScenePtr
	create_scene() const {
		ScenePtr s = make_shared<ParametricScene>(2);
		s->add_actor( 0, RIFAnchor() );
		return s;
	}

	virtual bool
	create_objectives(
		RifSceneObjectiveConfig const & config,
		std::vector<ObjectivePtr> & objectives,
		ObjectivePtr & packing_objective
	) const {

		for( int i_so = 0; i_so < config.rif_ptrs.size(); ++i_so ){
			if( i_so <= config.rif_ptrs.size()-1 ){
				shared_ptr< MySceneObjectiveRIF> objective = make_shared<MySceneObjectiveRIF>();
				objective->objective.template get_objective< MyScoreBBActorRIF >().set_rif( config.rif_ptrs[i_so] );
				if( i_so > 0 ){
					int i_tptg = std::min( 2, i_so-1 );
					// std::cout << "resl " << config.resolutions[i_so] << " using target_prox_grid " << config.resolutions[i_tptg] << std::endl;
					// use vdw grids as tgt proximity measure, use CH3 atom
					objective->objective.template get_objective<MyScoreBBActorRIF>().target_proximity_test_grid_ =
						config.target_bounding_by_atype->at(i_tptg).at(5);
				}
				objective->config = i_so;
				objectives.push_back( objective );
			}
		}
		// shared_ptr< MySceneObjectiveRIF> objective = make_shared<MySceneObjectiveRIF>();
		// objective->objective.template get_objective<MyScoreBBActorRIF>().set_rif( config.rif_ptrs.back() );
		// objective->config = config.rif_ptrs.size()-1;
		// objectives.push_back( objective );

		packing_objective = make_shared<MySceneObjectiveRIF>();
		dynamic_cast<MySceneObjectiveRIF&>(*packing_objective).objective.template get_objective<MyScoreBBActorRIF>().set_rif( config.rif_ptrs.back() );
		dynamic_cast<MySceneObjectiveRIF&>(*packing_objective).config = config.rif_ptrs.size()-1;
		if( config.require_satisfaction > 0 ){
			dynamic_cast<MySceneObjectiveRIF&>(*packing_objective).objective.template
				get_objective<MyScoreBBActorRIF>().n_sat_groups_ = config.n_sat_groups;
			dynamic_cast<MySceneObjectiveRIF&>(*packing_objective).objective.template
				get_objective<MyScoreBBActorRIF>().require_satisfaction_ = config.require_satisfaction;
		}


		for( auto op : objectives ){
			// dynamic_cast<MySceneObjectiveRIF&>(*op).objective.template get_objective<MyScoreBBActorRIF>().rotamer_energies_1b_ = config.local_onebody;
			// dynamic_cast<MySceneObjectiveRIF&>(*op).objective.template get_objective<MyScoreBBActorRIF>().scaffold_rotamers_ = config.local_rotamers;
			if( config.require_satisfaction > 0 ){
				dynamic_cast<MySceneObjectiveRIF&>(*op).objective.template get_objective<MyScoreBBActorRIF>().n_sat_groups_ = config.n_sat_groups;
				dynamic_cast<MySceneObjectiveRIF&>(*op).objective.template get_objective<MyScoreBBActorRIF>().require_satisfaction_ = config.require_satisfaction;
			}
			if (config.require_n_rifres > 0 ) {
				dynamic_cast<MySceneObjectiveRIF&>(*op).objective.template get_objective<MyScoreBBActorRIF>().require_n_rifres_ = config.require_n_rifres;
			}
		}
		// dynamic_cast<MySceneObjectiveRIF&>(*packing_objective).objective.template get_objective<MyScoreBBActorRIF>().rotamer_energies_1b_ = config.local_onebody;
		// dynamic_cast<MySceneObjectiveRIF&>(*packing_objective).objective.template get_objective<MyScoreBBActorRIF>().scaffold_rotamers_ = config.local_rotamers;

		// use 4.0A vdw grid for CH3 atoms as proximity test
		dynamic_cast<MySceneObjectiveRIF&>(*packing_objective).objective.template get_objective<MyScoreBBActorRIF>()
							.target_proximity_test_grid_ = config.target_bounding_by_atype->at(2).at(5);
		dynamic_cast<MySceneObjectiveRIF&>(*packing_objective).objective.template get_objective<MyScoreBBActorRIF>().init_for_packing(
			// *config.local_twobody,
			config.rot_index_p,
			*config.target_field_by_atype,
			*config.target_donors,
			*config.target_acceptors,
			*config.packopts
		);

		return true;

	}


};


shared_ptr<RifFactory>
create_rif_factory( RifFactoryConfig const & config )
{
	if( config.rif_type == "RotScore" )
	{
		typedef ::scheme::objective::storage::RotamerScore<> crfRotScore;
		typedef ::scheme::objective::storage::RotamerScores< 12, crfRotScore > crfXMapValue;
		BOOST_STATIC_ASSERT( sizeof( crfXMapValue ) == 24 );
		typedef ::scheme::objective::hash::XformMap<
				EigenXform,
				crfXMapValue,
				::scheme::objective::hash::XformHash_bt24_BCC6
			> crfXMap;
		BOOST_STATIC_ASSERT( sizeof( crfXMap::Map::value_type ) == 32 );

		return make_shared< RifFactoryImpl<crfXMap> >( config );
	}
	else if( config.rif_type == "RotScore64" )
	{
		typedef ::scheme::objective::storage::RotamerScore<> crfRotScore;
		typedef ::scheme::objective::storage::RotamerScores< 28, crfRotScore > crfXMapValue;
		BOOST_STATIC_ASSERT( sizeof( crfXMapValue ) == 56 );
		typedef ::scheme::objective::hash::XformMap<
				EigenXform,
				crfXMapValue,
				::scheme::objective::hash::XformHash_bt24_BCC6
			> crfXMap;
		BOOST_STATIC_ASSERT( sizeof( crfXMap::Map::value_type ) == 64 );

		return make_shared< RifFactoryImpl<crfXMap> >( config );
	}
	else if( config.rif_type == "RotScore128" )
	{
		typedef ::scheme::objective::storage::RotamerScore<> crfRotScore;
		typedef ::scheme::objective::storage::RotamerScores< 60, crfRotScore > crfXMapValue;
		BOOST_STATIC_ASSERT( sizeof( crfXMapValue ) == 120 );
		typedef ::scheme::objective::hash::XformMap<
				EigenXform,
				crfXMapValue,
				::scheme::objective::hash::XformHash_bt24_BCC6
			> crfXMap;
		BOOST_STATIC_ASSERT( sizeof( crfXMap::Map::value_type ) == 128 );

		return make_shared< RifFactoryImpl<crfXMap> >( config );
	}
	else if( config.rif_type == "RotScoreSat" )
	{
		typedef ::scheme::objective::storage::RotamerScoreSat<uint16_t, 9, -4> crfRotScore;
		typedef ::scheme::objective::storage::RotamerScores< 14, crfRotScore > crfXMapValue;
		BOOST_STATIC_ASSERT( sizeof( crfXMapValue ) == 56 );
		typedef ::scheme::objective::hash::XformMap<
				EigenXform,
				crfXMapValue,
				::scheme::objective::hash::XformHash_bt24_BCC6
			> crfXMap;
		BOOST_STATIC_ASSERT( sizeof( crfXMap::Map::value_type ) == 64 );

		return make_shared< RifFactoryImpl<crfXMap> >( config );
	}
	else if( config.rif_type == "RotScoreSat_1x16" )
	{
		using SatDatum = ::scheme::objective::storage::SatisfactionDatum<uint16_t>;
		typedef ::scheme::objective::storage::RotamerScoreSat<
					uint16_t, 9, -13, SatDatum, 1> crfRotScore;
		typedef ::scheme::objective::storage::RotamerScores< 14, crfRotScore > crfXMapValue;
		BOOST_STATIC_ASSERT( sizeof( crfXMapValue ) == 56 );
		typedef ::scheme::objective::hash::XformMap<
				EigenXform,
				crfXMapValue,
				::scheme::objective::hash::XformHash_bt24_BCC6
			> crfXMap;
		BOOST_STATIC_ASSERT( sizeof( crfXMap::Map::value_type ) == 64 );

		return make_shared< RifFactoryImpl<crfXMap> >( config );

	} else
	{
		utility_exit_with_message( "create_rif_factory_inner: unknown rif type "+config.rif_type );
	}
}







}}
