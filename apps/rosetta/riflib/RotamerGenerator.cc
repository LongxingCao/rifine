// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols

#include <riflib/RotamerGenerator.hh>
#include <riflib/util.hh>

#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/orbitals/OrbitalType.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/lkball/LK_BallInfo.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>
#include <core/pack/task/TaskFactory.hh>
#include <utility/graph/Graph.hh>
#include <core/pack/packer_neighbors.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/option.hh>

#include <scheme/actor/BackboneActor.hh>

#include <ObjexxFCL/format.hh>
#include <utility/io/ozstream.hh>

namespace devel {
namespace scheme {

using std::endl;
using std::cout;

struct RichardsonRotData {
	// ptp180   11      62   180    65  -175   14 17 10 13
	std::string name;
	float count, chi1, chi2, chi3, chi4, sd1, sd2, sd3, sd4;
	int nchi;
	typedef float F;
	typedef std::string S;
	RichardsonRotData( S n, F c, F c1, F c2, F c3, F c4, F s1, F s2, F s3, F s4 ) : name(n), count(c), chi1(c1), chi2(c2), chi3(c3), chi4(c4), sd1(s1), sd2(s2), sd3(s3), sd4(s4) { nchi=4; }
	RichardsonRotData( S n, F c, F c1, F c2, F c3, F s1, F s2, F s3 ) : name(n), count(c), chi1(c1), chi2(c2), chi3(c3), sd1(s1), sd2(s2), sd3(s3) { nchi=3; }
	RichardsonRotData( S n, F c, F c1, F c2, F s1, F s2 ) : name(n), count(c), chi1(c1), chi2(c2), sd1(s1), sd2(s2) { nchi=2; }
	RichardsonRotData( S n, F c, F c1, F s1 ) : name(n), count(c), chi1(c1), sd1(s1) { nchi=1; }

	std::vector<float>
	get_chi() const {
		std::vector<float> chi;
		if( nchi > 0 ) chi.push_back(chi1);
		if( nchi > 1 ) chi.push_back(chi2);
		if( nchi > 2 ) chi.push_back(chi3);
		if( nchi > 3 ) chi.push_back(chi4);
		return chi;
	}

	void
	get_ex_chis(
		bool ex1,
		bool ex2,
		bool ex3,
		bool ex4,
		std::vector<std::vector<float> > & exchis,
		bool include_primary = true
	) const {
		int nex1=1, nex2=1, nex3=1, nex4=1;
		if( nchi > 0 && ex1 ) nex1 = 3;
		if( nchi > 1 && ex2 ) nex2 = 3;
		if( nchi > 2 && ex3 ) nex3 = 3;
		if( nchi > 3 && ex4 ) nex4 = 3;

		for( int i1 = 0; i1 < nex1; ++i1 ){ float exchi1 = chi1; if( i1==1 ) exchi1 += sd1; else if( i1==2 ) exchi1 -= sd1;
		for( int i2 = 0; i2 < nex2; ++i2 ){ float exchi2 = chi2; if( i2==1 ) exchi2 += sd2;	else if( i2==2 ) exchi2 -= sd2;
		for( int i3 = 0; i3 < nex3; ++i3 ){	float exchi3 = chi3; if( i3==1 ) exchi3 += sd3;	else if( i3==2 ) exchi3 -= sd3;
		for( int i4 = 0; i4 < nex4; ++i4 ){	float exchi4 = chi4; if( i4==1 ) exchi4 += sd4; else if( i4==2 ) exchi4 -= sd4;
			if( include_primary || i1 || i2 || i3 || i4 ){
				std::vector<float> chi;
				if( nchi > 0 ) chi.push_back(exchi1);
				if( nchi > 1 ) chi.push_back(exchi2);
				if( nchi > 2 ) chi.push_back(exchi3);
				if( nchi > 3 ) chi.push_back(exchi4);
				exchis.push_back( chi );
			} else {
				// std::cout << "exclude primary " << i1 << " " << i2 << " " << i3 << " " << i4 << std::endl;
			}
		}}}}

	}

};


	void
	RosettaRotamerGenerator::get_atoms(
		std::string resname,
		std::vector<float> const & chi,
		std::vector<SchemeAtom> & atoms,
		std::vector<std::pair<SchemeAtom,SchemeAtom> > & hbonders,
		int & nheavyatoms,
		std::vector<HBondRay> & donors,
		std::vector<HBondRay> & acceptors
	){

		std::vector<int> atypemap = get_rif_atype_map();

		core::chemical::ResidueTypeSetCAP rts = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");
		core::chemical::ResidueType const & rtype = rts.lock()->name_map( resname );
		core::conformation::ResidueOP resop = core::conformation::ResidueFactory::create_residue( rtype );
		core::pose::Pose pose;
		pose.append_residue_by_jump(*resop,1);
		runtime_assert( chi.size() == resop->nchi() );
		for(int i = 0; i < chi.size(); ++i){
			pose.set_chi( i+1, 1, chi[i] );
		}

		Eigen::Vector3f N ( resop->xyz("N" ).x(), resop->xyz("N" ).y(), resop->xyz("N" ).z() );
		Eigen::Vector3f CA( resop->xyz("CA").x(), resop->xyz("CA").y(), resop->xyz("CA").z() );
		Eigen::Vector3f C ( resop->xyz("C" ).x(), resop->xyz("C" ).y(), resop->xyz("C" ).z() );
		// typedef Eigen::Transform<float,3,Eigen::AffineCompact> EigenXform;
		::scheme::actor::BackboneActor<EigenXform> bbactor( N, CA , C );

		nheavyatoms = resop->nheavyatoms()-1; // remove O
		for(int ia = 1; ia <= resop->natoms(); ++ia){
			if( resop->atom_name(ia)==" O  " ) continue; // remove bb O, shouldn't ever have OXT
			if( resop->atom_name(ia)==" H  " ) continue; // remove bb O, shouldn't ever have OXT
			SchemeAtom::Position p;
			p[0] = pose.residue(1).xyz(ia)[0];
			p[1] = pose.residue(1).xyz(ia)[1];
			p[2] = pose.residue(1).xyz(ia)[2];
			p = bbactor.position().inverse() * p;
				// Position const &    p,
				// int                 type     = 0,
				// std::string const & atomname = AtomData::default_atomname(),
				// std::string const & resname  = AtomData::default_resname(),
				// char                chain    = AtomData::default_chain(),
				// int                 resnum   = AtomData::default_resnum(),
				// int                 atomnum  = AtomData::default_atomnum(),
				// std::string const & elem     = AtomData::default_elem(),
				// bool                ishet    = AtomData::default_ishet(),
				// float               occ      = AtomData::default_occ(),
				// float               bfac     = AtomData::default_bfac()
			int iatype = atypemap[resop->atom_type_index(ia)];
			if( ia <= nheavyatoms+1 && iatype == 12345 ){
				--nheavyatoms; //  THIS IS A BIG HACK!!!!
				continue; // only pro NV??
			}

			SchemeAtom a(
				p,
				iatype,
											resop->atom_name(ia),
											resop->name3(),
											' ',
											1,
											ia,
											"",
											false,
											1.0,
											0.0
			);
			atoms.push_back(a);
		}

		HBRayOpts hbopt;
		hbopt.withbb = false;
		get_donor_rays   ( pose, 1, hbopt, donors );
		get_acceptor_rays( pose, 1, hbopt, acceptors );
		// std::cout << bbactor.position().rotation() << std::endl;
		// std::cout << bbactor.position().translation().transpose() << std::endl;
		for( int i = 0; i < donors.size(); ++i ){
			donors[i].horb_cen  = bbactor.position().inverse()            * donors[i].horb_cen;
			donors[i].direction = bbactor.position().inverse().rotation() * donors[i].direction;
		}
		for( int i = 0; i < acceptors.size(); ++i ){
			acceptors[i].horb_cen  = bbactor.position().inverse()            * acceptors[i].horb_cen;
			acceptors[i].direction = bbactor.position().inverse().rotation() * acceptors[i].direction;
		}


		// get acceptors
		int count = 0;
		for( auto i : rtype.accpt_pos_sc() ){
			for( auto j : rtype.bonded_orbitals( i ) ){
				if(
					// rtype.orbital_type(j).orbital_enum() == core::chemical::orbitals::C_pi_sp2    ||
				    // rtype.orbital_type(j).orbital_enum() == core::chemical::orbitals::N_pi_sp2    ||
				    rtype.orbital_type(j).orbital_enum() == core::chemical::orbitals::N_p_sp2     ||
				    // rtype.orbital_type(j).orbital_enum() == core::chemical::orbitals::O_pi_sp2    ||
				    rtype.orbital_type(j).orbital_enum() == core::chemical::orbitals::O_p_sp2     ||
				    // rtype.orbital_type(j).orbital_enum() == core::chemical::orbitals::S_p_sp3     ||
				    // rtype.orbital_type(j).orbital_enum() == core::chemical::orbitals::O_pi_sp2_bb ||
				    // rtype.orbital_type(j).orbital_enum() == core::chemical::orbitals::O_p_sp2_bb  ||
				    rtype.orbital_type(j).orbital_enum() == core::chemical::orbitals::O_p_sp3
				){
					// std::cout << rtype.name3() << " ACPTR " << rtype.atom_name(i) << " " << rtype.orbital_type(j).name() << " " << i << " " << j << std::endl;
					// std::cout << pose.residue(1).xyz(i) << std::endl;
					// std::cout << pose.residue(1).orbital_xyz(j) << std::endl;
					SchemeAtom::Position pa;
					pa[0] = pose.residue(1).xyz(i)[0];
					pa[1] = pose.residue(1).xyz(i)[1];
					pa[2] = pose.residue(1).xyz(i)[2];
					pa = bbactor.position().inverse() * pa;
					SchemeAtom a(
						pa, //pose.residue(1).xyz(i),          // Position const &    p,
						atypemap[resop->atom_type_index(i)],  // int                 type     = 0,
						resop->atom_name(i),             // std::string const & atomname = AtomData::default_atomname(),
						resop->name3(),                  // std::string const & resname  = AtomData::default_resname(),
						'O',                             // char                chain    = AtomData::default_chain(),
						1,                               // int                 resnum   = AtomData::default_resnum(),
						100+2*count,                    // int                 atomnum  = AtomData::default_atomnum(),
						"",                              // std::string const & elem     = AtomData::default_elem(),
						false,                           // bool                ishet    = AtomData::default_ishet(),
						1.0,                             // float               occ      = AtomData::default_occ(),
						0.0                              // float               bfac     = AtomData::default_bfac()
					);
					SchemeAtom::Position po;
					po[0] = pose.residue(1).orbital_xyz(j)[0];
					po[1] = pose.residue(1).orbital_xyz(j)[1];
					po[2] = pose.residue(1).orbital_xyz(j)[2];
					po = bbactor.position().inverse() * po;
					SchemeAtom o(
						po, //pose.residue(1).orbital_xyz(j),  // Position const &    p,
						-1,                              // int                 type     = 0,
						"ORBL",                          // std::string const & atomname = AtomData::default_atomname(),
						resop->name3(),                  // std::string const & resname  = AtomData::default_resname(),
						'O',                             // char                chain    = AtomData::default_chain(),
						1,                               // int                 resnum   = AtomData::default_resnum(),
						100+2*count+1,                  // int                 atomnum  = AtomData::default_atomnum(),
						"",                              // std::string const & elem     = AtomData::default_elem(),
						false,                           // bool                ishet    = AtomData::default_ishet(),
						1.0,                             // float               occ      = AtomData::default_occ(),
						0.0                              // float               bfac     = AtomData::default_bfac()
					);
					hbonders.push_back( std::make_pair( a, o ) );

					runtime_assert( o.type() == -1 );

					++count;
				}
			}
		}

		// get donors
		for( auto i : rtype.Hpos_polar_sc() ){
			core::Size ib = rtype.atom_base(i);
			// std::cout << rtype.name3() << " DONOR " << rtype.atom_name(i) << " " << rtype.atom_name(ib) << std::endl;
			SchemeAtom::Position pb;
			pb[0] = pose.residue(1).xyz(ib)[0];
			pb[1] = pose.residue(1).xyz(ib)[1];
			pb[2] = pose.residue(1).xyz(ib)[2];
			pb = bbactor.position().inverse() * pb;
			SchemeAtom b(
				pb, //pose.residue(1).xyz(ib),         // Position const &    p,
				atypemap[resop->atom_type_index(ib)], // int                 type     = 0,
				resop->atom_name(ib),            // std::string const & atomname = AtomData::default_atomname(),
				resop->name3(),                  // std::string const & resname  = AtomData::default_resname(),
				'H',                             // char                chain    = AtomData::default_chain(),
				1,                               // int                 resnum   = AtomData::default_resnum(),
				100+2*count,                     // int                 atomnum  = AtomData::default_atomnum(),
				"",                              // std::string const & elem     = AtomData::default_elem(),
				false,                           // bool                ishet    = AtomData::default_ishet(),
				1.0,                             // float               occ      = AtomData::default_occ(),
				0.0                              // float               bfac     = AtomData::default_bfac()
			);
			SchemeAtom::Position ph;
			ph[0] = pose.residue(1).xyz(i)[0];
			ph[1] = pose.residue(1).xyz(i)[1];
			ph[2] = pose.residue(1).xyz(i)[2];
			ph = bbactor.position().inverse() * ph;
			SchemeAtom h(
				ph, //pose.residue(1).xyz(i),  // Position const &    p,
				12345,  // int                 type     = 0,
				resop->atom_name(i),             // std::string const & atomname = AtomData::default_atomname(),
				resop->name3(),                  // std::string const & resname  = AtomData::default_resname(),
				'H',                             // char                chain    = AtomData::default_chain(),
				1,                               // int                 resnum   = AtomData::default_resnum(),
				100+2*count+1,                   // int                 atomnum  = AtomData::default_atomnum(),
				"",                              // std::string const & elem     = AtomData::default_elem(),
				false,                           // bool                ishet    = AtomData::default_ishet(),
				1.0,                             // float               occ      = AtomData::default_occ(),
				0.0                              // float               bfac     = AtomData::default_bfac()
			);
			hbonders.push_back( std::make_pair( b, h ) );

			++count;

		}

	}

core::pack::rotamer_set::RotamerSetOP
get_rosetta_rot_set(
	core::pose::Pose & pose,
	core::Size ir
 ){
	core::scoring::ScoreFunction dummy_sfxn;
	dummy_sfxn.score( pose );
	core::pack::task::PackerTaskOP dummy_task = core::pack::task::TaskFactory::create_packer_task( pose );
	dummy_task->initialize_from_command_line();
	dummy_task->initialize_extra_rotamer_flags_from_command_line();
	dummy_task->nonconst_residue_task( ir ).and_extrachi_cutoff(0);
	dummy_task->nonconst_residue_task( ir ).restrict_to_repacking();
	dummy_task->nonconst_residue_task( ir ).or_include_current( false ); //need to do this because the residue was built from internal coords and is probably crumpled up
	dummy_task->nonconst_residue_task( ir ).or_fix_his_tautomer( true ); //since we only want rotamers for the specified restype
	utility::graph::GraphOP dummy_png = core::pack::create_packer_graph( pose, dummy_sfxn, dummy_task );
	core::pack::rotamer_set::RotamerSetFactory rsf;
	core::pack::rotamer_set::RotamerSetOP rotset( rsf.create_rotamer_set( pose.residue( ir ) ) );
	rotset->set_resid( ir );
	rotset->build_rotamers( pose, dummy_sfxn, *dummy_task, dummy_png );
	return rotset;
 }




void
get_richardson_rot_data(
	std::map<std::string,std::vector<RichardsonRotData> > & rrdata,
	bool extra_chi2 = true
){
	bool const ER = extra_chi2;
	// if( extra_chi2 ){
	// 	std::cout << "get_richardson_rot_data: USING EXTRA ROTS" << std::endl;
	// 	std::cout << "extra rotamers aren't well tested yet, don't use" << std::endl;
	// 	std::exit(-1);
	// }

	rrdata.insert( std::make_pair( "ARG", std::vector<RichardsonRotData>() ) );

	rrdata["ARG"].push_back( RichardsonRotData( "ptp85b ",   3,      62,   180,    65,    85,   14, 17, 10, 13 ) );
	rrdata["ARG"].push_back( RichardsonRotData( "ptp180 ",  11,      62,   180,    65,  -175,   14, 17, 10, 13 ) );

	rrdata["ARG"].push_back( RichardsonRotData( "ptt85  ",  16,      62,   180,   180,    85,   13, 14, 13, 17 ) );
	rrdata["ARG"].push_back( RichardsonRotData( "ptt180 ",  16,      62,   180,   180,   180,   15, 13, 15, 19 ) );
	rrdata["ARG"].push_back( RichardsonRotData( "ptt-85 ",  15,      62,   180,   180,   -85,   15, 14, 12, 14 ) );

	rrdata["ARG"].push_back( RichardsonRotData( "ptm180 ",   6,      62,   180,   -65,   175,   13, 13, 12, 15 ) );
	rrdata["ARG"].push_back( RichardsonRotData( "ptm-85 ",   5,      62,   180,   -65,   -85,   13, 13, 12, 15 ) );

	rrdata["ARG"].push_back( RichardsonRotData( "tpp85  ",  11,    -177,    65,    65,    85,   13, 13, 12, 15 ) );
	rrdata["ARG"].push_back( RichardsonRotData( "tpp180 ",   8,    -177,    65,    65,  -175,   13, 13, 12, 15 ) );

	rrdata["ARG"].push_back( RichardsonRotData( "tpt85  ",  20,    -177,    65,   180,    85,   14, 13, 15, 14 ) );
	rrdata["ARG"].push_back( RichardsonRotData( "tpt180 ",  15,    -177,    65,   180,   180,   13, 17, 14, 17 ) );

	rrdata["ARG"].push_back( RichardsonRotData( "ttp85  ",  33,    -177,   180,    65,    85,   14, 17, 13, 15 ) );
	rrdata["ARG"].push_back( RichardsonRotData( "ttp180 ",  25,    -177,   180,    65,  -175,   14, 16, 14, 16 ) );
	rrdata["ARG"].push_back( RichardsonRotData( "ttp105 ",   9,    -177,   180,    65,  -105,   14, 16, 14, 16 ) );

	rrdata["ARG"].push_back( RichardsonRotData( "ttt85  ",  19,    -177,   180,   180,    85,   14, 14, 13, 14 ) );
	rrdata["ARG"].push_back( RichardsonRotData( "ttt180 ",  33,    -177,   180,   180,   180,   15, 13, 12, 27 ) );
	rrdata["ARG"].push_back( RichardsonRotData( "ttt-85 ",  26,    -177,   180,   180,   -85,   15, 14, 14, 15 ) );
	rrdata["ARG"].push_back( RichardsonRotData( "ttm105 ",  10,    -177,   180,   -65,   105,   15, 16, 15, 15 ) );
	rrdata["ARG"].push_back( RichardsonRotData( "ttm180 ",  13,    -177,   180,   -65,   175,   15, 12, 11, 15 ) );
	rrdata["ARG"].push_back( RichardsonRotData( "ttm-85 ",  28,    -177,   180,   -65,   -85,   14, 16, 15, 14 ) );
	rrdata["ARG"].push_back( RichardsonRotData( "mtp85  ",  22,     -67,   180,    65,    85,   13, 17, 13, 13 ) );

	if(ER) rrdata["ARG"].push_back( RichardsonRotData( "mtp180 ", -45,     -67,   161,    65,  -175,   12, 19, 13, 19 ) );
	       rrdata["ARG"].push_back( RichardsonRotData( "mtp180 ", -45,     -67,   180,    65,  -175,   12, 19, 13, 19 ) );
	if(ER) rrdata["ARG"].push_back( RichardsonRotData( "mtp180 ", -45,     -67,  -161,    65,  -175,   12, 19, 13, 19 ) );

	rrdata["ARG"].push_back( RichardsonRotData( "mtp-105",   7,     -67,   180,    65,  -105,   11, 15, 13, 15 ) );
	rrdata["ARG"].push_back( RichardsonRotData( "mtt85  ",  34,     -67,   180,   180,    85,   12, 19, 13, 19 ) );

	if(ER) rrdata["ARG"].push_back( RichardsonRotData( "mtt180 ", -89,     -67,   167,   180,   180,   14, 13, 13, 21 ) );
	       rrdata["ARG"].push_back( RichardsonRotData( "mtt180 ", -89,     -67,   180,   180,   180,   14, 13, 13, 21 ) );
	if(ER) rrdata["ARG"].push_back( RichardsonRotData( "mtt180 ", -89,     -67,  -167,   180,   180,   14, 13, 13, 21 ) );

	if(ER) rrdata["ARG"].push_back( RichardsonRotData( "mtt-85 ", -53,     -67,   167,   180,   -85,   13, 13, 13, 13 ) );
	       rrdata["ARG"].push_back( RichardsonRotData( "mtt-85 ", -53,     -67,   180,   180,   -85,   13, 13, 13, 13 ) );
	if(ER) rrdata["ARG"].push_back( RichardsonRotData( "mtt-85 ", -53,     -67,  -167,   180,   -85,   13, 13, 13, 13 ) );

	rrdata["ARG"].push_back( RichardsonRotData( "mtm105 ",  15,     -67,   180,   -65,   105,   12, 13, 13, 15 ) );
	rrdata["ARG"].push_back( RichardsonRotData( "mtm180 ",  48,     -67,   180,   -65,   175,   14, 17, 13, 29 ) );

	if(ER) rrdata["ARG"].push_back( RichardsonRotData( "mtm-85 ", -54,     -67,  -154,   -65,   -85,   14, 13, 13, 13 ) );
	       rrdata["ARG"].push_back( RichardsonRotData( "mtm-85 ", -54,     -67,  -167,   -65,   -85,   14, 13, 13, 13 ) );
	if(ER) rrdata["ARG"].push_back( RichardsonRotData( "mtm-85 ", -54,     -67,   180,   -65,   -85,   14, 13, 13, 13 ) );

	rrdata["ARG"].push_back( RichardsonRotData( "mmt85  ",   7,     -62,   -68,   180,    85,   13, 13, 10, 16 ) );
	rrdata["ARG"].push_back( RichardsonRotData( "mmt180 ",  18,     -62,   -68,   180,   180,   13, 13, 10, 16 ) );
	rrdata["ARG"].push_back( RichardsonRotData( "mmt85  ",  22,     -62,   -68,   180,   -85,   14, 13, 15, 13 ) );

	rrdata["ARG"].push_back( RichardsonRotData( "mmm180 ",  11,     -62,   -68,   -65,   175,   14, 15, 10, 13 ) );
	rrdata["ARG"].push_back( RichardsonRotData( "mmm-85 ",  22,     -62,   -68,   -65,   -85,   14, 13, 15, 13 ) );

			// Lysine
	rrdata["LYS"].push_back( RichardsonRotData( "ptpt",      7,      62,   180,    68,   180,   13, 14, 14, 11 ) );
	rrdata["LYS"].push_back( RichardsonRotData( "pttp",     13,      62,   180,   180,    65,   13, 14, 14, 11 ) );
	rrdata["LYS"].push_back( RichardsonRotData( "pttt",     29,      62,   180,   180,   180,   13, 13, 13, 10 ) );
	rrdata["LYS"].push_back( RichardsonRotData( "pttm",      8,      62,   180,   180,   -65,   13, 13, 13, 10 ) );
	rrdata["LYS"].push_back( RichardsonRotData( "ptmt",      5,      62,   180,   -68,   180,   13, 13, 13, 10 ) );
	rrdata["LYS"].push_back( RichardsonRotData( "tptp",     11,    -177,    68,   180,    65,   13, 12, 10, 11 ) );
	rrdata["LYS"].push_back( RichardsonRotData( "tptt",     32,    -177,    68,   180,   180,   10, 10, 13, 14 ) );
	rrdata["LYS"].push_back( RichardsonRotData( "tptm",      7,    -177,    68,   180,   -65,   14,  9, 12, 10 ) );
	rrdata["LYS"].push_back( RichardsonRotData( "ttpp",     12,    -177,   180,    68,    65,   14, 12, 14, 14 ) );
	rrdata["LYS"].push_back( RichardsonRotData( "ttpt",     25,    -177,   180,    68,   180,   14, 12, 14, 14 ) );
	rrdata["LYS"].push_back( RichardsonRotData( "tttp",     49,    -177,   180,   180,    65,   14, 13, 12, 12 ) );

	if(ER) rrdata["LYS"].push_back( RichardsonRotData( "tttt",   -162,    -177,   167,   180,   180,   13, 13, 15, 13 ) );
	       rrdata["LYS"].push_back( RichardsonRotData( "tttt",   -162,    -177,   180,   180,   180,   13, 13, 15, 13 ) );
	if(ER) rrdata["LYS"].push_back( RichardsonRotData( "tttt",   -162,    -177,  -167,   180,   180,   13, 13, 15, 13 ) );

	rrdata["LYS"].push_back( RichardsonRotData( "tttm",     37,    -177,   180,   180,   -65,   12, 13, 15, 13 ) );
	rrdata["LYS"].push_back( RichardsonRotData( "ttmt",     20,    -177,   180,   -68,   180,   14, 14, 10, 15 ) );
	rrdata["LYS"].push_back( RichardsonRotData( "ttmm",      5,    -177,   180,   -68,   -65,   14, 14, 10, 15 ) );
	rrdata["LYS"].push_back( RichardsonRotData( "mptt",      4,     -90,    68,   180,   180,   10,  9, 10, 13 ) );
	rrdata["LYS"].push_back( RichardsonRotData( "mtpp",     12,     -67,   180,    68,    65,   12, 13, 11,  9 ) );
	rrdata["LYS"].push_back( RichardsonRotData( "mtpt",     38,     -67,   180,    68,   180,   13, 13, 14, 14 ) );
	rrdata["LYS"].push_back( RichardsonRotData( "mttp",     42,     -67,   180,   180,    65,   13, 13, 14, 14 ) );

	if(ER) rrdata["LYS"].push_back( RichardsonRotData( "mttt",   -244,     -67,   167,   180,   180,   14, 13, 12, 14 ) );
	       rrdata["LYS"].push_back( RichardsonRotData( "mttt",   -244,     -67,   180,   180,   180,   14, 13, 12, 14 ) );
	if(ER) rrdata["LYS"].push_back( RichardsonRotData( "mttt",   -244,     -67,  -167,   180,   180,   14, 13, 12, 14 ) );

	rrdata["LYS"].push_back( RichardsonRotData( "mttm",     56,     -67,   180,   180,   -65,   13, 12, 13, 14 ) );
	rrdata["LYS"].push_back( RichardsonRotData( "mtmt",     40,     -67,   180,   -68,   180,   12, 13, 14, 13 ) );
	rrdata["LYS"].push_back( RichardsonRotData( "mtmm",     12,     -67,   180,   -68,   -65,   12, 12, 12, 11 ) );
	rrdata["LYS"].push_back( RichardsonRotData( "mmtp",      9,     -62,   -68,   180,    65,   12, 13, 13, 13 ) );
	rrdata["LYS"].push_back( RichardsonRotData( "mmtt",     77,     -62,   -68,   180,   180,   12, 13, 13, 13 ) );
	rrdata["LYS"].push_back( RichardsonRotData( "mmtm",     18,     -62,   -68,   180,   -65,   14, 12, 10, 15 ) );
	rrdata["LYS"].push_back( RichardsonRotData( "mmmt",     10,     -62,   -68,   -68,   180,   12, 13, 10, 15 ) );

			// Methionine
	rrdata["MET"].push_back( RichardsonRotData( "ptp",      12,      62,   180,    75,         11, 17, 12 ) );
	rrdata["MET"].push_back( RichardsonRotData( "ptm",      17,      62,   180,   -75,          9, 10,  9 ) );
	rrdata["MET"].push_back( RichardsonRotData( "tpp",      30,    -177,    65,    75,         10, 15, 15 ) );
	rrdata["MET"].push_back( RichardsonRotData( "tpt",       9,    -177,    65,   180,          9,  8,  9 ) );
	rrdata["MET"].push_back( RichardsonRotData( "ttp",      28,    -177,   180,    75,         10, 11, 11 ) );
	rrdata["MET"].push_back( RichardsonRotData( "ttt",      17,    -177,   180,   180,          9,  9, 19 ) );
	rrdata["MET"].push_back( RichardsonRotData( "ttm",      36,    -177,   180,   -75,         10, 10, 13 ) );
	rrdata["MET"].push_back( RichardsonRotData( "mtt",      43,     -67,   180,   180,         10, 13, 15 ) );
	rrdata["MET"].push_back( RichardsonRotData( "mmp",      15,     -65,   -65,   103,          9, 10, 10 ) );
	rrdata["MET"].push_back( RichardsonRotData( "mmt",      10,     -65,   -65,   180,         12, 14, 19 ) );

	if(ER) rrdata["MET"].push_back( RichardsonRotData( "mtm",     -58,     -67,   165,   -75,         12, 11, 16 ) );
	       rrdata["MET"].push_back( RichardsonRotData( "mtm",     -58,     -67,   180,   -75,         12, 11, 16 ) );
	if(ER) rrdata["MET"].push_back( RichardsonRotData( "mtm",     -58,     -67,  -165,   -75,         12, 11, 16 ) );
	if(ER) rrdata["MET"].push_back( RichardsonRotData( "mtp",     -92,     -67,   165,    75,         10, 12, 14 ) );
	       rrdata["MET"].push_back( RichardsonRotData( "mtp",     -92,     -67,   180,    75,         10, 12, 14 ) );
	if(ER) rrdata["MET"].push_back( RichardsonRotData( "mtp",     -92,     -67,  -165,    75,         10, 12, 14 ) );
	if(ER) rrdata["MET"].push_back( RichardsonRotData( "mmm",    -105,     -65,   -80,   -70,         11, 13, 16 ) );
	       rrdata["MET"].push_back( RichardsonRotData( "mmm",    -105,     -65,   -65,   -70,         11, 13, 16 ) );
	if(ER) rrdata["MET"].push_back( RichardsonRotData( "mmm",    -105,     -65,   -50,   -70,         11, 13, 16 ) );


			// Glutamate
	rrdata["GLU"].push_back( RichardsonRotData( "Spt-60",    27,     62,   180,   -60,         14, 13, 20 ) );//
	rrdata["GLU"].push_back( RichardsonRotData( "pt-20 ",    27,     62,   180,     0,         14, 13, 20 ) );//   -90 to 90
	rrdata["GLU"].push_back( RichardsonRotData( "Spt60 ",    27,     62,   180,    60,         14, 13, 20 ) );//
	rrdata["GLU"].push_back( RichardsonRotData( "pm0   ",    32,     70,   -80,     0,         14, 13, 17 ) );//   -50 to 50
	rrdata["GLU"].push_back( RichardsonRotData( "tp10  ",    91,   -177,    65,    10,         14, 13, 17 ) );//   -10 to 90
	rrdata["GLU"].push_back( RichardsonRotData( "tm20  ",    17,   -177,   -80,   -25,         13, 13, 15 ) );//   -50 to 10
	rrdata["GLU"].push_back( RichardsonRotData( "mp0   ",    88,    -65,    85,     0,         14, 13, 25 ) );//   -60 to 60

	if(ER) rrdata["GLU"].push_back( RichardsonRotData( "Stt-60",  -117,   -177,   166,   -60,         14, 14, 20 ) );//
	       rrdata["GLU"].push_back( RichardsonRotData( "Stt-60",  -117,   -177,   180,   -60,         14, 14, 20 ) );//
	if(ER) rrdata["GLU"].push_back( RichardsonRotData( "Stt-60",  -117,   -177,  -166,   -60,         14, 14, 20 ) );//

	if(ER) rrdata["GLU"].push_back( RichardsonRotData( "tt0   ",  -117,   -177,   166,     0,         14, 14, 20 ) );//   -90 to 90
	       rrdata["GLU"].push_back( RichardsonRotData( "tt0   ",  -117,   -177,   180,     0,         14, 14, 20 ) );//   -90 to 90
	if(ER) rrdata["GLU"].push_back( RichardsonRotData( "tt0   ",  -117,   -177,  -166,     0,         14, 14, 20 ) );//   -90 to 90

	if(ER) rrdata["GLU"].push_back( RichardsonRotData( "Stt60 ",  -117,   -177,   166,    60,         14, 14, 20 ) );//
	       rrdata["GLU"].push_back( RichardsonRotData( "Stt60 ",  -117,   -177,   180,    60,         14, 14, 20 ) );//
	if(ER) rrdata["GLU"].push_back( RichardsonRotData( "Stt60 ",  -117,   -177,  -166,    60,         14, 14, 20 ) );//

	if(ER) rrdata["GLU"].push_back( RichardsonRotData( "Smt60 ",  -161,    -67,   166,   -60,         13, 16, 20 ) );//
	       rrdata["GLU"].push_back( RichardsonRotData( "Smt60 ",  -161,    -67,   180,   -60,         13, 16, 20 ) );//
	if(ER) rrdata["GLU"].push_back( RichardsonRotData( "Smt60 ",  -161,    -67,  -166,   -60,         13, 16, 20 ) );//

	if(ER) rrdata["GLU"].push_back( RichardsonRotData( "mt10  ",  -161,    -67,   166,     0,         13, 16, 20 ) );//   -90 to 90
	       rrdata["GLU"].push_back( RichardsonRotData( "mt10  ",  -161,    -67,   180,     0,         13, 16, 20 ) );//   -90 to 90
	if(ER) rrdata["GLU"].push_back( RichardsonRotData( "mt10  ",  -161,    -67,  -166,     0,         13, 16, 20 ) );//   -90 to 90

	if(ER) rrdata["GLU"].push_back( RichardsonRotData( "Smt60 ",  -161,    -67,   166,    60,         13, 16, 20 ) );//
	       rrdata["GLU"].push_back( RichardsonRotData( "Smt60 ",  -161,    -67,   180,    60,         13, 16, 20 ) );//
	if(ER) rrdata["GLU"].push_back( RichardsonRotData( "Smt60 ",  -161,    -67,  -166,    60,         13, 16, 20 ) );//

	rrdata["GLU"].push_back( RichardsonRotData( "mm40  ",    98,    -65,   -65,   -60,         14, 14, 20 ) );//   -90 to 30
	rrdata["GLU"].push_back( RichardsonRotData( "Smm0  ",    98,    -65,   -75,     0,         14, 14, 20 ) );//


			// Glutamine
	rrdata["GLN"].push_back( RichardsonRotData( "Spt60",     12,      62,  180,   -60,         13, 14, 20 ) );//
	rrdata["GLN"].push_back( RichardsonRotData( "pt20 ",     12,      62,  180,     0,         13, 14, 20 ) );//   -90 to 90
	rrdata["GLN"].push_back( RichardsonRotData( "Spt60",     12,      62,  180,    60,         13, 14, 20 ) );//
	rrdata["GLN"].push_back( RichardsonRotData( "pm0  ",     15,      70,  -75,     0,         13, 14, 16 ) );//   -60 to 60
	rrdata["GLN"].push_back( RichardsonRotData( "tp100",     14,    -177,   65,  -100,         13, 14, 16 ) );//   -150 to 0
	rrdata["GLN"].push_back( RichardsonRotData( "Stt60",     47,    -177,  180,   -60,         14, 13, 20 ) );//
	rrdata["GLN"].push_back( RichardsonRotData( "tt0  ",     47,    -177,  180,     0,         14, 13, 20 ) );//   -90 to 90
	rrdata["GLN"].push_back( RichardsonRotData( "Stt60",     47,    -177,  180,    60,         14, 13, 20 ) );//
	rrdata["GLN"].push_back( RichardsonRotData( "mp0  ",     24,     -65,   85,     0,         14, 13, 29 ) );//   -60 to 60
	rrdata["GLN"].push_back( RichardsonRotData( "mm100",     22,     -65,  -65,   100,         16, 18, 26 ) );//      0 to 150

	if(ER) rrdata["GLN"].push_back( RichardsonRotData( "tp60 ",    -78,    -177,   50,    60,         14, 15, 24 ) );//   0 to 90
	       rrdata["GLN"].push_back( RichardsonRotData( "tp60 ",    -78,    -177,   65,    60,         14, 15, 24 ) );//   0 to 90
	if(ER) rrdata["GLN"].push_back( RichardsonRotData( "tp60 ",    -78,    -177,   80,    60,         14, 15, 24 ) );//   0 to 90

	if(ER) rrdata["GLN"].push_back( RichardsonRotData( "Smt60",   -101,     -67,  165,   -60,         16, 15, 20 ) );//
	       rrdata["GLN"].push_back( RichardsonRotData( "Smt60",   -101,     -67,  180,   -60,         16, 15, 20 ) );//
	if(ER) rrdata["GLN"].push_back( RichardsonRotData( "Smt60",   -101,     -67, -165,   -60,         16, 15, 20 ) );//

	if(ER) rrdata["GLN"].push_back( RichardsonRotData( "mt30 ",   -101,     -67,  165,     0,         16, 15, 20 ) );//    -90 to 90
	       rrdata["GLN"].push_back( RichardsonRotData( "mt30 ",   -101,     -67,  180,     0,         16, 15, 20 ) );//    -90 to 90
	if(ER) rrdata["GLN"].push_back( RichardsonRotData( "mt30 ",   -101,     -67, -165,     0,         16, 15, 20 ) );//    -90 to 90

	if(ER) rrdata["GLN"].push_back( RichardsonRotData( "Smt60",   -101,     -67,  165,    60,         16, 15, 20 ) );//
	       rrdata["GLN"].push_back( RichardsonRotData( "Smt60",   -101,     -67,  180,    60,         16, 15, 20 ) );//
	if(ER) rrdata["GLN"].push_back( RichardsonRotData( "Smt60",   -101,     -67, -165,    60,         16, 15, 20 ) );//

	if(ER) rrdata["GLN"].push_back( RichardsonRotData( "mm-40",   -127,     -65,  -83,   -40,         16, 18, 26 ) );//    -95 to 0
	       rrdata["GLN"].push_back( RichardsonRotData( "mm-40",   -127,     -65,  -65,   -40,         16, 18, 26 ) );//    -95 to 0
	if(ER) rrdata["GLN"].push_back( RichardsonRotData( "mm-40",   -127,     -65,  -47,   -40,         16, 18, 26 ) );//    -95 to 0


			// Aspartate range
	// rrdata["ASP"].push_back( RichardsonRotData( "Sp50",      60,      62,  -50,                9, 19 ) );//
	rrdata["ASP"].push_back( RichardsonRotData( "p10 ",     143,      62,  -30,          10/*9*/, 19 ) );//     -90 to 0
	rrdata["ASP"].push_back( RichardsonRotData( "p30 ",     194,      62,   30,          10/*8*/, 14 ) );//       0 to 90
	rrdata["ASP"].push_back( RichardsonRotData( "St30",     100,    -170,  -10,               12, 13 ) );//
	rrdata["ASP"].push_back( RichardsonRotData( "t0  ",     338,    -177,   10,               12, 13 ) );//     -50 to 50
	rrdata["ASP"].push_back( RichardsonRotData( "t70 ",     118,    -177,   65,               12, 18 ) );//      50 to 90
	rrdata["ASP"].push_back( RichardsonRotData( "m20 ",     888,     -70,  -15,               10, 16 ) );//     -90 to 20
	rrdata["ASP"].push_back( RichardsonRotData( "Sm60",     200,     -65,   60,               10, 16 ) );//


	// 		// Asparagine
	rrdata["ASN"].push_back( RichardsonRotData( "Sp50",      60,      62,  -60,           10/*8*/, 20 ) );//
	rrdata["ASN"].push_back( RichardsonRotData( "p10 ",      60,      62,    0,           10/*8*/, 20 ) );//     -90 to 0
	rrdata["ASN"].push_back( RichardsonRotData( "p30 ",     100,      62,   60,           10/*6*/, 20 ) );//     0 to 90
	rrdata["ASN"].push_back( RichardsonRotData( "St80",      77,    -174,  -60,           10/*5*/, 20 ) );//
	rrdata["ASN"].push_back( RichardsonRotData( "t20 ",     100,    -174,  -30,           10/*5*/, 20 ) );//    -120 to 0

	if(ER) rrdata["ASN"].push_back( RichardsonRotData( "t30 ",    -228,    -177,   15,              14, 15 ) );//    0 to 80
	       rrdata["ASN"].push_back( RichardsonRotData( "t30 ",    -228,    -177,   30,              14, 15 ) );//    0 to 80
	if(ER) rrdata["ASN"].push_back( RichardsonRotData( "t30 ",    -228,    -177,   45,              14, 15 ) );//    0 to 80

	if(ER) rrdata["ASN"].push_back( RichardsonRotData( "m20 ",    -580,     -65,    0,              10, 20 ) );//    -60 to 10
	       rrdata["ASN"].push_back( RichardsonRotData( "m20 ",    -580,     -65,  -20,              10, 20 ) );//    -60 to 10
	if(ER) rrdata["ASN"].push_back( RichardsonRotData( "m20 ",    -580,     -65,  -40,              10, 20 ) );//    -60 to 10

	rrdata["ASN"].push_back( RichardsonRotData( "m80 ",     118,     -65,  -75,              10, 10 ) );//   -100 to -60
	rrdata["ASN"].push_back( RichardsonRotData( "m120",      58,     -65,  120,              10, 18 ) );//     60 to 160

			// Isoleucine
	rrdata["ILE"].push_back( RichardsonRotData( "pp",        10,      62,  100,              10, 10 ) );
	rrdata["ILE"].push_back( RichardsonRotData( "pt",       216,      62,  170,              10, 10 ) );
	rrdata["ILE"].push_back( RichardsonRotData( "tp",        36,    -177,   66,              13, 11 ) );
	rrdata["ILE"].push_back( RichardsonRotData( "tt",       127,    -177,  165,              13, 11 ) );
	rrdata["ILE"].push_back( RichardsonRotData( "mp",        19,     -65,  100,              10, 10 ) );
	rrdata["ILE"].push_back( RichardsonRotData( "mt",       993,     -65,  170,              10, 10 ) );
	rrdata["ILE"].push_back( RichardsonRotData( "mm",       242,     -57,  -60,              10, 10 ) );


			// Leucine
	       rrdata["LEU"].push_back( RichardsonRotData( "pp",        21,      62,   80,             10, 10 ) );//

	if(ER) rrdata["LEU"].push_back( RichardsonRotData( "tp",      -750,    -177,   50,             10, 10 ) );// NOT FOLLOWING STATS
	       rrdata["LEU"].push_back( RichardsonRotData( "tp",      -750,    -177,   65,             10, 10 ) );//
	if(ER) rrdata["LEU"].push_back( RichardsonRotData( "tp",      -750,    -177,   80,             10, 10 ) );// NOT FOLLOWING STATS

	       rrdata["LEU"].push_back( RichardsonRotData( "tt",        49,    -172,  145,             10, 15 ) );//   120 to 180

	       rrdata["LEU"].push_back( RichardsonRotData( "mp",        63,     -85,   65,             11, 14 ) );//    45 to 105

	if(ER) rrdata["LEU"].push_back( RichardsonRotData( "mt",     -1548,     -65,  160,             11, 11 ) );// NOT FOLLOWING STATS
	       rrdata["LEU"].push_back( RichardsonRotData( "mt",     -1548,     -65,  175,             11, 11 ) );//
	if(ER) rrdata["LEU"].push_back( RichardsonRotData( "mt",     -1548,     -65, -170,             11, 11 ) );// NOT FOLLOWING STATS



			// Histidine
	rrdata["HIS"].push_back( RichardsonRotData( "p80 ",      51,      62,  -75,             10, 12 ) ); //  -120 to -50
	rrdata["HIS"].push_back( RichardsonRotData( "p80 ",      26,      62,   80,             13, 10 ) ); //  50 to 120
	rrdata["HIS"].push_back( RichardsonRotData( "t160",      31,    -177, -165,             12, 20 ) ); //   150 to -120
	rrdata["HIS"].push_back( RichardsonRotData( "t80 ",      64,    -177,  -80,             10, 22 ) ); //   -120 to -50

	if(ER) rrdata["HIS"].push_back( RichardsonRotData( "t60 ",     -94,    -177,   41,             13, 19 ) ); //   50 to 120
	       rrdata["HIS"].push_back( RichardsonRotData( "t60 ",     -94,    -177,   60,             13, 19 ) ); //   50 to 120
	if(ER) rrdata["HIS"].push_back( RichardsonRotData( "t60 ",     -94,    -177,   79,             13, 19 ) ); //   50 to 120

	if(ER) rrdata["HIS"].push_back( RichardsonRotData( "m70 ",    -174,     -65,  -47,             11, 23 ) ); //   -120 to -30
	       rrdata["HIS"].push_back( RichardsonRotData( "m70 ",    -174,     -65,  -70,             11, 23 ) ); //   -120 to -30
	if(ER) rrdata["HIS"].push_back( RichardsonRotData( "m70 ",    -174,     -65,  -93,             11, 23 ) ); //   -120 to -30

	rrdata["HIS"].push_back( RichardsonRotData( "m170",      44,     -65,  165,             10, 16 ) ); //   120 to -160

	if(ER) rrdata["HIS"].push_back( RichardsonRotData( "m80 ",     -78,     -65,   62,             11, 18 ) ); //   50 to 120
	       rrdata["HIS"].push_back( RichardsonRotData( "m80 ",     -78,     -65,   80,             11, 18 ) ); //   50 to 120
	if(ER) rrdata["HIS"].push_back( RichardsonRotData( "m80 ",     -78,     -65,   98,             11, 18 ) ); //   50 to 120

			rrdata["HIS_D"].push_back( RichardsonRotData( "p80 ",      51,      62,  -75,             10, 12 ) ); //  -120 to -50
			rrdata["HIS_D"].push_back( RichardsonRotData( "p80 ",      26,      62,   80,             13, 10 ) ); //  50 to 120
			rrdata["HIS_D"].push_back( RichardsonRotData( "t160",      31,    -177, -165,             12, 20 ) ); //   150 to -120
			rrdata["HIS_D"].push_back( RichardsonRotData( "t80 ",      64,    -177,  -80,             10, 22 ) ); //   -120 to -50

			if(ER) rrdata["HIS_D"].push_back( RichardsonRotData( "t60 ",     -94,    -177,   41,             13, 19 ) ); //   50 to 120
			       rrdata["HIS_D"].push_back( RichardsonRotData( "t60 ",     -94,    -177,   60,             13, 19 ) ); //   50 to 120
			if(ER) rrdata["HIS_D"].push_back( RichardsonRotData( "t60 ",     -94,    -177,   79,             13, 19 ) ); //   50 to 120

			if(ER) rrdata["HIS_D"].push_back( RichardsonRotData( "m70 ",    -174,     -65,  -47,             11, 23 ) ); //   -120 to -30
			       rrdata["HIS_D"].push_back( RichardsonRotData( "m70 ",    -174,     -65,  -70,             11, 23 ) ); //   -120 to -30
			if(ER) rrdata["HIS_D"].push_back( RichardsonRotData( "m70 ",    -174,     -65,  -93,             11, 23 ) ); //   -120 to -30

			rrdata["HIS_D"].push_back( RichardsonRotData( "m170",      44,     -65,  165,             10, 16 ) ); //   120 to -160

			if(ER) rrdata["HIS_D"].push_back( RichardsonRotData( "m80 ",     -78,     -65,   62,             11, 18 ) ); //   50 to 120
			       rrdata["HIS_D"].push_back( RichardsonRotData( "m80 ",     -78,     -65,   80,             11, 18 ) ); //   50 to 120
			if(ER) rrdata["HIS_D"].push_back( RichardsonRotData( "m80 ",     -78,     -65,   98,             11, 18 ) ); //   50 to 120



			// Tryptophan // added ex2 "std dev" ex1 handled as "child" rotamer
	if(ER) rrdata["TRP"].push_back( RichardsonRotData( "p90 ",     -67 ,     62, -100,             12, 10 ) );//   -130 to -60
	       rrdata["TRP"].push_back( RichardsonRotData( "p90 ",     -67 ,     62,  -90,             12, 10 ) );//   -130 to -60
	if(ER) rrdata["TRP"].push_back( RichardsonRotData( "p90 ",     -67 ,     62,  -80,             12, 10 ) );//   -130 to -60

	if(ER) rrdata["TRP"].push_back( RichardsonRotData( "p90 ",     -34 ,     62,   70,             12,  8 ) );//   60 to 130
	       rrdata["TRP"].push_back( RichardsonRotData( "p90 ",     -34 ,     62,   90,             12,  8 ) );//   60 to 130
	if(ER) rrdata["TRP"].push_back( RichardsonRotData( "p90 ",     -34 ,     62,  110,             12,  8 ) );//   60 to 130

	if(ER) rrdata["TRP"].push_back( RichardsonRotData( "t105",    -100 ,   -177, -119,             16, 14 ) );//   -130 to -60
	       rrdata["TRP"].push_back( RichardsonRotData( "t105",    -100 ,   -177, -105,             16, 14 ) );//   -130 to -60
	if(ER) rrdata["TRP"].push_back( RichardsonRotData( "t105",    -100 ,   -177,  -92,             16, 14 ) );//   -130 to -60

	if(ER) rrdata["TRP"].push_back( RichardsonRotData( "t90 ",    -109 ,   -177,   79,             10, 11 ) );//   0 to 100
	       rrdata["TRP"].push_back( RichardsonRotData( "t90 ",    -109 ,   -177,   90,             10, 11 ) );//   0 to 100
	if(ER) rrdata["TRP"].push_back( RichardsonRotData( "t90 ",    -109 ,   -177,  101,             10, 11 ) );//   0 to 100

	if(ER) rrdata["TRP"].push_back( RichardsonRotData( "m0  ",     -48 ,    -65,  -40,              9, 20 ) );//   -40 to 20
	       rrdata["TRP"].push_back( RichardsonRotData( "m0  ",     -48 ,    -65,   -5,              9, 20 ) );//   -40 to 20
	       rrdata["TRP"].push_back( RichardsonRotData( "m0  ",     -48 ,    -56,  -51,              9, 20 ) );//TODO: remove extra rotamer for 1a22 TRP
	
	if(ER) rrdata["TRP"].push_back( RichardsonRotData( "m0  ",     -48 ,    -65,   20,              9, 20 ) );//   -40 to 20

	if(ER) rrdata["TRP"].push_back( RichardsonRotData( "m95 ",    -195 ,    -65,   76,             11, 19 ) );//   60 to 130
	       rrdata["TRP"].push_back( RichardsonRotData( "m95 ",    -195 ,    -65,   95,             11, 19 ) );//   60 to 130
	if(ER) rrdata["TRP"].push_back( RichardsonRotData( "m95 ",    -195 ,    -65,  114,             11, 19 ) );//   60 to 130


			// Tyrosine
	for( float chi3 = 0.0f; chi3 < 181.0f; chi3 += 180.0f){
		if(ER) rrdata["TYR"].push_back( RichardsonRotData( "p90 ",    -182,      62,   73,     chi3,        13, 13, 0 ) );//   60 to 90, -90 to-60
		       rrdata["TYR"].push_back( RichardsonRotData( "p90 ",    -182,      62,   90,     chi3,        13, 13, 0 ) );//   60 to 90, -90 to-60
		if(ER) rrdata["TYR"].push_back( RichardsonRotData( "p90 ",    -182,      62,  107,     chi3,        13, 13, 0 ) );//   60 to 90, -90 to-60

		if(ER) rrdata["TYR"].push_back( RichardsonRotData( "t80 ",    -486,    -177,   63,     chi3,        11, 14, 0 ) );//   20 to 90, -90 to -75
		       rrdata["TYR"].push_back( RichardsonRotData( "t80 ",    -486,    -177,   80,     chi3,        11, 14, 0 ) );//   20 to 90, -90 to -75
		if(ER) rrdata["TYR"].push_back( RichardsonRotData( "t80 ",    -486,    -177,   97,     chi3,        11, 14, 0 ) );//   20 to 90, -90 to -75

		if(ER) rrdata["TYR"].push_back( RichardsonRotData( "m85 ",    -618,     -65, -106,     chi3,        11, 21, 0 ) );//   50 to 90, -90 to -50
		       rrdata["TYR"].push_back( RichardsonRotData( "m85 ",    -618,     -65,  -85,     chi3,        11, 21, 0 ) );//   50 to 90, -90 to -50
		if(ER) rrdata["TYR"].push_back( RichardsonRotData( "m85 ",    -618,     -65,  -64,     chi3,        11, 21, 0 ) );//   50 to 90, -90 to -50

		rrdata["TYR"].push_back( RichardsonRotData( "m30 ",      84,     -65,  -30,     chi3,        11, 18, 0 ) );//   -50 to 0, 0 to 50
		rrdata["TYR"].push_back( RichardsonRotData( "Sm35",      40,     -85,   30,     chi3,        11, 18, 0 ) );//
	}

			// Phenylalanine
	if(ER) rrdata["PHE"].push_back( RichardsonRotData( "p90 ",    -202,      62,   79,             11, 11 ) );//   60 to 90, -90 to-60
	       rrdata["PHE"].push_back( RichardsonRotData( "p90 ",    -202,      62,   90,             11, 11 ) );//   60 to 90, -90 to-60
	if(ER) rrdata["PHE"].push_back( RichardsonRotData( "p90 ",    -202,      62,  101,             11, 11 ) );//   60 to 90, -90 to-60

	if(ER) rrdata["PHE"].push_back( RichardsonRotData( "t80 ",    -522,    -177,   63,             13, 17 ) );//   20 to 90, -90 to -75
	       rrdata["PHE"].push_back( RichardsonRotData( "t80 ",    -522,    -177,   80,             13, 17 ) );//   20 to 90, -90 to -75
	if(ER) rrdata["PHE"].push_back( RichardsonRotData( "t80 ",    -522,    -177,   97,             13, 17 ) );//   20 to 90, -90 to -75

	if(ER) rrdata["PHE"].push_back( RichardsonRotData( "m85 ",    -697,     -65, -102,             12, 17 ) );//   50 to 90, -90 to -50
	       rrdata["PHE"].push_back( RichardsonRotData( "m85 ",    -697,     -65,  -85,             12, 17 ) );//   50 to 90, -90 to -50
	if(ER) rrdata["PHE"].push_back( RichardsonRotData( "m85 ",    -697,     -65,  -68,             12, 17 ) );//   50 to 90, -90 to -50

	if(ER) rrdata["PHE"].push_back( RichardsonRotData( "m30 ",    -100,     -65,  -50,              9, 20 ) );//   -50 to 0, 0 to 50
	       rrdata["PHE"].push_back( RichardsonRotData( "m30 ",    -100,     -65,  -30,              9, 20 ) );//   -50 to 0, 0 to 50
	if(ER) rrdata["PHE"].push_back( RichardsonRotData( "m30 ",    -100,     -65,  -10,              9, 20 ) );//   -50 to 0, 0 to 50

	if(ER) rrdata["PHE"].push_back( RichardsonRotData( "Sm35",     -49,     -85,   10,              9, 20 ) );//
	       rrdata["PHE"].push_back( RichardsonRotData( "Sm35",     -49,     -85,   30,              9, 20 ) );//
	if(ER) rrdata["PHE"].push_back( RichardsonRotData( "Sm35",     -49,     -85,   50,              9, 20 ) );//

			// Proline
	// rrdata["PRO"].push_back( RichardsonRotData( "Cyendo",   379,     30,                    7  ) );//        15 to 60
	// rrdata["PRO"].push_back( RichardsonRotData( "Cyexo ",   372,    -30,                    6  ) );//        -60 to -15
	// rrdata["PRO"].push_back( RichardsonRotData( "cis   ",    56,     30,                    5  ) );//      15 to 60

			// Threonine
	rrdata["THR"].push_back( RichardsonRotData( "p",       1200,     62,                 12/*10*/ ) );
	rrdata["THR"].push_back( RichardsonRotData( "t",        169,   -175,                 12/* 6*/ ) );
	rrdata["THR"].push_back( RichardsonRotData( "m",       1062,    -65,                 12/* 7*/ ) );

			// Valine
	rrdata["VAL"].push_back( RichardsonRotData( "p",        169,     63,                  12/*8*/) );
	rrdata["VAL"].push_back( RichardsonRotData( "t",       1931,    175,                  12/*8*/) );
	rrdata["VAL"].push_back( RichardsonRotData( "m",        526,    -60,                  12/*7*/) );

			// Serine
	rrdata["SER"].push_back( RichardsonRotData( "p",       1201,     62,                  12/*10*/ ) );
	rrdata["SER"].push_back( RichardsonRotData( "t",        541,   -177,                  12/*11*/ ) );
	rrdata["SER"].push_back( RichardsonRotData( "m",        714,    -65,                  12/* 9*/ ) );

			// Cysteine
			// p         64     62                   14
			// t         74   -177                   10
			// m        142    -65                   11

			// Disulfideg range range range
			// mmm       70  -61 -81 -75 -95 to -30 10 to -50 20 to -40
			// ppp       15  63 85 85 30 to 90 55 to 115 60 to 115
			// mpp       33  -65 100 85 00 to -35 70 to 130 50 to 115
			// pmm        6  90 -91 -64 60 to 120 20 to -60 -95 to -35
			// mpm       11  -86 102 -102 20 to -40 70 to 130 60 to -90
			// mmt       19  -92 -90 -149 20 to -30 20 to -60 0 to -120
			// ppt       16  52 82 180 30 to 95 50 to 110 0 to -120
			// mpt        5  -68 96 147 00 to -40 65 to 125 15 to 175
			// tmt        6  172 -83 -168 0 to -160 15 to -55 0 to -140
			// tpt        1  122 87 163 95 to 150 55 to 115 0 to


}

//foward decleration
void
get_rotamer_spec_default(
	::scheme::chemical::RotamerIndexSpec & rot_index,
	bool extra_rotamers,
	bool extra_primary_rotamers
);
std::shared_ptr<RotamerIndex>
get_rotamer_index(
	
	::scheme::chemical::RotamerIndexSpec const& rot_index_spec

){
	std::shared_ptr<RotamerIndex> rot_index = std::make_shared<RotamerIndex>();

	
	//rot_index_spec.fill_rotamer_index(rot_index);
	//std::cout << "finish fill rotamer index now try loading specific rotamers ...  " << std::endl;
	
	//utility::io::ozstream outfile("test_out.text");
	//rot_index_spec.save(outfile);
	//outfile.close();

	//utility::io::izstream infile("test_in.text");
	//rot_index_spec.load(infile);
	//infile.close();
	
	std::cout << "finish builidng rotamer_index now try fill rotamer index ... " << std::endl;
	rot_index_spec.fill_rotamer_index(*rot_index);
	return rot_index;
}
std::shared_ptr<RotamerIndex>
get_rotamer_index(
	std::string cachefile

){
	std::shared_ptr<RotamerIndex> rot_index = std::make_shared<RotamerIndex>();
	::scheme::chemical::RotamerIndexSpec rot_index_spec;
	//load rotamer spec from cache file
	utility::io::izstream infile(cachefile);
	if ( ! infile ) {
		utility_exit_with_message("Rotamer index spec file: " + cachefile +
			" could not be loaded.\nIf you\'re on the digs and you\'re using an old rifgen, add this line to your flags:\n" 
			+ "-rif_dock:rot_spec_fname /home/bcov/sc/random/old_rif_rotamer_index_spec.txt");
	}
	rot_index_spec.load(infile);
	rot_index_spec.fill_rotamer_index(*rot_index);
	return rot_index;
}


void
get_rotamer_spec_default(
	::scheme::chemical::RotamerIndexSpec & rot_index,
	bool extra_rotamers,
	bool extra_primary_rotamers	
){
	std::cout << "get_rotamer_index" << std::endl;

	// todo cache
	std::vector<std::string> resnames;
		resnames.push_back("ALA");
		resnames.push_back("CYS");
		resnames.push_back("ASP");
		resnames.push_back("GLU");
		resnames.push_back("PHE");
		resnames.push_back("GLY");
		resnames.push_back("HIS");
		resnames.push_back("HIS_D");
		resnames.push_back("ILE");
		resnames.push_back("LYS");
		resnames.push_back("LEU");
		resnames.push_back("MET");
		resnames.push_back("ASN");
		resnames.push_back("PRO");
		resnames.push_back("GLN");
		resnames.push_back("ARG");
		resnames.push_back("SER");
		resnames.push_back("THR");
		resnames.push_back("VAL");
		resnames.push_back("TRP");
		resnames.push_back("TYR");

	std::vector<std::string> use_rosetta_rots;
		use_rosetta_rots.push_back("SER");
		use_rosetta_rots.push_back("THR");
		use_rosetta_rots.push_back("CYS");
		use_rosetta_rots.push_back("PRO");


	std::map<std::string,std::vector<RichardsonRotData> > rrdata;
	get_richardson_rot_data( rrdata, extra_primary_rotamers );
	std::map<int,std::pair<std::string,RichardsonRotData const *>> rrd_map;

	rot_index.add_rotamer( "ALA", std::vector<float>(), 0 );
	rot_index.add_rotamer( "GLY", std::vector<float>(), 0 );

	for( auto resname : resnames ){

		core::chemical::ResidueTypeSetCAP rts = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");
		core::chemical::ResidueType const & rtype = rts.lock()->name_map( resname );
		size_t n_proton_chi = rtype.n_proton_chi();

		if( std::count( use_rosetta_rots.begin(), use_rosetta_rots.end(), resname ) ){ // TYR now done through richardson rots

			core::conformation::ResidueOP resop = core::conformation::ResidueFactory::create_residue( rtype );
			core::pose::Pose pose;
			pose.append_residue_by_jump( *resop, 1 );
			// if( resname=="THR" || resname == "ILE" || resname == "VAL" || resname == "TYR" ){ // use sheet
			// 	// parallel 119, +113 / antiparallel -139, +135
			//  	pose.set_phi(1,129);
			//  	pose.set_psi(1,124);
			// } else {
			//  	pose.set_phi(1,-57);
			//  	pose.set_psi(1,-47);
			// }
			core::pack::rotamer_set::RotamerSetOP rotset = get_rosetta_rot_set( pose, 1 );
			int nex = 1, nchi = rotset->rotamer(1)->nchi();
			using namespace basic::options::OptionKeys::packing;
			using namespace basic::options;

				// NO_EXTRA_CHI_SAMPLES          0          original dihedral only; same as using no flag at all
				// EX_ONE_STDDEV                 1 Default  +/- one standard deviation (sd); 3 samples
				// EX_ONE_HALF_STEP_STDDEV       2          +/- 0.5 sd; 3 samples
				// EX_TWO_FULL_STEP_STDDEVS      3          +/- 1 & 2 sd; 5 samples
				// EX_TWO_HALF_STEP_STDDEVS      4          +/- 0.5 & 1 sd; 5 samples
				// EX_FOUR_HALF_STEP_STDDEVS     5          +/- 0.5, 1, 1.5 & 2 sd; 9 samples
				// EX_THREE_THIRD_STEP_STDDEVS   6          +/- 0.33, 0.67, 1 sd; 7 samples
				// EX_SIX_QUARTER_STEP_STDDEVS   7          +/- 0.25, 0.5, 0.75, 1, 1.25 & 1.5 sd; 13 samples",

			if( 1 <= nchi-n_proton_chi && option[ex1::ex1]() && option[ex1::level]()==1 ) nex *= 3;

			if( 2 <= nchi-n_proton_chi && option[ex2::ex2]() && option[ex2::level]()==1 ) nex *= 3;

			if( 3 <= nchi-n_proton_chi && option[ex3::ex3]() && option[ex3::level]()==1 ) nex *= 3;

			if( 4 <= nchi-n_proton_chi && option[ex4::ex4]() && option[ex4::level]()==1 ) nex *= 3;

			// protol chi
			if( resname=="SER" ) nex *= 18;
			if( resname=="TYR" ) nex *=  2;
			if( resname=="THR" ) nex *= 18;
			if( resname=="CYS" ) nex *=  3;

			// std::cerr << option[ex1::ex1]() << " " << option[ex1::level]() << std::endl;
			// std::cerr << option[ex2::ex2]() << " " << option[ex2::level]() << std::endl;
			// std::cerr << option[ex3::ex3]() << " " << option[ex3::level]() << std::endl;
			// std::cerr << option[ex4::ex4]() << " " << option[ex4::level]() << std::endl;
			// std::cerr << nex << std::endl;

			size_t prev_parent_key = -1;
			for(int irot = 1; irot <= rotset->num_rotamers(); ++irot){
				// bool is_primary = (irot-1)%(nex) == 0;
				std::vector<float> chi;
				// if( is_primary ) std::cerr << "PRIMARY ";
				// else             std::cerr << "        ";
				// std::cerr << irot-1 << " " << resname;
				for(int ichi = 1; ichi <= rotset->rotamer(irot)->nchi(); ++ichi){
					chi.push_back( rotset->rotamer(irot)->chi(ichi) );
					// std::cerr << " " << rotset->rotamer(irot)->chi(ichi);
				}
				// std::cerr << " " << prev_parent_key << std::endl;
				// if( is_primary ){
					// prev_parent_key = rot_index.add_rotamer( resname, chi, n_proton_chi );
				// } else {
					// rot_index.add_rotamer( resname, chi, n_proton_chi, prev_parent_key );
				// }
				rot_index.add_rotamer( resname, chi, n_proton_chi );
			}

		} else {

			for( auto const & rrd : rrdata[resname] ){
				int key = rot_index.add_rotamer( resname, rrd.get_chi(), n_proton_chi );
				rrd_map[key] = std::make_pair(resname,&rrd);
			}

		}


	}

	if( extra_rotamers ){
		int const orig_rotindex_size = rot_index.size();
		// std::cerr << "adding extra rotamers up to " << orig_rotindex_size << std::endl;

		// add exrots here...
		for( int irot = 0; irot < orig_rotindex_size; ++irot ){
			if( rrd_map.find(irot) == rrd_map.end() ) continue;
			std::string resname = rot_index.resname(irot);
			int n_proton_chi = 0;
			if( resname=="TYR" ) n_proton_chi = 1;
			bool ex1 = true;
			bool ex2 = true;
			std::string rrdresname = rrd_map[irot].first;
		 	RichardsonRotData const & rrd = *rrd_map[irot].second;
			if( extra_primary_rotamers && rrd.count < 0) ex2 = false; // this is already an ex2 rotamer
			runtime_assert(rot_index.resname(irot) == rrdresname);

			std::vector< std::vector<float> > exchis;
			rrd.get_ex_chis(
			    ex1, //option[ex1::ex1](),
			    ex2, //option[ex2::ex2](),
			    false, //option[ex3::ex3](),
			    false, //option[ex4::ex4](),
			    exchis,
			    false // is_primary
			);
			for( int iex = 0; iex < exchis.size(); ++iex ){
				// std::cerr << irot << " " << resname << " " << iex << std::endl;
				rot_index.add_rotamer( resname, exchis[iex], n_proton_chi, irot /*parent*/ );
			}

		}
	}

	//rot_index.build_index();
	//std::cout << "done building rotamer index, size: " << rot_index.size() << ", nprimary: " << rot_index.n_primary_rotamers() << std::endl;

	// for( std::string resname : resnames ){
	// 	std::cerr << resname << " " <<  rot_index.index_bounds( resname ).second - rot_index.index_bounds( resname ).first << std::endl;
	// }
	// std::cerr << rot_index.size() << std::endl;

	// for( int irot = 0; irot < rot_index.n_primary_rotamers(); ++irot ){
	// 	utility::io::ozstream out( "rotamer_" + rot_index.resname(irot) + str(irot) + ".pdb" );
	// 	rot_index.dump_pdb_with_children(out,irot);
	// 	out.close();
	// }

	// for( int irot = 0; irot < rot_index.size(); ++irot ){
	//  	utility::io::ozstream out( "rotalnstruct_" + rot_index.resname(irot) + str(irot) + ".pdb" );
	//  	rot_index.dump_pdb(out, irot, rot_index.to_structural_parent_frame_.at(irot) );
	//  	out.close();
	// }

	// for( int isp : rot_index.structural_parents_ ){
	// 	utility::io::ozstream out( "rotstructgroup_" + rot_index.resname(isp) + str(isp) + ".pdb" );
	// 	rot_index.dump_pdb_by_structure(out,isp);
	// 	out.close();
	// }

	// utility_exit_with_message("testing out extra rotamers");
}









std::set<std::pair<int,int>>
get_satisfied_atoms(core::pose::Pose pose, float ethresh){
 	std::set<std::pair<int,int>> sat;

 	core::scoring::ScoreFunctionOP sf = core::scoring::get_score_function(true);
	core::scoring::methods::EnergyMethodOptions myopt = sf->energy_method_options();
	myopt.hbond_options().decompose_bb_hb_into_pair_energies(true);
	sf->set_energy_method_options(myopt);
 	sf->score(pose);
	core::scoring::hbonds::HBondSet hbset;
	core::scoring::hbonds::fill_hbond_set(pose,false,hbset);

	for(core::Size ihb = 1; ihb <= hbset.nhbonds(); ++ihb){
		core::scoring::hbonds::HBond const & hb(hbset.hbond(ihb));
		if( hb.energy() <= ethresh ){
			int dr = hb.don_res();
			int dh = hb.don_hatm();
			int ar = hb.acc_res();
			int aa = hb.acc_atm();
			sat.insert(std::make_pair(dh,dr));
			sat.insert(std::make_pair(aa,ar));
		}
	}

	return sat;
}



void get_donor_rays(core::pose::Pose const& pose, int ir, HBRayOpts const & opt, std::vector<HBondRay>& donors )
{
    std::vector<std::pair<int,std::string>> tmp;
	get_donor_rays(pose,ir,opt,donors,tmp);
}

void get_acceptor_rays(core::pose::Pose const& pose, int ir, HBRayOpts const & opt, std::vector<HBondRay>& donors )
{
    std::vector<std::pair<int,std::string>> tmp;
	get_acceptor_rays(pose,ir,opt,donors,tmp);
}


void get_donor_rays( core::pose::Pose const & pose, int ir, HBRayOpts const & opt,
	std::vector<HBondRay> & donors, std::vector<std::pair<int,std::string>> & anames )
{
	core::chemical::AtomIndices const * positions = &pose.residue_type(ir).Hpos_polar_sc();
	if( opt.withbb )                    positions = &pose.residue_type(ir).Hpos_polar();
	for( auto ih : *positions ){
		if( opt.satisfied_atoms.count(std::make_pair(int(ih),int(ir))) ){
			continue; // is satisfied internally
		}
		core::Size ib = pose.residue_type(ir).atom_base(ih);
		HBondRay hray;
		::Eigen::Vector3f base;
		for(int k = 0; k < 3; ++k) hray.horb_cen[k] = pose.residue(ir).xyz(ih)[k];
		for(int k = 0; k < 3; ++k)          base[k] = pose.residue(ir).xyz(ib)[k];
		hray.direction = hray.horb_cen - base;
		hray.direction.normalize();
		// hray.horb_cen += hray.direction * ( 2.7 - 1.01 - 0.61 )/2.0; // move so position is at cen of ideal hbond
		// hray.horb_cen += hray.direction * 0.1;
		donors.push_back( hray );
		anames.push_back( std::make_pair(ir,pose.residue(ir).atom_name(ih)));
	}
}




void get_acceptor_rays_lkball( core::pose::Pose const & pose, int ir, HBRayOpts const & opt,
	std::vector<HBondRay> & acceptors, std::vector<std::pair<int,std::string>> & anames )
{
	core::conformation::Residue const & rsd = pose.residue(ir);
	core::chemical::AtomIndices const * positions = &rsd.accpt_pos_sc();
	if( opt.withbb )                    positions = &rsd.accpt_pos();
	auto lkbinfo = core::scoring::lkball::LKB_ResidueInfo( rsd );
	for( auto iacc : *positions ){

		// if is involdev in N-N bond, skip
		if(rsd.atom_type(iacc).element()=="N"){
			bool isNN = false;
			for( int ia : rsd.bonded_neighbor(iacc) ){
				if( rsd.atom_type(ia).element()=="N" ) isNN = true;
			}
			if(isNN) continue;
		}

		if( opt.satisfied_atoms.count(std::make_pair(int(iacc),int(ir))) )
		{
			continue; // is satisfied internally
		}
		else if( rsd.is_protein() && rsd.atom_is_backbone(iacc) )
		{
			// std::cerr << "IS BB " << ir << " " << iacc << std::endl;
			std::string aname = rsd.atom_name(iacc);
			auto caxyz = rsd.xyz("CA");
			auto cxyz = rsd.xyz("C");
			if( aname==" O  " || (rsd.is_upper_terminus() && aname==" OXT")){
				auto oxyz = rsd.xyz("O");
				if( aname==" OXT" ) oxyz = rsd.xyz("OXT");
				auto const dir1 = (cxyz-caxyz).normalized();
				auto const dir2 = -(dir1+(cxyz-oxyz).normalized()).normalized();
				auto const dir3 = (oxyz-cxyz).normalized();
				auto const cen1 = oxyz + dir1*opt.orblen;
				auto const cen2 = oxyz + dir2*opt.orblen;
				auto const cen3 = oxyz + dir3*opt.orblen;
				HBondRay hr1, hr2, hr3;
				for( int k = 0; k < 3; ++k) hr1.direction[k] = dir1[k];
				for( int k = 0; k < 3; ++k) hr1.horb_cen[k] = cen1[k];
				for( int k = 0; k < 3; ++k) hr2.direction[k] = dir2[k];
				for( int k = 0; k < 3; ++k) hr2.horb_cen[k] = cen2[k];
				for( int k = 0; k < 3; ++k) hr3.direction[k] = dir3[k];
				for( int k = 0; k < 3; ++k) hr3.horb_cen[k] = cen3[k];
				acceptors.push_back(hr1);
				anames.push_back( std::make_pair(ir,rsd.atom_name(iacc)));
				if(opt.add_acceptor_mid){
					acceptors.push_back(hr3);
					anames.push_back( std::make_pair(ir,rsd.atom_name(iacc)));
				}
				acceptors.push_back(hr2);
				anames.push_back( std::make_pair(ir,rsd.atom_name(iacc)));
			} else {
				std::cout << "WARNING: unknown backbone acceptor '" << rsd.atom_name(iacc) << "'" << std::endl;
			}
		}
		else if( rsd.is_DNA() && rsd.atom_is_backbone(iacc) )
		{
			if(rsd.is_terminus()){
				std::cout << "need to implement DNA backbone termini " << ir << " " << rsd.atom_name(iacc) << std::endl;
			} else {
				std::string aname = rsd.atom_name(iacc);
				if(aname==" OP1" || aname==" OP2"){
					auto opxyz = rsd.xyz(iacc);
					auto pxyz = rsd.xyz("P");
					auto oxyz = rsd.xyz("O5'");
					auto op_p = (opxyz-pxyz).normalized();
					auto o_p = oxyz-pxyz;
					o_p -= o_p.dot(op_p)*op_p;
					o_p.normalize();
					runtime_assert( fabs(o_p.dot(op_p)) < 0.001 );
					auto dir0 = 0.5*op_p + sqrt(3.0)/2.0*o_p;
					HBondRay hr0;
					for( int k = 0; k < 3; ++k) hr0.direction[k] = dir0[k];
					for( int k = 0; k < 3; ++k) hr0.horb_cen[k] = (opxyz + dir0*opt.orblen)[k];
					acceptors.push_back(hr0);
					anames.push_back( std::make_pair(ir,rsd.atom_name(iacc)));
					auto rot = numeric::rotation_matrix_degrees(opxyz-pxyz,60.0);
					for(int j = 1; j < 6; ++j){
						dir0 = rot*dir0; // rotate by 60 degrees
						HBondRay hr;
						for( int k = 0; k < 3; ++k) hr.direction[k] = dir0[k];
						for( int k = 0; k < 3; ++k) hr.horb_cen[k] = (opxyz + dir0*opt.orblen)[k];
						acceptors.push_back(hr);
						anames.push_back( std::make_pair(ir,rsd.atom_name(iacc)));
					}
				}
				else if(aname==" O4'") // C1'		O4'		C4'
				{
					auto c1 = rsd.xyz("C1'");
					auto o4 = rsd.xyz("O4'");
					auto c4 = rsd.xyz("C4'");
					auto ocen = (o4 - (c1+c4)/2.0).normalized();
					auto operp = (o4-c1);
					operp -= operp.dot(ocen)*ocen;
					operp.normalize();
					operp = numeric::rotation_matrix_degrees(ocen,90.0)*operp;
					float sn = sin(numeric::conversions::radians((180.0-109.5)/2.0));
					float cs = cos(numeric::conversions::radians((180.0-109.5)/2.0));
					auto dir1 = sn*ocen + cs*operp;
					auto dir2 = sn*ocen - cs*operp;
					HBondRay hr0, hr1;
					for( int k = 0; k < 3; ++k) hr0.direction[k] = dir1[k];
					for( int k = 0; k < 3; ++k) hr0.horb_cen[k] = (o4 + dir1*opt.orblen)[k];
					for( int k = 0; k < 3; ++k) hr1.direction[k] = dir2[k];
					for( int k = 0; k < 3; ++k) hr1.horb_cen[k] = (o4 + dir2*opt.orblen)[k];
					acceptors.push_back(hr0);
					acceptors.push_back(hr1);
					anames.push_back( std::make_pair(ir,rsd.atom_name(iacc)));
					anames.push_back( std::make_pair(ir,rsd.atom_name(iacc)));
				}
				else if( aname==" O5'" || (aname==" O3'" && ir < pose.size() && pose.residue(ir+1).is_DNA()) )
				{
					auto cxyz = rsd.xyz("C5'");
					auto oxyz = rsd.xyz("O5'");
					auto pxyz = rsd.xyz("P");
					if( aname==" O3'" ){
						cxyz = rsd.xyz("C3'");
						oxyz = rsd.xyz("O3'");
						pxyz = pose.residue(ir+1).xyz("P");
					}
					auto dir = ( oxyz - (cxyz+pxyz)/2.0 ).normalized();
					HBondRay hr;
					for( int k = 0; k < 3; ++k) hr.direction[k] = dir[k];
					for( int k = 0; k < 3; ++k) hr.horb_cen[k] = (oxyz + dir*opt.orblen)[k];
					acceptors.push_back(hr);
					anames.push_back( std::make_pair(ir,rsd.atom_name(iacc)));
				}
				else
				{
					std::cout << "WARNING: unknown DNA backbone acceptor '" << rsd.atom_name(iacc) << "'" << std::endl;
				}
			}

		}
		else if( rsd.is_RNA() && rsd.atom_is_backbone(iacc) ){
			std::cout << "need to implement RNA backbone" << std::endl;
		}
		else // use lkball for sidechains
		{
			auto waters = lkbinfo.waters()[iacc];
			auto basexyz = rsd.xyz(iacc);
			std::vector<HBondRay> rays_this_atom;
			for( auto watxyz : waters ){
				bool is_donor_wat = false;
				for( core::Size ih = rsd.attached_H_begin(iacc); ih <= rsd.attached_H_end(iacc); ++ih){
					auto hxyz = rsd.xyz(ih);
					if( (hxyz-basexyz).normalized().dot((watxyz-basexyz).normalized()) > 0.8 ){
						is_donor_wat = true;
					}
				}
				if( is_donor_wat ) continue;
				HBondRay hray;
				::Eigen::Vector3f base, wat;
				for(int k = 0; k < 3; ++k) wat[k] = watxyz[k];
				for(int k = 0; k < 3; ++k) base[k] = basexyz[k];
				hray.direction = wat - base;
				hray.direction.normalize();
				hray.horb_cen = base + hray.direction*opt.orblen; // 0.61 seems to be what the rosetta orbitals use
				rays_this_atom.push_back(hray);
			}
			if( rays_this_atom.size()==2 && opt.add_acceptor_mid ){
				HBondRay hray, ray1=rays_this_atom[0], ray2=rays_this_atom[1];
				hray.direction = (ray1.direction + ray2.direction).normalized();
				hray.horb_cen = ray1.horb_cen - ray1.direction*opt.orblen + hray.direction*opt.orblen;
				acceptors.push_back(rays_this_atom[0]);
				acceptors.push_back(hray);
				acceptors.push_back(rays_this_atom[1]);
				anames.push_back( std::make_pair(ir,rsd.atom_name(iacc)));
				anames.push_back( std::make_pair(ir,rsd.atom_name(iacc)));
				anames.push_back( std::make_pair(ir,rsd.atom_name(iacc)));
			} else {
				for(auto hray: rays_this_atom){
					acceptors.push_back(hray);
					anames.push_back( std::make_pair(ir,rsd.atom_name(iacc)));
				}
			}
			// std::cout << ir << " " << rsd.name3() << " " << rsd.atom_name(iacc) << " " << rays_this_atom.size() << std::endl;
		}
	}

}

void get_acceptor_rays_rosetta_orbitals( core::pose::Pose const & pose, int ir, HBRayOpts const & opt,
	std::vector<HBondRay> & acceptors, std::vector<std::pair<int,std::string>> & anames )
{
	core::conformation::Residue const & rsd = pose.residue(ir);
	core::chemical::AtomIndices const * positions = &rsd.accpt_pos_sc();
	if( opt.withbb )                    positions = &rsd.accpt_pos();
	std::vector< ::Eigen::Vector3f> seenit;
	for( auto iacc : *positions ){
		if( opt.satisfied_atoms.count(std::make_pair(int(iacc),int(ir))) ){
			continue; // is satisfied internally
		}
		for( auto jorb : rsd.bonded_orbitals( iacc ) ){
			if( // rsd.orbital_type(jorb).orbital_enum() == core::chemical::orbitals::C_pi_sp2    ||
			    // rsd.orbital_type(jorb).orbital_enum() == core::chemical::orbitals::N_pi_sp2    ||
			    rsd.orbital_type(jorb).orbital_enum() == core::chemical::orbitals::N_p_sp2     ||
			    // rsd.orbital_type(jorb).orbital_enum() == core::chemical::orbitals::O_pi_sp2    ||
			    rsd.orbital_type(jorb).orbital_enum() == core::chemical::orbitals::O_p_sp2     ||
			    // rsd.orbital_type(jorb).orbital_enum() == core::chemical::orbitals::S_p_sp3     ||
			    // rsd.orbital_type(jorb).orbital_enum() == core::chemical::orbitals::O_pi_sp2_bb ||
			    // rsd.orbital_type(jorb).orbital_enum() == core::chemical::orbitals::O_p_sp2_bb  ||
			    rsd.orbital_type(jorb).orbital_enum() == core::chemical::orbitals::O_p_sp3
			){
				HBondRay hray;
				::Eigen::Vector3f base;
				for(int k = 0; k < 3; ++k) hray.horb_cen[k] = rsd.orbital_xyz(jorb)[k];
				for(int k = 0; k < 3; ++k)          base[k] = rsd.        xyz(iacc)[k];
				// sometimes rosetta generates duplicate orbitals???
				bool isnew = true;
				for( int i = 0; i < seenit.size(); ++i ){
					if( (seenit[i]-hray.horb_cen).squaredNorm() < 0.001 ) isnew = false;
				}
				if( !isnew ){
					std::cout << "DUPLICATE ORBITAL! WHY DOES ROSETTA DO THIS??? " << rsd.name3() << " "
				              << rsd.atom_name(iacc) << " orb " << rsd.orbital_name(jorb)
				              << " orb# = " << jorb << std::endl;
				} else {
					seenit.push_back( hray.horb_cen );
					hray.direction = hray.horb_cen - base;
					hray.direction.normalize();
					// hray.horb_cen += hray.direction * ( 2.7 - 1.01 - 0.61 )/2.0; // move so position is at cen of ideal hbond
					// hray.horb_cen += hray.direction * 0.1;
					// std::cout << "ACCEPTOR " << rsd.atom_name(iacc) << std::endl;
					acceptors.push_back(hray);
					anames.push_back( std::make_pair(ir,rsd.atom_name(iacc)));
				}
			}
		}
	}

}

void get_acceptor_rays( core::pose::Pose const & pose, int ir, HBRayOpts const & opt,
	std::vector<HBondRay> & acceptors, std::vector<std::pair<int,std::string>> & anames ){
		if(opt.lkball) get_acceptor_rays_lkball(pose, ir, opt, acceptors, anames);
		else get_acceptor_rays_rosetta_orbitals(pose, ir, opt, acceptors, anames);
}


void dump_hbond_rays( std::ostream & out, std::vector<HBondRay> hbonders, bool isdonor ){
	int anum=0, rnum=0;
	for( auto hbr : hbonders ){
		for( int ibase = 0; ibase < 2; ibase++ ){
			float x = hbr.horb_cen[0] + ( ibase? 0.0f : hbr.direction[0] );
			float y = hbr.horb_cen[1] + ( ibase? 0.0f : hbr.direction[1] );
			float z = hbr.horb_cen[2] + ( ibase? 0.0f : hbr.direction[2] );
			char buf[128];
			std::string aname;
			if( isdonor ) aname = ibase? "DBSE" : "DHPL";
			else          aname = ibase? "ABSE" : "AORB";
			snprintf(buf,128,"%s%5i %4s %3s %c%4i    %8.3f%8.3f%8.3f%6.2f%6.2f %11s\n",
				"HETATM",
				++anum,
				aname.c_str(),
				isdonor? "DON" : "ACC",
				isdonor? "D" : "A",
				rnum,
				x,y,z,
				1.0,
				1.0,
				"HB"
			);
			out << buf;
		}
		++rnum;
	}
}




}}

