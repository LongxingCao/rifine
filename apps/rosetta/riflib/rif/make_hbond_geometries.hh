// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.



#ifndef INCLUDED_riflib_rif_make_hbond_geometries_hh
#define INCLUDED_riflib_rif_make_hbond_geometries_hh

#include <riflib/util.hh>
#include <riflib/RotamerGenerator.hh>
#include <utility/vector1.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <numeric/xyzTransform.hh>
#include <map>



namespace devel {
namespace scheme {
namespace rif {


struct RelRotPos {
	numeric::xyzVector<float> n, ca, c, stub1, stub2, stub3; // 6*3*4/8 = 9 words
	float score;
	int16_t rotamer;
	int16_t hbonding_atom_num;
	int16_t hbonding_atom_num_other;
	bool operator==(RelRotPos const & o) const {
		return n==o.n && ca==o.ca && c==o.c && stub1==o.stub1 && stub2==o.stub2 && stub3==o.stub3 && score==o.score && rotamer==o.rotamer &&
					hbonding_atom_num==o.hbonding_atom_num && hbonding_atom_num_other==o.hbonding_atom_num_other;
	}
}; // size 10 words = 80 bytes


struct MakeHbondGeomOpts {
	double tip_tol_deg, rot_samp_resl, rot_samp_range;
};

void make_hbond_geometries(
	RotamerIndex const & rot_index,
	std::string resn1,
	std::string resn2,
	bool fix_donor,
	bool fix_acceptor,
	std::map<std::string,core::pose::Pose> const & exemplars,
	utility::vector1< RelRotPos > & hbonders,
	MakeHbondGeomOpts const & opts
 );



}
}
}

#endif
