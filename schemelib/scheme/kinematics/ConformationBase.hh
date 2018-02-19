#ifndef INCLUDED_kinematics_Conformation_HH
#define INCLUDED_kinematics_Conformation_HH

#include "scheme/types.hh"

#ifdef CEREAL
#include <cereal/access.hpp>
#include <cereal/types/memory.hpp>
#endif


#include <vector>

namespace scheme {
namespace kinematics {


	struct ConformationBase {
		ConformationBase(){}

		// This is only here to make the class polymorphic
		virtual ~ConformationBase() {};
	};



}
}

#endif
