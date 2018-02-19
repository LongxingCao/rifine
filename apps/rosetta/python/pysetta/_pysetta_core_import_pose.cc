#include <pybind11/pybind11.h>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

using namespace ::core;
using namespace ::core::import_pose;

namespace py = pybind11;


std::shared_ptr<pose::Pose> my_pose_from_file( std::string filename ){
	return pose_from_file( filename );
}

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

PYBIND11_PLUGIN(_pysetta_core_import_pose) {
    py::module m( "_pysetta_core_import_pose", "rosetta core::import_pose" );

    m.def( "pose_from_file", &my_pose_from_file, "read a pose" );

    return m.ptr();
}

