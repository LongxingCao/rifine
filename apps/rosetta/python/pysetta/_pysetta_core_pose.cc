#include <pybind11/pybind11.h>

#include <core/pose/Pose.hh>

using namespace ::core::pose;

namespace py = pybind11;

struct TEST : public std::enable_shared_from_this<TEST> {
	TEST(){}
	void do_something() const {
		std::cout << "something" << std::endl;
	}
};


PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

PYBIND11_PLUGIN(_pysetta_core_pose) {
    py::module m("_pysetta_core_pose", "rosetta core::pose");

	py::class_< TEST, std::shared_ptr<TEST> >(m, "TEST")
		.def( py::init<>() )
		.def( "do_something", &TEST::do_something )
	;

	py::class_< Pose, std::shared_ptr<Pose> >(m, "Pose")
		.def( py::init<>() )
		.def( "dump_pdb", (bool (Pose::*)( std::string const & file_name, std::string const & tag) const)
			 &Pose::dump_pdb, "dump Pose to pdb file", py::arg("file_name"), py::arg("tag") = "1" )
	;



    return m.ptr();
}

/*
import os
from pysetta import devel
from pysetta.core.pose import Pose
from pysetta.core.import_pose import pose_from_file

args = "dummy -database "+os.environ['CMAKE_ROSETTA_PATH']+"/database"
devel.init( args.split() )
p = Pose()
p.dump_pdb("test.pdb")
p = pose_from_file("/work/sheffler/1ffw_native.pdb")
p.dump_pdb("test2.pdb")
*/