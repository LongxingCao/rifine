#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <devel/init.hh>

namespace py = pybind11;

void rosetta_init( std::vector<std::string> args ) {

	// copied from core::init::init
	std::cout << "options from python:" << std::endl;
	int argc = args.size();
	char **argv = new char*[ argc ];
	for ( int ii = 0; ii < argc; ++ii ) {
		argv[ ii ] = new char[ args[ii].size()+1 ];
		strncpy( argv[ii], args[ii].c_str(), args[ii].size() );
		argv[ ii ][ args[ii].size() ] = 0; // ensure null termination
		std::cout << "    '" << args[ii] << "'" << std::endl;
		// std::cout << "(c) '" << argv[ii] << "'" << std::endl;
	}
	::devel::init( argc, argv );

	// std::cout << "initializing rosetta with args:" << std::endl;
	// char** c_args = new char*[args.size()];
	// int i = 0;
	// for( auto const & arg : args ){
	// 	// if( arg[0] == '-' ) std::cout << std::endl << "    ";
	// 	c_args[i] = new char[arg.size()+1];
	// 	for( int j = 0; j < arg.size(); ++j ){
	// 		c_args[i][j] = arg[j];
	// 	}
	// 	c_args[i][arg.size()] = 0;
	// 	std::cout << " in: '" <<   arg     << "' " << std::endl;
	// 	std::cout << "  c: '" << c_args[i] << "' " << std::endl;		
	// 	++i;
	// }
	// std::cout << std::endl;

}

PYBIND11_PLUGIN(_pysetta_devel) {
    py::module m("_pysetta_devel", "rosetta devel namespace");
    m.def("init", &rosetta_init, "init rosetta from 'cli args'");
    return m.ptr();
}

