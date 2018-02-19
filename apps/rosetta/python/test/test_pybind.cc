#include <pybind11/pybind11.h>
#include <memory>
#include <iostream>

namespace py = pybind11;

struct Foo : public std::enable_shared_from_this<Foo> {
	Foo(){}
	void bar() const {
		std::cout << "bar" << std::endl;
	}
};

std::shared_ptr<Foo> make_Foo(){
	return std::make_shared<Foo>();
}

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

PYBIND11_PLUGIN(test_pybind) {
    py::module m("test_pybind", "test pybind11");

	py::class_< Foo, std::shared_ptr<Foo> >(m, "Foo")
		.def( py::init<>() )
		.def( "bar", &Foo::bar )
	;

 	m.def( "make_Foo", &make_Foo );

    return m.ptr();
}

