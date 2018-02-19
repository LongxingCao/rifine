
  add_library(test_pybind SHARED test_pybind.cc )
  target_link_libraries( test_pybind ${PYTHON_LIB_NAME} ${ALL_ROSETTA_LIBS} )
  set_target_properties( test_pybind  PROPERTIES PREFIX "")
  install ( TARGETS test_pybind LIBRARY DESTINATION lib/python2 )

  set_target_properties(test_pybind  PROPERTIES SUFFIX ".so")
