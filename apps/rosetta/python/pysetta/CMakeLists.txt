# Create the binding library

FILE(GLOB PYSETTA_SOURCES "./_pysetta_*.cc")

foreach( SRC_PATH ${PYSETTA_SOURCES} )

  get_filename_component( PYSETTA_SRC ${SRC_PATH} NAME_WE )

  message( "process pysetta source file: " ${PYSETTA_SRC}.cc )

  add_library(${PYSETTA_SRC} SHARED ${PYSETTA_SRC}.cc )

  target_link_libraries( ${PYSETTA_SRC} ${PYTHON_LIB_NAME} ${ALL_ROSETTA_LIBS} )

  set_target_properties( ${PYSETTA_SRC}  PROPERTIES PREFIX "")
  install ( TARGETS ${PYSETTA_SRC} LIBRARY DESTINATION lib/python${PYTHON_VERSION} )

  set_target_properties(${PYSETTA_SRC}  PROPERTIES SUFFIX ".so")

  # Strip unnecessary sections of the binary on Linux/Mac OS
  if(APPLE)
    # message("using APPLE settings")
    set_target_properties(${PYSETTA_SRC}  PROPERTIES MACOSX_RPATH ".")
    set_target_properties(${PYSETTA_SRC}  PROPERTIES LINK_FLAGS "-undefined dynamic_lookup ")
    if (NOT ${U_CMAKE_BUILD_TYPE} MATCHES DEBUG)
      add_custom_command(TARGET ${PYSETTA_SRC}  POST_BUILD
                         COMMAND strip -u -r ${PROJECT_BINARY_DIR}/pysetta/${PYSETTA_SRC}.so)
    endif()
  else()
    if (NOT ${U_CMAKE_BUILD_TYPE} MATCHES DEBUG)
      add_custom_command(TARGET ${PYSETTA_SRC}  POST_BUILD
                         COMMAND strip ${PROJECT_BINARY_DIR}/pysetta/${PYSETTA_SRC}.so)
    endif()
  endif()

endforeach( SRC_PATH )
