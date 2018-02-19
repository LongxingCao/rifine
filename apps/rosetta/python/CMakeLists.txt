
project(pysetta)

include(CheckCXXCompilerFlag)

# Set a default build configuration if none is specified.
# 'MinSizeRel' produces the smallest binaries
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'MinSizeRel' as none was specified.")
  set(CMAKE_BUILD_TYPE MinSizeRel CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()
string(TOUPPER "${CMAKE_BUILD_TYPE}" U_CMAKE_BUILD_TYPE)


if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  # Enable link time optimization and set the default symbol
  # visibility to hidden (very important to obtain small binaries)
  if (NOT ${U_CMAKE_BUILD_TYPE} MATCHES DEBUG)
    # Default symbol visibility
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility=hidden")

    # Check for Link Time Optimization support
    CHECK_CXX_COMPILER_FLAG("-flto" HAS_LTO_FLAG)
    if (HAS_LTO_FLAG)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flto")
    endif()
  endif()
endif()

# Include path for Python header files
include_directories(../../../external/pybind11/include)

#####################################################################################

if(${PYTHON_VERSION} MATCHES 2)
  find_path( PYTHON_INCLUDE_DIR python2.7/Python.h )
  find_library( PYTHON_LIB_FILE python2.7 )
  get_filename_component( PYTHON_LIB_DIR ${PYTHON_LIB_FILE} DIRECTORY )
  set( PYTHON_LIB_NAME python2.7 )
  include_directories( ${PYTHON_INCLUDE_DIR}/${PYTHON_LIB_NAME} )
  link_directories( ${PYTHON_LIB_DIR} )
  install(DIRECTORY pysetta/ DESTINATION lib/python2/pysetta FILES_MATCHING PATTERN "*.py")
else()
  find_path( PYTHON_INCLUDE_DIR python3.5m/Python.h )
  find_library( PYTHON_LIB_FILE python3.5m )
  get_filename_component( PYTHON_LIB_DIR ${PYTHON_LIB_FILE} DIRECTORY )
  set( PYTHON_LIB_NAME python3.5m )  
  include_directories( ${PYTHON_INCLUDE_DIR}/${PYTHON_LIB_NAME} )
  link_directories( ${PYTHON_LIB_DIR} )
  install(DIRECTORY pysetta/ DESTINATION lib/python3/pysetta FILES_MATCHING PATTERN "*.py")
endif()
message( "using python include .......... " ${PYTHON_INCLUDE_DIR}/${PYTHON_LIB_NAME} )
message( "using python lib dir .......... " ${PYTHON_LIB_DIR} )
message( "using python lib .............. " -l${PYTHON_LIB_NAME} )
message( "install python libs to ........ " ${CMAKE_INSTALL_PREFIX}/lib/python${PYTHON_VERSION} )




add_subdirectory(pysetta)
add_subdirectory(test)

