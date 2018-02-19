# #### devtoolset-3:
# COLLECT_GCC=/opt/rh/devtoolset-3/root/usr/bin/gcc
# COLLECT_LTO_WRAPPER=/opt/rh/devtoolset-3/root/usr/libexec/gcc/x86_64-redhat-linux/4.9.2/lto-wrapper
# Target: x86_64-redhat-linux
# Configured with: ../configure --prefix=/opt/rh/devtoolset-3/root/usr --mandir=/opt/rh/devtoolset-3/root/usr/share/man --infodir=/opt/rh/devtoolset-3/root/usr/share/info --with-bugurl=http://bugzilla.redhat.com/bugzilla --enable-bootstrap --enable-shared --enable-threads=posix --enable-checking=release --enable-multilib --with-system-zlib --enable-__cxa_atexit --disable-libunwind-exceptions --enable-gnu-unique-object --enable-linker-build-id --enable-languages=c,c++,fortran,lto --enable-plugin --with-linker-hash-style=gnu --enable-initfini-array --disable-libgcj --with-isl=/builddir/build/BUILD/gcc-4.9.2-20150212/obj-x86_64-redhat-linux/isl-install --with-cloog=/builddir/build/BUILD/gcc-4.9.2-20150212/obj-x86_64-redhat-linux/cloog-install --enable-gnu-indirect-function --with-tune=generic --with-arch_32=i686 --build=x86_64-redhat-linux
# Thread model: posix
# gcc version 4.9.2 20150212 (Red Hat 4.9.2-6) (GCC)

# #### my gcc:
# COLLECT_GCC=g++
# COLLECT_LTO_WRAPPER=/scratch/USERS/sheffler/root/gcc61/libexec/gcc/x86_64-pc-linux-gnu/6.1.0/lto-wrapper
# Target: x86_64-pc-linux-gnu
# Configured with: ../gcc-6.1.0/configure --prefix=/scratch/USERS/sheffler/root/gcc61 --enable-languages=c,c++,fortran --disable-multilib
# Thread model: posix
# gcc version 6.1.0 (GCC)

cmake_minimum_required (VERSION 3.3)

set( PYTHON_VERSION $ENV{PYVERSION} )
set( CXX_STD 11 )
message( "build type .................... " ${CMAKE_BUILD_TYPE} )
message( "python version ................ " ${PYTHON_VERSION} )
message( "c++ standard .................. " "c++${CXX_STD}" )


if( ${CXX_STD} MATCHES 11 )
  set( CXX_DIR cxx11_omp )
  add_definitions( -DCXX11 )
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -w -fno-omit-frame-pointer -Wall -Wextra -Werror" )
elseif( ${CXX_STD} MATCHES 14 )
	set( CXX_DIR cxx14_omp )
  add_definitions( -DCXX11 )
  add_definitions( -DCXX14 )
	set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -w -fno-omit-frame-pointer -Wall -Wextra -Werror" )
else()
	message( FATAL_ERROR "bad CXX_STD ${CXX_STD}" )
endif()

string(TOUPPER "${CMAKE_BUILD_TYPE}" U_CMAKE_BUILD_TYPE)
if( ${U_CMAKE_BUILD_TYPE} MATCHES DEBUG )

    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Og" )
    # set( ROSETTA_BIN $ENV{CMAKE_ROSETTA_PATH}/source/cmake/build_${CXX_DIR}_omp_debug )
    set( ROSETTA_BIN $ENV{CMAKE_ROSETTA_PATH}/source/cmake/build_${CXX_DIR}_debug )
    if (NOT $ENV{CMAKE_FINAL_ROSETTA_PATH} STREQUAL "")
        set( ROSETTA_BIN $ENV{CMAKE_FINAL_ROSETTA_PATH} )
    endif()

else()
    # set( ROSETTA_BIN $ENV{CMAKE_ROSETTA_PATH}/source/cmake/build_${CXX_DIR}_omp )
    set( ROSETTA_BIN $ENV{CMAKE_ROSETTA_PATH}/source/cmake/build_${CXX_DIR} )
    if (NOT $ENV{CMAKE_FINAL_ROSETTA_PATH} STREQUAL "")
        set( ROSETTA_BIN $ENV{CMAKE_FINAL_ROSETTA_PATH} )
    endif()
endif()
message( "path to rosetta binaries ...... " ${ROSETTA_BIN} )



set( CMAKE_INSTALL_PREFIX $ENV{SCHEME_INSTALL_PATH} )
message( "install prefix ................ " ${CMAKE_INSTALL_PREFIX} )


######################### external dependencies ####################################

#find_path( EIGEN3_INCLUDE_DIR eigen3/Eigen )
#if( ${EIGEN3_INCLUDE_DIR} MATCHES EIGEN3_INCLUDE_DIR-NOTFOUND )
#	message( FATAL_ERROR "NO eigen3/Eigen found")
#endif()
set( EIGEN3_INCLUDE_DIR external )
include_directories( ${EIGEN3_INCLUDE_DIR}/eigen3 )

add_subdirectory(external/gmock)
include_directories("external/gmock/gtest/include")
include_directories("external/gmock/include")
list( APPEND EXTRA_LIBS gmock )

include_directories( "external/include" )
include_directories( "schemelib" )
include_directories( "/usr/local/include" )
link_directories( "/usr/local/lib" )
link_directories( "/usr/local/lib64" )

####################################################################################

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")

# # remove sublime text turds
# message( "removing ._* turds from " ${CMAKE_SOURCE_DIR} )
# execute_process( COMMAND find ${CMAKE_SOURCE_DIR} -regex .*[.]_.* -exec rm {} \; )


if ( ${USEHDF5} )
    add_definitions("-DUSEHDF5")
endif()



add_subdirectory(schemelib)
add_subdirectory(apps)

