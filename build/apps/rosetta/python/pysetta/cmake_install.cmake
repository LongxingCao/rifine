# Install script for directory: /home/longxing/rosetta_gitclone/rifine/apps/rosetta/python/pysetta

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python/_pysetta_core_import_pose.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python/_pysetta_core_import_pose.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python/_pysetta_core_import_pose.so"
         RPATH "/usr/local/lib:/usr/local/lib64:/software/rosetta-master/source/cmake/build_cxx11_omp")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/python" TYPE SHARED_LIBRARY FILES "/home/longxing/rosetta_gitclone/rifine/build/apps/rosetta/python/pysetta/_pysetta_core_import_pose.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python/_pysetta_core_import_pose.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python/_pysetta_core_import_pose.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python/_pysetta_core_import_pose.so"
         OLD_RPATH "/usr/local/lib:/usr/local/lib64:/software/rosetta-master/source/cmake/build_cxx11_omp:"
         NEW_RPATH "/usr/local/lib:/usr/local/lib64:/software/rosetta-master/source/cmake/build_cxx11_omp")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python/_pysetta_core_import_pose.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python/_pysetta_core_pose.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python/_pysetta_core_pose.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python/_pysetta_core_pose.so"
         RPATH "/usr/local/lib:/usr/local/lib64:/software/rosetta-master/source/cmake/build_cxx11_omp")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/python" TYPE SHARED_LIBRARY FILES "/home/longxing/rosetta_gitclone/rifine/build/apps/rosetta/python/pysetta/_pysetta_core_pose.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python/_pysetta_core_pose.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python/_pysetta_core_pose.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python/_pysetta_core_pose.so"
         OLD_RPATH "/usr/local/lib:/usr/local/lib64:/software/rosetta-master/source/cmake/build_cxx11_omp:"
         NEW_RPATH "/usr/local/lib:/usr/local/lib64:/software/rosetta-master/source/cmake/build_cxx11_omp")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python/_pysetta_core_pose.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python/_pysetta_devel.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python/_pysetta_devel.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python/_pysetta_devel.so"
         RPATH "/usr/local/lib:/usr/local/lib64:/software/rosetta-master/source/cmake/build_cxx11_omp")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/python" TYPE SHARED_LIBRARY FILES "/home/longxing/rosetta_gitclone/rifine/build/apps/rosetta/python/pysetta/_pysetta_devel.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python/_pysetta_devel.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python/_pysetta_devel.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python/_pysetta_devel.so"
         OLD_RPATH "/usr/local/lib:/usr/local/lib64:/software/rosetta-master/source/cmake/build_cxx11_omp:"
         NEW_RPATH "/usr/local/lib:/usr/local/lib64:/software/rosetta-master/source/cmake/build_cxx11_omp")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python/_pysetta_devel.so")
    endif()
  endif()
endif()

