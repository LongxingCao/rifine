# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/longxing/rosetta_gitclone/rifine

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/longxing/rosetta_gitclone/rifine/build

# Include any dependencies generated for this target.
include apps/rosetta/python/pysetta/CMakeFiles/_pysetta_devel.dir/depend.make

# Include the progress variables for this target.
include apps/rosetta/python/pysetta/CMakeFiles/_pysetta_devel.dir/progress.make

# Include the compile flags for this target's objects.
include apps/rosetta/python/pysetta/CMakeFiles/_pysetta_devel.dir/flags.make

apps/rosetta/python/pysetta/CMakeFiles/_pysetta_devel.dir/_pysetta_devel.cc.o: apps/rosetta/python/pysetta/CMakeFiles/_pysetta_devel.dir/flags.make
apps/rosetta/python/pysetta/CMakeFiles/_pysetta_devel.dir/_pysetta_devel.cc.o: ../apps/rosetta/python/pysetta/_pysetta_devel.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/longxing/rosetta_gitclone/rifine/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object apps/rosetta/python/pysetta/CMakeFiles/_pysetta_devel.dir/_pysetta_devel.cc.o"
	cd /home/longxing/rosetta_gitclone/rifine/build/apps/rosetta/python/pysetta && /usr/local/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/_pysetta_devel.dir/_pysetta_devel.cc.o -c /home/longxing/rosetta_gitclone/rifine/apps/rosetta/python/pysetta/_pysetta_devel.cc

apps/rosetta/python/pysetta/CMakeFiles/_pysetta_devel.dir/_pysetta_devel.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/_pysetta_devel.dir/_pysetta_devel.cc.i"
	cd /home/longxing/rosetta_gitclone/rifine/build/apps/rosetta/python/pysetta && /usr/local/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/longxing/rosetta_gitclone/rifine/apps/rosetta/python/pysetta/_pysetta_devel.cc > CMakeFiles/_pysetta_devel.dir/_pysetta_devel.cc.i

apps/rosetta/python/pysetta/CMakeFiles/_pysetta_devel.dir/_pysetta_devel.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/_pysetta_devel.dir/_pysetta_devel.cc.s"
	cd /home/longxing/rosetta_gitclone/rifine/build/apps/rosetta/python/pysetta && /usr/local/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/longxing/rosetta_gitclone/rifine/apps/rosetta/python/pysetta/_pysetta_devel.cc -o CMakeFiles/_pysetta_devel.dir/_pysetta_devel.cc.s

apps/rosetta/python/pysetta/CMakeFiles/_pysetta_devel.dir/_pysetta_devel.cc.o.requires:

.PHONY : apps/rosetta/python/pysetta/CMakeFiles/_pysetta_devel.dir/_pysetta_devel.cc.o.requires

apps/rosetta/python/pysetta/CMakeFiles/_pysetta_devel.dir/_pysetta_devel.cc.o.provides: apps/rosetta/python/pysetta/CMakeFiles/_pysetta_devel.dir/_pysetta_devel.cc.o.requires
	$(MAKE) -f apps/rosetta/python/pysetta/CMakeFiles/_pysetta_devel.dir/build.make apps/rosetta/python/pysetta/CMakeFiles/_pysetta_devel.dir/_pysetta_devel.cc.o.provides.build
.PHONY : apps/rosetta/python/pysetta/CMakeFiles/_pysetta_devel.dir/_pysetta_devel.cc.o.provides

apps/rosetta/python/pysetta/CMakeFiles/_pysetta_devel.dir/_pysetta_devel.cc.o.provides.build: apps/rosetta/python/pysetta/CMakeFiles/_pysetta_devel.dir/_pysetta_devel.cc.o


# Object files for target _pysetta_devel
_pysetta_devel_OBJECTS = \
"CMakeFiles/_pysetta_devel.dir/_pysetta_devel.cc.o"

# External object files for target _pysetta_devel
_pysetta_devel_EXTERNAL_OBJECTS =

apps/rosetta/python/pysetta/_pysetta_devel.so: apps/rosetta/python/pysetta/CMakeFiles/_pysetta_devel.dir/_pysetta_devel.cc.o
apps/rosetta/python/pysetta/_pysetta_devel.so: apps/rosetta/python/pysetta/CMakeFiles/_pysetta_devel.dir/build.make
apps/rosetta/python/pysetta/_pysetta_devel.so: apps/rosetta/python/pysetta/CMakeFiles/_pysetta_devel.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/longxing/rosetta_gitclone/rifine/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX shared library _pysetta_devel.so"
	cd /home/longxing/rosetta_gitclone/rifine/build/apps/rosetta/python/pysetta && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/_pysetta_devel.dir/link.txt --verbose=$(VERBOSE)
	cd /home/longxing/rosetta_gitclone/rifine/build/apps/rosetta/python/pysetta && strip /home/longxing/rosetta_gitclone/rifine/build/apps/rosetta/python/pysetta/_pysetta_devel.so

# Rule to build all files generated by this target.
apps/rosetta/python/pysetta/CMakeFiles/_pysetta_devel.dir/build: apps/rosetta/python/pysetta/_pysetta_devel.so

.PHONY : apps/rosetta/python/pysetta/CMakeFiles/_pysetta_devel.dir/build

apps/rosetta/python/pysetta/CMakeFiles/_pysetta_devel.dir/requires: apps/rosetta/python/pysetta/CMakeFiles/_pysetta_devel.dir/_pysetta_devel.cc.o.requires

.PHONY : apps/rosetta/python/pysetta/CMakeFiles/_pysetta_devel.dir/requires

apps/rosetta/python/pysetta/CMakeFiles/_pysetta_devel.dir/clean:
	cd /home/longxing/rosetta_gitclone/rifine/build/apps/rosetta/python/pysetta && $(CMAKE_COMMAND) -P CMakeFiles/_pysetta_devel.dir/cmake_clean.cmake
.PHONY : apps/rosetta/python/pysetta/CMakeFiles/_pysetta_devel.dir/clean

apps/rosetta/python/pysetta/CMakeFiles/_pysetta_devel.dir/depend:
	cd /home/longxing/rosetta_gitclone/rifine/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/longxing/rosetta_gitclone/rifine /home/longxing/rosetta_gitclone/rifine/apps/rosetta/python/pysetta /home/longxing/rosetta_gitclone/rifine/build /home/longxing/rosetta_gitclone/rifine/build/apps/rosetta/python/pysetta /home/longxing/rosetta_gitclone/rifine/build/apps/rosetta/python/pysetta/CMakeFiles/_pysetta_devel.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : apps/rosetta/python/pysetta/CMakeFiles/_pysetta_devel.dir/depend

