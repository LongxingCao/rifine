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
include apps/rosetta/CMakeFiles/rif_dock_test.dir/depend.make

# Include the progress variables for this target.
include apps/rosetta/CMakeFiles/rif_dock_test.dir/progress.make

# Include the compile flags for this target's objects.
include apps/rosetta/CMakeFiles/rif_dock_test.dir/flags.make

apps/rosetta/CMakeFiles/rif_dock_test.dir/rif_dock_test.cc.o: apps/rosetta/CMakeFiles/rif_dock_test.dir/flags.make
apps/rosetta/CMakeFiles/rif_dock_test.dir/rif_dock_test.cc.o: ../apps/rosetta/rif_dock_test.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/longxing/rosetta_gitclone/rifine/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object apps/rosetta/CMakeFiles/rif_dock_test.dir/rif_dock_test.cc.o"
	cd /home/longxing/rosetta_gitclone/rifine/build/apps/rosetta && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/rif_dock_test.dir/rif_dock_test.cc.o -c /home/longxing/rosetta_gitclone/rifine/apps/rosetta/rif_dock_test.cc

apps/rosetta/CMakeFiles/rif_dock_test.dir/rif_dock_test.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/rif_dock_test.dir/rif_dock_test.cc.i"
	cd /home/longxing/rosetta_gitclone/rifine/build/apps/rosetta && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/longxing/rosetta_gitclone/rifine/apps/rosetta/rif_dock_test.cc > CMakeFiles/rif_dock_test.dir/rif_dock_test.cc.i

apps/rosetta/CMakeFiles/rif_dock_test.dir/rif_dock_test.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/rif_dock_test.dir/rif_dock_test.cc.s"
	cd /home/longxing/rosetta_gitclone/rifine/build/apps/rosetta && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/longxing/rosetta_gitclone/rifine/apps/rosetta/rif_dock_test.cc -o CMakeFiles/rif_dock_test.dir/rif_dock_test.cc.s

apps/rosetta/CMakeFiles/rif_dock_test.dir/rif_dock_test.cc.o.requires:

.PHONY : apps/rosetta/CMakeFiles/rif_dock_test.dir/rif_dock_test.cc.o.requires

apps/rosetta/CMakeFiles/rif_dock_test.dir/rif_dock_test.cc.o.provides: apps/rosetta/CMakeFiles/rif_dock_test.dir/rif_dock_test.cc.o.requires
	$(MAKE) -f apps/rosetta/CMakeFiles/rif_dock_test.dir/build.make apps/rosetta/CMakeFiles/rif_dock_test.dir/rif_dock_test.cc.o.provides.build
.PHONY : apps/rosetta/CMakeFiles/rif_dock_test.dir/rif_dock_test.cc.o.provides

apps/rosetta/CMakeFiles/rif_dock_test.dir/rif_dock_test.cc.o.provides.build: apps/rosetta/CMakeFiles/rif_dock_test.dir/rif_dock_test.cc.o


# Object files for target rif_dock_test
rif_dock_test_OBJECTS = \
"CMakeFiles/rif_dock_test.dir/rif_dock_test.cc.o"

# External object files for target rif_dock_test
rif_dock_test_EXTERNAL_OBJECTS =

apps/rosetta/rif_dock_test: apps/rosetta/CMakeFiles/rif_dock_test.dir/rif_dock_test.cc.o
apps/rosetta/rif_dock_test: apps/rosetta/CMakeFiles/rif_dock_test.dir/build.make
apps/rosetta/rif_dock_test: apps/rosetta/riflib/libriflib.so
apps/rosetta/rif_dock_test: apps/rosetta/CMakeFiles/rif_dock_test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/longxing/rosetta_gitclone/rifine/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable rif_dock_test"
	cd /home/longxing/rosetta_gitclone/rifine/build/apps/rosetta && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/rif_dock_test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
apps/rosetta/CMakeFiles/rif_dock_test.dir/build: apps/rosetta/rif_dock_test

.PHONY : apps/rosetta/CMakeFiles/rif_dock_test.dir/build

apps/rosetta/CMakeFiles/rif_dock_test.dir/requires: apps/rosetta/CMakeFiles/rif_dock_test.dir/rif_dock_test.cc.o.requires

.PHONY : apps/rosetta/CMakeFiles/rif_dock_test.dir/requires

apps/rosetta/CMakeFiles/rif_dock_test.dir/clean:
	cd /home/longxing/rosetta_gitclone/rifine/build/apps/rosetta && $(CMAKE_COMMAND) -P CMakeFiles/rif_dock_test.dir/cmake_clean.cmake
.PHONY : apps/rosetta/CMakeFiles/rif_dock_test.dir/clean

apps/rosetta/CMakeFiles/rif_dock_test.dir/depend:
	cd /home/longxing/rosetta_gitclone/rifine/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/longxing/rosetta_gitclone/rifine /home/longxing/rosetta_gitclone/rifine/apps/rosetta /home/longxing/rosetta_gitclone/rifine/build /home/longxing/rosetta_gitclone/rifine/build/apps/rosetta /home/longxing/rosetta_gitclone/rifine/build/apps/rosetta/CMakeFiles/rif_dock_test.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : apps/rosetta/CMakeFiles/rif_dock_test.dir/depend

