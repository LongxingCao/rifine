# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/longxing/devel/rifine

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/longxing/devel/rifine/build

# Include any dependencies generated for this target.
include apps/rosetta/CMakeFiles/rifgen.dir/depend.make

# Include the progress variables for this target.
include apps/rosetta/CMakeFiles/rifgen.dir/progress.make

# Include the compile flags for this target's objects.
include apps/rosetta/CMakeFiles/rifgen.dir/flags.make

apps/rosetta/CMakeFiles/rifgen.dir/rifgen.cc.o: apps/rosetta/CMakeFiles/rifgen.dir/flags.make
apps/rosetta/CMakeFiles/rifgen.dir/rifgen.cc.o: ../apps/rosetta/rifgen.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/longxing/devel/rifine/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object apps/rosetta/CMakeFiles/rifgen.dir/rifgen.cc.o"
	cd /home/longxing/devel/rifine/build/apps/rosetta && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/rifgen.dir/rifgen.cc.o -c /home/longxing/devel/rifine/apps/rosetta/rifgen.cc

apps/rosetta/CMakeFiles/rifgen.dir/rifgen.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/rifgen.dir/rifgen.cc.i"
	cd /home/longxing/devel/rifine/build/apps/rosetta && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/longxing/devel/rifine/apps/rosetta/rifgen.cc > CMakeFiles/rifgen.dir/rifgen.cc.i

apps/rosetta/CMakeFiles/rifgen.dir/rifgen.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/rifgen.dir/rifgen.cc.s"
	cd /home/longxing/devel/rifine/build/apps/rosetta && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/longxing/devel/rifine/apps/rosetta/rifgen.cc -o CMakeFiles/rifgen.dir/rifgen.cc.s

apps/rosetta/CMakeFiles/rifgen.dir/rifgen.cc.o.requires:

.PHONY : apps/rosetta/CMakeFiles/rifgen.dir/rifgen.cc.o.requires

apps/rosetta/CMakeFiles/rifgen.dir/rifgen.cc.o.provides: apps/rosetta/CMakeFiles/rifgen.dir/rifgen.cc.o.requires
	$(MAKE) -f apps/rosetta/CMakeFiles/rifgen.dir/build.make apps/rosetta/CMakeFiles/rifgen.dir/rifgen.cc.o.provides.build
.PHONY : apps/rosetta/CMakeFiles/rifgen.dir/rifgen.cc.o.provides

apps/rosetta/CMakeFiles/rifgen.dir/rifgen.cc.o.provides.build: apps/rosetta/CMakeFiles/rifgen.dir/rifgen.cc.o


# Object files for target rifgen
rifgen_OBJECTS = \
"CMakeFiles/rifgen.dir/rifgen.cc.o"

# External object files for target rifgen
rifgen_EXTERNAL_OBJECTS =

apps/rosetta/rifgen: apps/rosetta/CMakeFiles/rifgen.dir/rifgen.cc.o
apps/rosetta/rifgen: apps/rosetta/CMakeFiles/rifgen.dir/build.make
apps/rosetta/rifgen: apps/rosetta/riflib/libriflib.so
apps/rosetta/rifgen: apps/rosetta/CMakeFiles/rifgen.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/longxing/devel/rifine/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable rifgen"
	cd /home/longxing/devel/rifine/build/apps/rosetta && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/rifgen.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
apps/rosetta/CMakeFiles/rifgen.dir/build: apps/rosetta/rifgen

.PHONY : apps/rosetta/CMakeFiles/rifgen.dir/build

apps/rosetta/CMakeFiles/rifgen.dir/requires: apps/rosetta/CMakeFiles/rifgen.dir/rifgen.cc.o.requires

.PHONY : apps/rosetta/CMakeFiles/rifgen.dir/requires

apps/rosetta/CMakeFiles/rifgen.dir/clean:
	cd /home/longxing/devel/rifine/build/apps/rosetta && $(CMAKE_COMMAND) -P CMakeFiles/rifgen.dir/cmake_clean.cmake
.PHONY : apps/rosetta/CMakeFiles/rifgen.dir/clean

apps/rosetta/CMakeFiles/rifgen.dir/depend:
	cd /home/longxing/devel/rifine/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/longxing/devel/rifine /home/longxing/devel/rifine/apps/rosetta /home/longxing/devel/rifine/build /home/longxing/devel/rifine/build/apps/rosetta /home/longxing/devel/rifine/build/apps/rosetta/CMakeFiles/rifgen.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : apps/rosetta/CMakeFiles/rifgen.dir/depend

