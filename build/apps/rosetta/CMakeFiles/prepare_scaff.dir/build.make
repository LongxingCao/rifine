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
CMAKE_SOURCE_DIR = /home/longxing/devel/rifine-devel

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/longxing/devel/rifine-devel/build

# Include any dependencies generated for this target.
include apps/rosetta/CMakeFiles/prepare_scaff.dir/depend.make

# Include the progress variables for this target.
include apps/rosetta/CMakeFiles/prepare_scaff.dir/progress.make

# Include the compile flags for this target's objects.
include apps/rosetta/CMakeFiles/prepare_scaff.dir/flags.make

apps/rosetta/CMakeFiles/prepare_scaff.dir/prepare_scaff.cc.o: apps/rosetta/CMakeFiles/prepare_scaff.dir/flags.make
apps/rosetta/CMakeFiles/prepare_scaff.dir/prepare_scaff.cc.o: ../apps/rosetta/prepare_scaff.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/longxing/devel/rifine-devel/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object apps/rosetta/CMakeFiles/prepare_scaff.dir/prepare_scaff.cc.o"
	cd /home/longxing/devel/rifine-devel/build/apps/rosetta && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/prepare_scaff.dir/prepare_scaff.cc.o -c /home/longxing/devel/rifine-devel/apps/rosetta/prepare_scaff.cc

apps/rosetta/CMakeFiles/prepare_scaff.dir/prepare_scaff.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/prepare_scaff.dir/prepare_scaff.cc.i"
	cd /home/longxing/devel/rifine-devel/build/apps/rosetta && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/longxing/devel/rifine-devel/apps/rosetta/prepare_scaff.cc > CMakeFiles/prepare_scaff.dir/prepare_scaff.cc.i

apps/rosetta/CMakeFiles/prepare_scaff.dir/prepare_scaff.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/prepare_scaff.dir/prepare_scaff.cc.s"
	cd /home/longxing/devel/rifine-devel/build/apps/rosetta && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/longxing/devel/rifine-devel/apps/rosetta/prepare_scaff.cc -o CMakeFiles/prepare_scaff.dir/prepare_scaff.cc.s

apps/rosetta/CMakeFiles/prepare_scaff.dir/prepare_scaff.cc.o.requires:

.PHONY : apps/rosetta/CMakeFiles/prepare_scaff.dir/prepare_scaff.cc.o.requires

apps/rosetta/CMakeFiles/prepare_scaff.dir/prepare_scaff.cc.o.provides: apps/rosetta/CMakeFiles/prepare_scaff.dir/prepare_scaff.cc.o.requires
	$(MAKE) -f apps/rosetta/CMakeFiles/prepare_scaff.dir/build.make apps/rosetta/CMakeFiles/prepare_scaff.dir/prepare_scaff.cc.o.provides.build
.PHONY : apps/rosetta/CMakeFiles/prepare_scaff.dir/prepare_scaff.cc.o.provides

apps/rosetta/CMakeFiles/prepare_scaff.dir/prepare_scaff.cc.o.provides.build: apps/rosetta/CMakeFiles/prepare_scaff.dir/prepare_scaff.cc.o


# Object files for target prepare_scaff
prepare_scaff_OBJECTS = \
"CMakeFiles/prepare_scaff.dir/prepare_scaff.cc.o"

# External object files for target prepare_scaff
prepare_scaff_EXTERNAL_OBJECTS =

apps/rosetta/prepare_scaff: apps/rosetta/CMakeFiles/prepare_scaff.dir/prepare_scaff.cc.o
apps/rosetta/prepare_scaff: apps/rosetta/CMakeFiles/prepare_scaff.dir/build.make
apps/rosetta/prepare_scaff: apps/rosetta/riflib/libriflib.so
apps/rosetta/prepare_scaff: apps/rosetta/CMakeFiles/prepare_scaff.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/longxing/devel/rifine-devel/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable prepare_scaff"
	cd /home/longxing/devel/rifine-devel/build/apps/rosetta && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/prepare_scaff.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
apps/rosetta/CMakeFiles/prepare_scaff.dir/build: apps/rosetta/prepare_scaff

.PHONY : apps/rosetta/CMakeFiles/prepare_scaff.dir/build

apps/rosetta/CMakeFiles/prepare_scaff.dir/requires: apps/rosetta/CMakeFiles/prepare_scaff.dir/prepare_scaff.cc.o.requires

.PHONY : apps/rosetta/CMakeFiles/prepare_scaff.dir/requires

apps/rosetta/CMakeFiles/prepare_scaff.dir/clean:
	cd /home/longxing/devel/rifine-devel/build/apps/rosetta && $(CMAKE_COMMAND) -P CMakeFiles/prepare_scaff.dir/cmake_clean.cmake
.PHONY : apps/rosetta/CMakeFiles/prepare_scaff.dir/clean

apps/rosetta/CMakeFiles/prepare_scaff.dir/depend:
	cd /home/longxing/devel/rifine-devel/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/longxing/devel/rifine-devel /home/longxing/devel/rifine-devel/apps/rosetta /home/longxing/devel/rifine-devel/build /home/longxing/devel/rifine-devel/build/apps/rosetta /home/longxing/devel/rifine-devel/build/apps/rosetta/CMakeFiles/prepare_scaff.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : apps/rosetta/CMakeFiles/prepare_scaff.dir/depend

