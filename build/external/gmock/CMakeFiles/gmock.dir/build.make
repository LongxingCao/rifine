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
include external/gmock/CMakeFiles/gmock.dir/depend.make

# Include the progress variables for this target.
include external/gmock/CMakeFiles/gmock.dir/progress.make

# Include the compile flags for this target's objects.
include external/gmock/CMakeFiles/gmock.dir/flags.make

external/gmock/CMakeFiles/gmock.dir/gtest/src/gtest-all.cc.o: external/gmock/CMakeFiles/gmock.dir/flags.make
external/gmock/CMakeFiles/gmock.dir/gtest/src/gtest-all.cc.o: ../external/gmock/gtest/src/gtest-all.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/longxing/devel/rifine-devel/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object external/gmock/CMakeFiles/gmock.dir/gtest/src/gtest-all.cc.o"
	cd /home/longxing/devel/rifine-devel/build/external/gmock && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/gmock.dir/gtest/src/gtest-all.cc.o -c /home/longxing/devel/rifine-devel/external/gmock/gtest/src/gtest-all.cc

external/gmock/CMakeFiles/gmock.dir/gtest/src/gtest-all.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/gmock.dir/gtest/src/gtest-all.cc.i"
	cd /home/longxing/devel/rifine-devel/build/external/gmock && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/longxing/devel/rifine-devel/external/gmock/gtest/src/gtest-all.cc > CMakeFiles/gmock.dir/gtest/src/gtest-all.cc.i

external/gmock/CMakeFiles/gmock.dir/gtest/src/gtest-all.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/gmock.dir/gtest/src/gtest-all.cc.s"
	cd /home/longxing/devel/rifine-devel/build/external/gmock && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/longxing/devel/rifine-devel/external/gmock/gtest/src/gtest-all.cc -o CMakeFiles/gmock.dir/gtest/src/gtest-all.cc.s

external/gmock/CMakeFiles/gmock.dir/gtest/src/gtest-all.cc.o.requires:

.PHONY : external/gmock/CMakeFiles/gmock.dir/gtest/src/gtest-all.cc.o.requires

external/gmock/CMakeFiles/gmock.dir/gtest/src/gtest-all.cc.o.provides: external/gmock/CMakeFiles/gmock.dir/gtest/src/gtest-all.cc.o.requires
	$(MAKE) -f external/gmock/CMakeFiles/gmock.dir/build.make external/gmock/CMakeFiles/gmock.dir/gtest/src/gtest-all.cc.o.provides.build
.PHONY : external/gmock/CMakeFiles/gmock.dir/gtest/src/gtest-all.cc.o.provides

external/gmock/CMakeFiles/gmock.dir/gtest/src/gtest-all.cc.o.provides.build: external/gmock/CMakeFiles/gmock.dir/gtest/src/gtest-all.cc.o


external/gmock/CMakeFiles/gmock.dir/src/gmock-all.cc.o: external/gmock/CMakeFiles/gmock.dir/flags.make
external/gmock/CMakeFiles/gmock.dir/src/gmock-all.cc.o: ../external/gmock/src/gmock-all.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/longxing/devel/rifine-devel/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object external/gmock/CMakeFiles/gmock.dir/src/gmock-all.cc.o"
	cd /home/longxing/devel/rifine-devel/build/external/gmock && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/gmock.dir/src/gmock-all.cc.o -c /home/longxing/devel/rifine-devel/external/gmock/src/gmock-all.cc

external/gmock/CMakeFiles/gmock.dir/src/gmock-all.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/gmock.dir/src/gmock-all.cc.i"
	cd /home/longxing/devel/rifine-devel/build/external/gmock && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/longxing/devel/rifine-devel/external/gmock/src/gmock-all.cc > CMakeFiles/gmock.dir/src/gmock-all.cc.i

external/gmock/CMakeFiles/gmock.dir/src/gmock-all.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/gmock.dir/src/gmock-all.cc.s"
	cd /home/longxing/devel/rifine-devel/build/external/gmock && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/longxing/devel/rifine-devel/external/gmock/src/gmock-all.cc -o CMakeFiles/gmock.dir/src/gmock-all.cc.s

external/gmock/CMakeFiles/gmock.dir/src/gmock-all.cc.o.requires:

.PHONY : external/gmock/CMakeFiles/gmock.dir/src/gmock-all.cc.o.requires

external/gmock/CMakeFiles/gmock.dir/src/gmock-all.cc.o.provides: external/gmock/CMakeFiles/gmock.dir/src/gmock-all.cc.o.requires
	$(MAKE) -f external/gmock/CMakeFiles/gmock.dir/build.make external/gmock/CMakeFiles/gmock.dir/src/gmock-all.cc.o.provides.build
.PHONY : external/gmock/CMakeFiles/gmock.dir/src/gmock-all.cc.o.provides

external/gmock/CMakeFiles/gmock.dir/src/gmock-all.cc.o.provides.build: external/gmock/CMakeFiles/gmock.dir/src/gmock-all.cc.o


# Object files for target gmock
gmock_OBJECTS = \
"CMakeFiles/gmock.dir/gtest/src/gtest-all.cc.o" \
"CMakeFiles/gmock.dir/src/gmock-all.cc.o"

# External object files for target gmock
gmock_EXTERNAL_OBJECTS =

external/gmock/libgmock.a: external/gmock/CMakeFiles/gmock.dir/gtest/src/gtest-all.cc.o
external/gmock/libgmock.a: external/gmock/CMakeFiles/gmock.dir/src/gmock-all.cc.o
external/gmock/libgmock.a: external/gmock/CMakeFiles/gmock.dir/build.make
external/gmock/libgmock.a: external/gmock/CMakeFiles/gmock.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/longxing/devel/rifine-devel/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX static library libgmock.a"
	cd /home/longxing/devel/rifine-devel/build/external/gmock && $(CMAKE_COMMAND) -P CMakeFiles/gmock.dir/cmake_clean_target.cmake
	cd /home/longxing/devel/rifine-devel/build/external/gmock && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/gmock.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
external/gmock/CMakeFiles/gmock.dir/build: external/gmock/libgmock.a

.PHONY : external/gmock/CMakeFiles/gmock.dir/build

external/gmock/CMakeFiles/gmock.dir/requires: external/gmock/CMakeFiles/gmock.dir/gtest/src/gtest-all.cc.o.requires
external/gmock/CMakeFiles/gmock.dir/requires: external/gmock/CMakeFiles/gmock.dir/src/gmock-all.cc.o.requires

.PHONY : external/gmock/CMakeFiles/gmock.dir/requires

external/gmock/CMakeFiles/gmock.dir/clean:
	cd /home/longxing/devel/rifine-devel/build/external/gmock && $(CMAKE_COMMAND) -P CMakeFiles/gmock.dir/cmake_clean.cmake
.PHONY : external/gmock/CMakeFiles/gmock.dir/clean

external/gmock/CMakeFiles/gmock.dir/depend:
	cd /home/longxing/devel/rifine-devel/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/longxing/devel/rifine-devel /home/longxing/devel/rifine-devel/external/gmock /home/longxing/devel/rifine-devel/build /home/longxing/devel/rifine-devel/build/external/gmock /home/longxing/devel/rifine-devel/build/external/gmock/CMakeFiles/gmock.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : external/gmock/CMakeFiles/gmock.dir/depend

