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
include schemelib/test/CMakeFiles/quick_test_libscheme.dir/depend.make

# Include the progress variables for this target.
include schemelib/test/CMakeFiles/quick_test_libscheme.dir/progress.make

# Include the compile flags for this target's objects.
include schemelib/test/CMakeFiles/quick_test_libscheme.dir/flags.make

schemelib/test/CMakeFiles/quick_test_libscheme.dir/quick_test.cc.o: schemelib/test/CMakeFiles/quick_test_libscheme.dir/flags.make
schemelib/test/CMakeFiles/quick_test_libscheme.dir/quick_test.cc.o: ../schemelib/test/quick_test.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/longxing/devel/rifine/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object schemelib/test/CMakeFiles/quick_test_libscheme.dir/quick_test.cc.o"
	cd /home/longxing/devel/rifine/build/schemelib/test && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/quick_test_libscheme.dir/quick_test.cc.o -c /home/longxing/devel/rifine/schemelib/test/quick_test.cc

schemelib/test/CMakeFiles/quick_test_libscheme.dir/quick_test.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/quick_test_libscheme.dir/quick_test.cc.i"
	cd /home/longxing/devel/rifine/build/schemelib/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/longxing/devel/rifine/schemelib/test/quick_test.cc > CMakeFiles/quick_test_libscheme.dir/quick_test.cc.i

schemelib/test/CMakeFiles/quick_test_libscheme.dir/quick_test.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/quick_test_libscheme.dir/quick_test.cc.s"
	cd /home/longxing/devel/rifine/build/schemelib/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/longxing/devel/rifine/schemelib/test/quick_test.cc -o CMakeFiles/quick_test_libscheme.dir/quick_test.cc.s

schemelib/test/CMakeFiles/quick_test_libscheme.dir/quick_test.cc.o.requires:

.PHONY : schemelib/test/CMakeFiles/quick_test_libscheme.dir/quick_test.cc.o.requires

schemelib/test/CMakeFiles/quick_test_libscheme.dir/quick_test.cc.o.provides: schemelib/test/CMakeFiles/quick_test_libscheme.dir/quick_test.cc.o.requires
	$(MAKE) -f schemelib/test/CMakeFiles/quick_test_libscheme.dir/build.make schemelib/test/CMakeFiles/quick_test_libscheme.dir/quick_test.cc.o.provides.build
.PHONY : schemelib/test/CMakeFiles/quick_test_libscheme.dir/quick_test.cc.o.provides

schemelib/test/CMakeFiles/quick_test_libscheme.dir/quick_test.cc.o.provides.build: schemelib/test/CMakeFiles/quick_test_libscheme.dir/quick_test.cc.o


# Object files for target quick_test_libscheme
quick_test_libscheme_OBJECTS = \
"CMakeFiles/quick_test_libscheme.dir/quick_test.cc.o"

# External object files for target quick_test_libscheme
quick_test_libscheme_EXTERNAL_OBJECTS =

schemelib/test/quick_test_libscheme: schemelib/test/CMakeFiles/quick_test_libscheme.dir/quick_test.cc.o
schemelib/test/quick_test_libscheme: schemelib/test/CMakeFiles/quick_test_libscheme.dir/build.make
schemelib/test/quick_test_libscheme: schemelib/test/libscheme.a
schemelib/test/quick_test_libscheme: external/gmock/libgmock.a
schemelib/test/quick_test_libscheme: schemelib/test/CMakeFiles/quick_test_libscheme.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/longxing/devel/rifine/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable quick_test_libscheme"
	cd /home/longxing/devel/rifine/build/schemelib/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/quick_test_libscheme.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
schemelib/test/CMakeFiles/quick_test_libscheme.dir/build: schemelib/test/quick_test_libscheme

.PHONY : schemelib/test/CMakeFiles/quick_test_libscheme.dir/build

schemelib/test/CMakeFiles/quick_test_libscheme.dir/requires: schemelib/test/CMakeFiles/quick_test_libscheme.dir/quick_test.cc.o.requires

.PHONY : schemelib/test/CMakeFiles/quick_test_libscheme.dir/requires

schemelib/test/CMakeFiles/quick_test_libscheme.dir/clean:
	cd /home/longxing/devel/rifine/build/schemelib/test && $(CMAKE_COMMAND) -P CMakeFiles/quick_test_libscheme.dir/cmake_clean.cmake
.PHONY : schemelib/test/CMakeFiles/quick_test_libscheme.dir/clean

schemelib/test/CMakeFiles/quick_test_libscheme.dir/depend:
	cd /home/longxing/devel/rifine/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/longxing/devel/rifine /home/longxing/devel/rifine/schemelib/test /home/longxing/devel/rifine/build /home/longxing/devel/rifine/build/schemelib/test /home/longxing/devel/rifine/build/schemelib/test/CMakeFiles/quick_test_libscheme.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : schemelib/test/CMakeFiles/quick_test_libscheme.dir/depend

