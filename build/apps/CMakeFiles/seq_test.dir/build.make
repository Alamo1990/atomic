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
CMAKE_SOURCE_DIR = /home/alamo/Documents/arquitecturaComputadors/atomic

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/alamo/Documents/arquitecturaComputadors/atomic/build

# Include any dependencies generated for this target.
include apps/CMakeFiles/seq_test.dir/depend.make

# Include the progress variables for this target.
include apps/CMakeFiles/seq_test.dir/progress.make

# Include the compile flags for this target's objects.
include apps/CMakeFiles/seq_test.dir/flags.make

apps/CMakeFiles/seq_test.dir/seq_app.cpp.o: apps/CMakeFiles/seq_test.dir/flags.make
apps/CMakeFiles/seq_test.dir/seq_app.cpp.o: ../apps/seq_app.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/alamo/Documents/arquitecturaComputadors/atomic/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object apps/CMakeFiles/seq_test.dir/seq_app.cpp.o"
	cd /home/alamo/Documents/arquitecturaComputadors/atomic/build/apps && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/seq_test.dir/seq_app.cpp.o -c /home/alamo/Documents/arquitecturaComputadors/atomic/apps/seq_app.cpp

apps/CMakeFiles/seq_test.dir/seq_app.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/seq_test.dir/seq_app.cpp.i"
	cd /home/alamo/Documents/arquitecturaComputadors/atomic/build/apps && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/alamo/Documents/arquitecturaComputadors/atomic/apps/seq_app.cpp > CMakeFiles/seq_test.dir/seq_app.cpp.i

apps/CMakeFiles/seq_test.dir/seq_app.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/seq_test.dir/seq_app.cpp.s"
	cd /home/alamo/Documents/arquitecturaComputadors/atomic/build/apps && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/alamo/Documents/arquitecturaComputadors/atomic/apps/seq_app.cpp -o CMakeFiles/seq_test.dir/seq_app.cpp.s

apps/CMakeFiles/seq_test.dir/seq_app.cpp.o.requires:

.PHONY : apps/CMakeFiles/seq_test.dir/seq_app.cpp.o.requires

apps/CMakeFiles/seq_test.dir/seq_app.cpp.o.provides: apps/CMakeFiles/seq_test.dir/seq_app.cpp.o.requires
	$(MAKE) -f apps/CMakeFiles/seq_test.dir/build.make apps/CMakeFiles/seq_test.dir/seq_app.cpp.o.provides.build
.PHONY : apps/CMakeFiles/seq_test.dir/seq_app.cpp.o.provides

apps/CMakeFiles/seq_test.dir/seq_app.cpp.o.provides.build: apps/CMakeFiles/seq_test.dir/seq_app.cpp.o


# Object files for target seq_test
seq_test_OBJECTS = \
"CMakeFiles/seq_test.dir/seq_app.cpp.o"

# External object files for target seq_test
seq_test_EXTERNAL_OBJECTS =

seq_test: apps/CMakeFiles/seq_test.dir/seq_app.cpp.o
seq_test: apps/CMakeFiles/seq_test.dir/build.make
seq_test: apps/CMakeFiles/seq_test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/alamo/Documents/arquitecturaComputadors/atomic/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../seq_test"
	cd /home/alamo/Documents/arquitecturaComputadors/atomic/build/apps && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/seq_test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
apps/CMakeFiles/seq_test.dir/build: seq_test

.PHONY : apps/CMakeFiles/seq_test.dir/build

apps/CMakeFiles/seq_test.dir/requires: apps/CMakeFiles/seq_test.dir/seq_app.cpp.o.requires

.PHONY : apps/CMakeFiles/seq_test.dir/requires

apps/CMakeFiles/seq_test.dir/clean:
	cd /home/alamo/Documents/arquitecturaComputadors/atomic/build/apps && $(CMAKE_COMMAND) -P CMakeFiles/seq_test.dir/cmake_clean.cmake
.PHONY : apps/CMakeFiles/seq_test.dir/clean

apps/CMakeFiles/seq_test.dir/depend:
	cd /home/alamo/Documents/arquitecturaComputadors/atomic/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/alamo/Documents/arquitecturaComputadors/atomic /home/alamo/Documents/arquitecturaComputadors/atomic/apps /home/alamo/Documents/arquitecturaComputadors/atomic/build /home/alamo/Documents/arquitecturaComputadors/atomic/build/apps /home/alamo/Documents/arquitecturaComputadors/atomic/build/apps/CMakeFiles/seq_test.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : apps/CMakeFiles/seq_test.dir/depend

