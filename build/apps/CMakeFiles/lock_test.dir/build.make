# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.7

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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.7.1/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.7.1/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/blayhem/GitHub/CompArch/atomic

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/blayhem/GitHub/CompArch/atomic/build

# Include any dependencies generated for this target.
include apps/CMakeFiles/lock_test.dir/depend.make

# Include the progress variables for this target.
include apps/CMakeFiles/lock_test.dir/progress.make

# Include the compile flags for this target's objects.
include apps/CMakeFiles/lock_test.dir/flags.make

apps/CMakeFiles/lock_test.dir/lock_app_base.cpp.o: apps/CMakeFiles/lock_test.dir/flags.make
apps/CMakeFiles/lock_test.dir/lock_app_base.cpp.o: ../apps/lock_app_base.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/blayhem/GitHub/CompArch/atomic/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object apps/CMakeFiles/lock_test.dir/lock_app_base.cpp.o"
	cd /Users/blayhem/GitHub/CompArch/atomic/build/apps && /usr/local/bin/g++-6   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/lock_test.dir/lock_app_base.cpp.o -c /Users/blayhem/GitHub/CompArch/atomic/apps/lock_app_base.cpp

apps/CMakeFiles/lock_test.dir/lock_app_base.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lock_test.dir/lock_app_base.cpp.i"
	cd /Users/blayhem/GitHub/CompArch/atomic/build/apps && /usr/local/bin/g++-6  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/blayhem/GitHub/CompArch/atomic/apps/lock_app_base.cpp > CMakeFiles/lock_test.dir/lock_app_base.cpp.i

apps/CMakeFiles/lock_test.dir/lock_app_base.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lock_test.dir/lock_app_base.cpp.s"
	cd /Users/blayhem/GitHub/CompArch/atomic/build/apps && /usr/local/bin/g++-6  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/blayhem/GitHub/CompArch/atomic/apps/lock_app_base.cpp -o CMakeFiles/lock_test.dir/lock_app_base.cpp.s

apps/CMakeFiles/lock_test.dir/lock_app_base.cpp.o.requires:

.PHONY : apps/CMakeFiles/lock_test.dir/lock_app_base.cpp.o.requires

apps/CMakeFiles/lock_test.dir/lock_app_base.cpp.o.provides: apps/CMakeFiles/lock_test.dir/lock_app_base.cpp.o.requires
	$(MAKE) -f apps/CMakeFiles/lock_test.dir/build.make apps/CMakeFiles/lock_test.dir/lock_app_base.cpp.o.provides.build
.PHONY : apps/CMakeFiles/lock_test.dir/lock_app_base.cpp.o.provides

apps/CMakeFiles/lock_test.dir/lock_app_base.cpp.o.provides.build: apps/CMakeFiles/lock_test.dir/lock_app_base.cpp.o


# Object files for target lock_test
lock_test_OBJECTS = \
"CMakeFiles/lock_test.dir/lock_app_base.cpp.o"

# External object files for target lock_test
lock_test_EXTERNAL_OBJECTS =

lock_test: apps/CMakeFiles/lock_test.dir/lock_app_base.cpp.o
lock_test: apps/CMakeFiles/lock_test.dir/build.make
lock_test: apps/CMakeFiles/lock_test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/blayhem/GitHub/CompArch/atomic/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../lock_test"
	cd /Users/blayhem/GitHub/CompArch/atomic/build/apps && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/lock_test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
apps/CMakeFiles/lock_test.dir/build: lock_test

.PHONY : apps/CMakeFiles/lock_test.dir/build

apps/CMakeFiles/lock_test.dir/requires: apps/CMakeFiles/lock_test.dir/lock_app_base.cpp.o.requires

.PHONY : apps/CMakeFiles/lock_test.dir/requires

apps/CMakeFiles/lock_test.dir/clean:
	cd /Users/blayhem/GitHub/CompArch/atomic/build/apps && $(CMAKE_COMMAND) -P CMakeFiles/lock_test.dir/cmake_clean.cmake
.PHONY : apps/CMakeFiles/lock_test.dir/clean

apps/CMakeFiles/lock_test.dir/depend:
	cd /Users/blayhem/GitHub/CompArch/atomic/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/blayhem/GitHub/CompArch/atomic /Users/blayhem/GitHub/CompArch/atomic/apps /Users/blayhem/GitHub/CompArch/atomic/build /Users/blayhem/GitHub/CompArch/atomic/build/apps /Users/blayhem/GitHub/CompArch/atomic/build/apps/CMakeFiles/lock_test.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : apps/CMakeFiles/lock_test.dir/depend

