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
CMAKE_SOURCE_DIR = /home/blayhem/Github/atomic

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/blayhem/Github/atomic/build

# Include any dependencies generated for this target.
include apps/CMakeFiles/atomic_test.dir/depend.make

# Include the progress variables for this target.
include apps/CMakeFiles/atomic_test.dir/progress.make

# Include the compile flags for this target's objects.
include apps/CMakeFiles/atomic_test.dir/flags.make

apps/CMakeFiles/atomic_test.dir/atomic_app_base.cpp.o: apps/CMakeFiles/atomic_test.dir/flags.make
apps/CMakeFiles/atomic_test.dir/atomic_app_base.cpp.o: ../apps/atomic_app_base.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/blayhem/Github/atomic/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object apps/CMakeFiles/atomic_test.dir/atomic_app_base.cpp.o"
	cd /home/blayhem/Github/atomic/build/apps && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/atomic_test.dir/atomic_app_base.cpp.o -c /home/blayhem/Github/atomic/apps/atomic_app_base.cpp

apps/CMakeFiles/atomic_test.dir/atomic_app_base.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/atomic_test.dir/atomic_app_base.cpp.i"
	cd /home/blayhem/Github/atomic/build/apps && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/blayhem/Github/atomic/apps/atomic_app_base.cpp > CMakeFiles/atomic_test.dir/atomic_app_base.cpp.i

apps/CMakeFiles/atomic_test.dir/atomic_app_base.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/atomic_test.dir/atomic_app_base.cpp.s"
	cd /home/blayhem/Github/atomic/build/apps && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/blayhem/Github/atomic/apps/atomic_app_base.cpp -o CMakeFiles/atomic_test.dir/atomic_app_base.cpp.s

apps/CMakeFiles/atomic_test.dir/atomic_app_base.cpp.o.requires:

.PHONY : apps/CMakeFiles/atomic_test.dir/atomic_app_base.cpp.o.requires

apps/CMakeFiles/atomic_test.dir/atomic_app_base.cpp.o.provides: apps/CMakeFiles/atomic_test.dir/atomic_app_base.cpp.o.requires
	$(MAKE) -f apps/CMakeFiles/atomic_test.dir/build.make apps/CMakeFiles/atomic_test.dir/atomic_app_base.cpp.o.provides.build
.PHONY : apps/CMakeFiles/atomic_test.dir/atomic_app_base.cpp.o.provides

apps/CMakeFiles/atomic_test.dir/atomic_app_base.cpp.o.provides.build: apps/CMakeFiles/atomic_test.dir/atomic_app_base.cpp.o


# Object files for target atomic_test
atomic_test_OBJECTS = \
"CMakeFiles/atomic_test.dir/atomic_app_base.cpp.o"

# External object files for target atomic_test
atomic_test_EXTERNAL_OBJECTS =

atomic_test: apps/CMakeFiles/atomic_test.dir/atomic_app_base.cpp.o
atomic_test: apps/CMakeFiles/atomic_test.dir/build.make
atomic_test: apps/CMakeFiles/atomic_test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/blayhem/Github/atomic/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../atomic_test"
	cd /home/blayhem/Github/atomic/build/apps && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/atomic_test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
apps/CMakeFiles/atomic_test.dir/build: atomic_test

.PHONY : apps/CMakeFiles/atomic_test.dir/build

apps/CMakeFiles/atomic_test.dir/requires: apps/CMakeFiles/atomic_test.dir/atomic_app_base.cpp.o.requires

.PHONY : apps/CMakeFiles/atomic_test.dir/requires

apps/CMakeFiles/atomic_test.dir/clean:
	cd /home/blayhem/Github/atomic/build/apps && $(CMAKE_COMMAND) -P CMakeFiles/atomic_test.dir/cmake_clean.cmake
.PHONY : apps/CMakeFiles/atomic_test.dir/clean

apps/CMakeFiles/atomic_test.dir/depend:
	cd /home/blayhem/Github/atomic/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/blayhem/Github/atomic /home/blayhem/Github/atomic/apps /home/blayhem/Github/atomic/build /home/blayhem/Github/atomic/build/apps /home/blayhem/Github/atomic/build/apps/CMakeFiles/atomic_test.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : apps/CMakeFiles/atomic_test.dir/depend

