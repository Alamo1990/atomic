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
include apps/CMakeFiles/atomic_adv.dir/depend.make

# Include the progress variables for this target.
include apps/CMakeFiles/atomic_adv.dir/progress.make

# Include the compile flags for this target's objects.
include apps/CMakeFiles/atomic_adv.dir/flags.make

apps/CMakeFiles/atomic_adv.dir/atomic_app_adv.cpp.o: apps/CMakeFiles/atomic_adv.dir/flags.make
apps/CMakeFiles/atomic_adv.dir/atomic_app_adv.cpp.o: ../apps/atomic_app_adv.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/blayhem/GitHub/CompArch/atomic/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object apps/CMakeFiles/atomic_adv.dir/atomic_app_adv.cpp.o"
	cd /Users/blayhem/GitHub/CompArch/atomic/build/apps && /usr/local/bin/g++-6   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/atomic_adv.dir/atomic_app_adv.cpp.o -c /Users/blayhem/GitHub/CompArch/atomic/apps/atomic_app_adv.cpp

apps/CMakeFiles/atomic_adv.dir/atomic_app_adv.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/atomic_adv.dir/atomic_app_adv.cpp.i"
	cd /Users/blayhem/GitHub/CompArch/atomic/build/apps && /usr/local/bin/g++-6  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/blayhem/GitHub/CompArch/atomic/apps/atomic_app_adv.cpp > CMakeFiles/atomic_adv.dir/atomic_app_adv.cpp.i

apps/CMakeFiles/atomic_adv.dir/atomic_app_adv.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/atomic_adv.dir/atomic_app_adv.cpp.s"
	cd /Users/blayhem/GitHub/CompArch/atomic/build/apps && /usr/local/bin/g++-6  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/blayhem/GitHub/CompArch/atomic/apps/atomic_app_adv.cpp -o CMakeFiles/atomic_adv.dir/atomic_app_adv.cpp.s

apps/CMakeFiles/atomic_adv.dir/atomic_app_adv.cpp.o.requires:

.PHONY : apps/CMakeFiles/atomic_adv.dir/atomic_app_adv.cpp.o.requires

apps/CMakeFiles/atomic_adv.dir/atomic_app_adv.cpp.o.provides: apps/CMakeFiles/atomic_adv.dir/atomic_app_adv.cpp.o.requires
	$(MAKE) -f apps/CMakeFiles/atomic_adv.dir/build.make apps/CMakeFiles/atomic_adv.dir/atomic_app_adv.cpp.o.provides.build
.PHONY : apps/CMakeFiles/atomic_adv.dir/atomic_app_adv.cpp.o.provides

apps/CMakeFiles/atomic_adv.dir/atomic_app_adv.cpp.o.provides.build: apps/CMakeFiles/atomic_adv.dir/atomic_app_adv.cpp.o


# Object files for target atomic_adv
atomic_adv_OBJECTS = \
"CMakeFiles/atomic_adv.dir/atomic_app_adv.cpp.o"

# External object files for target atomic_adv
atomic_adv_EXTERNAL_OBJECTS =

atomic_adv: apps/CMakeFiles/atomic_adv.dir/atomic_app_adv.cpp.o
atomic_adv: apps/CMakeFiles/atomic_adv.dir/build.make
atomic_adv: apps/CMakeFiles/atomic_adv.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/blayhem/GitHub/CompArch/atomic/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../atomic_adv"
	cd /Users/blayhem/GitHub/CompArch/atomic/build/apps && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/atomic_adv.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
apps/CMakeFiles/atomic_adv.dir/build: atomic_adv

.PHONY : apps/CMakeFiles/atomic_adv.dir/build

apps/CMakeFiles/atomic_adv.dir/requires: apps/CMakeFiles/atomic_adv.dir/atomic_app_adv.cpp.o.requires

.PHONY : apps/CMakeFiles/atomic_adv.dir/requires

apps/CMakeFiles/atomic_adv.dir/clean:
	cd /Users/blayhem/GitHub/CompArch/atomic/build/apps && $(CMAKE_COMMAND) -P CMakeFiles/atomic_adv.dir/cmake_clean.cmake
.PHONY : apps/CMakeFiles/atomic_adv.dir/clean

apps/CMakeFiles/atomic_adv.dir/depend:
	cd /Users/blayhem/GitHub/CompArch/atomic/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/blayhem/GitHub/CompArch/atomic /Users/blayhem/GitHub/CompArch/atomic/apps /Users/blayhem/GitHub/CompArch/atomic/build /Users/blayhem/GitHub/CompArch/atomic/build/apps /Users/blayhem/GitHub/CompArch/atomic/build/apps/CMakeFiles/atomic_adv.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : apps/CMakeFiles/atomic_adv.dir/depend

