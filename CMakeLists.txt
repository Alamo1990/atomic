cmake_minimum_required (VERSION 2.8)

set(PROJECT_NAME_STR bbuffer)
project (${PROJECT_NAME_STR})

file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/output")
# Local project include directory
include_directories("${CMAKE_SOURCE_DIR}/include")

# Build binary directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")

#list(APPEND CMAKE_CXX_FLAGS "-std=c++11 -stdlib=libc++ -Wall -Wextra -Wno-deprecated -Werror -pedantic -pedantic-errors -pthread -fsanitize=thread")
list(APPEND CMAKE_CXX_FLAGS "-std=c++11 -Wall -Wextra -Wno-deprecated -Werror -pedantic -pedantic-errors -pthread")

#Process subdirectories
add_subdirectory(apps)
