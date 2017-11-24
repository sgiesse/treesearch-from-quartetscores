cmake_minimum_required (VERSION 2.6)
project (treesearch)


set (CMAKE_BUILD_TYPE RELEASE)

# --------------------------------------------------------------------------------------------------
#   Compiler and Linker Options
# --------------------------------------------------------------------------------------------------

# Set the warn flags to a very high level - except unknown pragmas, as this causes needless
# warnings with OpenMP and other pragma-based techniques.
set (WARN_FLAGS              "-Wall -Wextra -Wno-unknown-pragmas -pedantic -pedantic-errors")
set (CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -std=c++11 ${WARN_FLAGS}")

set (CMAKE_CXX_FLAGS_DEBUG   "-O2 -g -ggdb3 -DDEBUG -D_GLIBCXX_DEBUG")
set (CMAKE_CXX_FLAGS_RELEASE "-O2       -DNDEBUG")



# Add Genesis as dependency. You need to adapt the path to Genesis as needed.
add_subdirectory(${PROJECT_SOURCE_DIR}/genesis)
include_directories( ${GENESIS_INCLUDE_DIR} )

# Use all flags, linker options etc that Genesis exports.
add_definitions( ${GENESIS_DEFINITIONS} )
set( CMAKE_C_FLAGS          "${CMAKE_C_FLAGS}          ${GENESIS_C_FLAGS}")
set( CMAKE_CXX_FLAGS        "${CMAKE_CXX_FLAGS}        ${GENESIS_CXX_FLAGS}" )
set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${GENESIS_EXE_LINKER_FLAGS}" )


# Add QuartetScores
include_directories ("${PROJECT_SOURCE_DIR}/QuartetScores/src")

include_directories( "tclap/include" )


include_directories( "src" )

add_executable(treesearch "src/main.cpp")
target_link_libraries (treesearch ${GENESIS_LINK_LIBRARIES} )
