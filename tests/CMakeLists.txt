set(UNIT_TEST_SOURCE_LIST
  tests.cpp)

set(TARGET_NAME tests)

# Add Genesis as dependency. You need to adapt the path to Genesis as needed.
include_directories(../genesis)
include_directories( ${GENESIS_INCLUDE_DIR} )
# Use all flags, linker options etc that Genesis exports.
add_definitions( ${GENESIS_DEFINITIONS} )
set( CMAKE_C_FLAGS          "${CMAKE_C_FLAGS}          ${GENESIS_C_FLAGS}")
set( CMAKE_CXX_FLAGS        "${CMAKE_CXX_FLAGS}        ${GENESIS_CXX_FLAGS}" )
set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${GENESIS_EXE_LINKER_FLAGS}" )
# Add QuartetScores
include_directories ("${PROJECT_SOURCE_DIR}/QuartetScores/src")

add_executable(${TARGET_NAME}
  main.cpp
  ${UNIT_TEST_SOURCE_LIST} )
target_link_libraries (${TARGET_NAME} ${GENESIS_LINK_LIBRARIES} )

target_include_directories(${TARGET_NAME}
  PUBLIC ../src/)

target_include_directories(${TARGET_NAME}
  PUBLIC ../externals/catch2/)

add_test(
  NAME ${TARGET_NAME}
  COMMAND ${TARGET_NAME} -s -o report.xml
  WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests")
