cmake_minimum_required (VERSION 2.6)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")
set(CMAKE_INSTALL_PREFIX "/home/einstein/grad/cec24/local")
project (Emulator)

SET(${project_name}_major_version 1)
SET(${project_name}_minor_version 0)


# the RPATH to be used when installing
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

find_package(GSL REQUIRED)
message("GSL_INCLUDE_DIRS: ${GSL_INCLUDE_DIR}")
message("CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")

INCLUDE_DIRECTORIES(${GSL_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/libEmu)
MESSAGE("include dirs: ${CMAKE_include_directories}")
SET(libs ${libs} ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})


ADD_DEFINITIONS("-std=c99 -pg -O0 -g -Wextra -pedantic")

add_subdirectory(src)

INSTALL_PROGRAMS("/bin" FILES src/estimator src/emulator)
INSTALL_FILES("/lib" FILES src/libEmu/libEmu.so src/libEmu/libLBFGS.so src/libRBIND.so)

