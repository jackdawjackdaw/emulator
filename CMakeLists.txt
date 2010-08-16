cmake_minimum_required (VERSION 2.6)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")
set(CMAKE_INSTALL_PREFIX "/home/einstein/grad/cec24/local")
project (Emulator)

SET(${project_name}_major_version 1)
SET(${project_name}_minor_version 0)


# the RPATH to be used when installing
#SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib/libEmu")

find_package(GSL REQUIRED)
message("GSL_INCLUDE_DIRS: ${GSL_INCLUDE_DIR}")
message("CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
message("CMAKE_INSTALL_PATH: ${CMAKE_INSTALL_PREFIX}")

INCLUDE_DIRECTORIES(${GSL_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/libEmu)
MESSAGE("include dirs: ${CMAKE_include_directories}")
SET(libs ${libs} ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})


ADD_DEFINITIONS("-std=c99 -pg -O0 -g -Wextra -pedantic")

add_subdirectory(src)

## this is dangerous, its very creaky right now.
INSTALL_PROGRAMS("/bin" FILES src/estimator src/emulator)
INSTALL_FILES("/lib" FILES src/libEmu/libEmu.so src/libEmu/libLBFGS.so src/libRbind/libRBIND.so)
INSTALL_FILES("/include/libEmu/" FILES src/libEmu/emulate-fns.h src/libEmu/emulator.h src/libEmu/estimate_threaded.h src/libEmu/estimator.h src/libEmu/lbfgs.h src/libEmu/maxlbfgs.h src/libEmu/regression.h)
INSTALL_FILES("/include/libRbind/" FILES src/libRbind/rbind.h)


