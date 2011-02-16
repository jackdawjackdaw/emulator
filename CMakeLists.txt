cmake_minimum_required (VERSION 2.6)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")
set(CMAKE_INSTALL_PREFIX "/home/einstein/grad/cec24/local")
project (Emulator)

SET(${project_name}_major_version 1)
SET(${project_name}_minor_version 0)

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 

# the RPATH to be used when installing
SET(CMAKE_INSTALL_RPATH "")

# don't add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)



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
## plus the rpaths get all broken

# found this in the cmake docs
# install(TARGETS myExe mySharedLib myStaticLib
#          RUNTIME DESTINATION bin
#        LIBRARY DESTINATION lib
#          ARCHIVE DESTINATION lib/static)
#  install(TARGETS mySharedLib DESTINATION /some/full/path)
				

# INSTALL_PROGRAMS("~/local/bin" FILES src/estimator src/emulator)
# INSTALL_FILES("~/local/lib" FILES src/libEmu/libEmu.so src/libEmu/libLBFGS.so src/libRbind/libRBIND.so)
# INSTALL_FILES("~/local/include/libEmu/" FILES src/libEmu/emulate-fns.h src/libEmu/emulator.h src/libEmu/estimate_threaded.h src/libEmu/estimator.h src/libEmu/lbfgs.h src/libEmu/maxlbfgs.h src/libEmu/regression.h)
# INSTALL_FILES("~/local/include/libRbind/" FILES src/libRbind/rbind.h)


