cmake_minimum_required (VERSION 2.6)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")
project (Emulator)

################### 
# the install prefix is where the bins and libs get dropped. 
#

message("Setting Install prefix to Home folder: $ENV{HOME}")
set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/local")


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


# make compilation verbose
#SET(CMAKE_VERBOSE_MAKEFILE TRUE)

find_package(GSL REQUIRED)
message("GSL_INCLUDE_DIRS: ${GSL_INCLUDE_DIR}")
message("CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
message("CMAKE_INSTALL_PATH: ${CMAKE_INSTALL_PREFIX}")

INCLUDE_DIRECTORIES(${GSL_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/libEmu)
MESSAGE("include dirs: ${CMAKE_include_directories}")
SET(libs ${libs} ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})

if(APPLE)
# use clang
SET(CMAKE_C_COMPILER /usr/bin/clang)
endif(APPLE)

# make a debug release
#set(CMAKE_BUILD_TYPE debug)

SET(CMAKE_C_FLAGS "-std=c99 -pg -O4")
SET(CMAKE_C_FLAGS_DEBUG "-std=c99 -pg -O0 -g -Wall -Wextra -pedantic")


MESSAGE("C Compiler: ${CMAKE_C_COMPILER}")

if(${CMAKE_BUILD_TYPE} MATCHES "debug")
 MESSAGE("C DEBUG FLAGS: ${CMAKE_C_FLAGS_DEBUG}")
else()
 MESSAGE("C FLAGS: ${CMAKE_C_FLAGS}")
endif(${CMAKE_BUILD_TYPE} MATCHES "debug")

add_subdirectory(src)

# explicitly install the headers we need, and the r-files				
INSTALL_FILES("/include/libEmu/" FILES src/libEmu/emulate-fns.h src/libEmu/emulator.h src/libEmu/estimate_threaded.h src/libEmu/estimator.h src/libEmu/lbfgs.h src/libEmu/maxlbfgs.h src/libEmu/regression.h)
# # ccs: added the R files we need to use rbind
# we have moved the contents of multivar.R into EmuRbind.R
INSTALL_FILES("/include/libRbind/" FILES src/libRbind/rbind.h src/libRbind/EmuRbind.R )



