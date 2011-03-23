cmake_minimum_required (VERSION 2.6)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")
project (Emulator)

################### 
# the install prefix is where the bins and libs get dropped. 
#
set(CMAKE_INSTALL_PREFIX "/home/einstein/grad/cec24/local")
if(APPLE)
	set(CMAKE_INSTALL_PREFIX "/Users/chris/local")
endif(APPLE)										 


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

# use clang
#SET(CMAKE_C_COMPILER /usr/bin/clang)


# make a debug release
set(CMAKE_BUILD_TYPE debug)

SET(CMAKE_C_FLAGS "-std=c99 -pg -O3")
SET(CMAKE_C_FLAGS_DEBUG "-std=c99 -pg -O0 -g -Wall -Wextra -pedantic")


MESSAGE("C Compiler: ${CMAKE_C_COMPILER}")

if(${CMAKE_BUILD_TYPE} MATCHES "debug")
 MESSAGE("C DEBUG FLAGS: ${CMAKE_C_FLAGS_DEBUG}")
else()
 MESSAGE("C FLAGS: ${CMAKE_C_FLAGS}")
endif(${CMAKE_BUILD_TYPE} MATCHES "debug")

add_subdirectory(src)


## do all the installations at the target level, saves worrying about the
## actual file names etc

#INSTALL_PROGRAMS("/bin" FILES src/estimator src/emulator)
# on os-x we get .dylib instead of .so
# there is certainly a better way to do this
# if(APPLE)
# 	INSTALL_FILES("/lib" FILES src/libEmu/libEmu.dylib src/libEmu/libLBFGS.dylib src/libRbind/libRBIND.dylib)
# else()
# 	INSTALL_FILES("/lib" FILES src/libEmu/libEmu.so src/libEmu/libLBFGS.so src/libRbind/libRBIND.so)				
# endif(APPLE)
				
# INSTALL_FILES("/include/libEmu/" FILES src/libEmu/emulate-fns.h src/libEmu/emulator.h src/libEmu/estimate_threaded.h src/libEmu/estimator.h src/libEmu/lbfgs.h src/libEmu/maxlbfgs.h src/libEmu/regression.h)
# INSTALL_FILES("/include/libRbind/" FILES src/libRbind/rbind.h)


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



