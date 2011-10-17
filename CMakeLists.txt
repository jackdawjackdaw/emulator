cmake_minimum_required (VERSION 2.6)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")
project (Emulator)


################### 
# the install prefix is where the bins and libs get dropped. 
message("Setting Install prefix to Home folder: $ENV{HOME}")
set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/local")


# lets get a version number from git
execute_process(
	#COMMAND git describe --tags --long
	COMMAND ../version.sh
	RESULT_VARIABLE v_result_var
	OUTPUT_VARIABLE v_output_var)

if(${v_result_var} MATCHES "0")
	SET(VERSION_NUMBER ${v_output_var})
	message("Building version: ${VERSION_NUMBER}")
else()
	SET(VERSION_NUMBER 1)
endif(${v_result_var} MATCHES "0")

configure_file (
	"${PROJECT_SOURCE_DIR}/src/buildConfig.h.in"
	"${PROJECT_SOURCE_DIR}/src/buildConfig.h"
)

include_directories("${PROJECT_BINARY_DIR}")

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


# make a debug release
#set(CMAKE_BUILD_TYPE debug)

if(APPLE)
	SET(CMAKE_C_FLAGS "-std=gnu99 -O2 -DAPPLE")
else()
	SET(CMAKE_C_FLAGS "-std=gnu99 -O2")
endif(APPLE)

SET(CMAKE_C_FLAGS_DEBUG "-std=gnu99 -O0 -g -pg -Wall -Wextra -pedantic")

if(APPLE)
# we need to use the same compiler as R CMD SHLIB would use
# or we're fucked
SET(CMAKE_C_COMPILER "gcc-4.2")
endif(APPLE)

MESSAGE("C Compiler: ${CMAKE_C_COMPILER}")

if(${CMAKE_BUILD_TYPE} MATCHES "debug")
 MESSAGE("C DEBUG FLAGS: ${CMAKE_C_FLAGS_DEBUG}")
else()
 MESSAGE("C FLAGS: ${CMAKE_C_FLAGS}")
endif(${CMAKE_BUILD_TYPE} MATCHES "debug")

add_subdirectory(src)

# explicitly install the headers we need, and the r-files				
INSTALL_FILES("/include/libEmu/" FILES src/libEmu/emulate-fns.h src/libEmu/emulator.h src/libEmu/estimate_threaded.h src/libEmu/estimator.h src/libEmu/regression.h)
# # ccs: added the R files we need to use rbind
# we have moved the contents of multivar.R into EmuRbind.R
INSTALL_FILES("/include/libRbind/" FILES src/libRbind/rbind.h src/libRbind/EmuRbind.R 
																	 src/libRbind/emuOverDesign.R src/libRbind/implausOverDesign.R
																	 src/libRbind/testEst.R)



