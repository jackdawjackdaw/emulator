# cmake lists for libEmu
# created by ccs: 15-may 2010
# modified for multimin oct-2011
cmake_minimum_required (VERSION 2.6)
#enable_language (Fortran)

include_directories(${CMAKE_SOURCE_DIR})

# using multimin
set(emudeps emulator.c estimate_threaded.c estimator-fns.c maxmultimin.c regression.c emulator.h estimate_threaded.h estimator-fns.h regression.h ../useful.c ../optstruct.c ../modelstruct.c ../resultstruct.c emulate-fns.c emulate-fns.h maxmultimin.h )

# we don't need this on the multimin branch
# FFLAGS depend on the compiler
# get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

# if (Fortran_COMPILER_NAME STREQUAL "gfortran")
#   # gfortran
#   set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3")
#   set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g")
# elseif (Fortran_COMPILER_NAME STREQUAL "ifort")
#   # ifort (untested)
#   set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3")
#   set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g")
# elseif (Fortran_COMPILER_NAME STREQUAL "g77")
#   # g77
#   set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -m32")
#   set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -m32")
# else (Fortran_COMPILER_NAME STREQUAL "gfortran")
#   message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
#   message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
#   message ("No optimized Fortran compiler flags are known, we just try -O2...")
#   set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
#   set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
# endif (Fortran_COMPILER_NAME STREQUAL "gfortran")

# ADD_LIBRARY(LBFGS SHARED routines.f)


IF(APPLE)

	SET(CMAKE_C_COMPILER "gcc")

	SET(CMAKE_SKIP_BUILD_RPATH  TRUE)
	SET(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE) 
	SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

	# this is cheating but at the end we export -DAPPLE
	# which will turn on the lines in the header we need that use the R framework
	SET(CMAKE_C_FLAGS "-g -O2 -fPIC")

	SET(EXTRA_LIBS "-undefined dynamic_lookup -single_module -multiply_defined suppress -Wl,-framework -Wl,CoreFoundation")
ENDIF(APPLE)



ADD_LIBRARY(Emu SHARED ${emudeps})
# this seems to be important
TARGET_LINK_LIBRARIES(Emu ${libs} pthread ${EXTRA_LIBS})
#add_dependencies(Emu )

install(TARGETS Emu  
 				RUNTIME DESTINATION bin 
 				LIBRARY DESTINATION lib)
