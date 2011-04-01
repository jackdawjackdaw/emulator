# cmake lists for libEmu
# created by ccs: 15-may 2010
cmake_minimum_required (VERSION 2.6)
enable_language (Fortran)

include_directories(${CMAKE_SOURCE_DIR})
set(emudeps emulator.c estimate_threaded.c estimator.c lbfgs.c maxlbfgs.c regression.c emulator.h estimate_threaded.h estimator.h lbfgs.h  maxlbfgs.h regression.h ../useful.c ../optstruct.c ../modelstruct.c ../resultstruct.c emulate-fns.c emulate-fns.h)

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME STREQUAL "gfortran")
  # gfortran
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g")
elseif (Fortran_COMPILER_NAME STREQUAL "ifort")
  # ifort (untested)
  set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g")
elseif (Fortran_COMPILER_NAME STREQUAL "g77")
  # g77
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -m32")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -m32")
else (Fortran_COMPILER_NAME STREQUAL "gfortran")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME STREQUAL "gfortran")

ADD_LIBRARY(LBFGS SHARED routines.f)

ADD_LIBRARY(Emu SHARED ${emudeps})
TARGET_LINK_LIBRARIES(Emu LBFGS ${libs})
add_dependencies(Emu LBFGS)


install(TARGETS Emu LBFGS 
 				RUNTIME DESTINATION bin 
 				LIBRARY DESTINATION lib)
