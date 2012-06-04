# cmake lists for libEmu
# created by ccs: 15-may 2010
# modified for multimin oct-2011
cmake_minimum_required (VERSION 2.6)


include_directories(${CMAKE_SOURCE_DIR})

# using multimin
set(emudeps emulator.c estimate_threaded.c estimator-fns.c maxmultimin.c regression.c emulator.h estimate_threaded.h estimator-fns.h regression.h ../useful.c ../optstruct.c ../modelstruct.c ../resultstruct.c emulate-fns.c emulate-fns.h maxmultimin.h  ../emulator_struct.c ../emulator_struct.h)

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
