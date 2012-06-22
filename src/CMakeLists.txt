# cmake lists
# created by ccs: 15-may 2010
cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")

get_filename_component(default_prefix ".." ABSOLUTE)
set(CMAKE_INSTALL_PREFIX $(default_prefix))

add_subdirectory("./libEmu")
add_subdirectory("./libRbind")

MESSAGE("link libs: ${libs}")

set(common_SOURCES ioread useful optstruct resultstruct modelstruct multi_modelstruct multivar_support bin_support)
set(emulator_SOURCES emulator ${common_SOURCES}) 
set(estimator_SOURCES estimator ${common_SOURCES})
set(interactive_SOURCES interactive_emulator ${common_SOURCES} emulator_struct)

#set(lhood_SOURCES likelyhoody ioread useful)
# on my deb6 vmware machine we get an error about cannot 
# find pthread_create, pthread should be explicitly linked?
# 
link_directories(${PROJECT_SOURCE_DIR}/src/libEmu)
link_libraries(${LIBS} pthread)

# add_executable(estimator ${estimator_SOURCES})
# target_link_libraries(estimator Emu ${libs} pthread)
# add_dependencies(estimator Emu)

# add_executable(emulator ${emulator_SOURCES})
# target_link_libraries(emulator Emu ${libs} pthread)
# add_dependencies(emulator Emu)

add_executable(interactive_emulator ${interactive_SOURCES})
target_link_libraries(interactive_emulator Emu ${libs} pthread)
add_dependencies(interactive_emulator Emu)

## the emuplusplus library
set(emuplusplus_SOURCES EmuPlusPlus ${common_SOURCES} emulator_struct)
add_library(EmuPlusPlus SHARED ${emuplusplus_SOURCES})
target_link_libraries(EmuPlusPlus Emu ${libs} pthread)
add_dependencies(EmuPlusPlus Emu)

install(TARGETS EmuPlusPlus
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib)

install(FILES EmuPlusPlus.h
	DESTINATION include)


install (TARGETS interactive_emulator
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib)

