# cmake lists
# created by ccs: 15-may 2010

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")

get_filename_component(default_prefix ".." ABSOLUTE)
set(CMAKE_INSTALL_PREFIX $(default_prefix))

include_directories("./libEmu")

add_subdirectory("./libEmu")
add_subdirectory("./libRbind")

MESSAGE("link libs: ${libs}")

set(emulator_SOURCES emulator-main ioread useful optstruct resultstruct modelstruct bin_support)
set(estimator_SOURCES main ioread useful optstruct modelstruct resultstruct bin_support)

#set(lhood_SOURCES likelyhoody ioread useful)

link_directories(${PROJECT_SOURCE_DIR}/src/libEmu)
link_libraries(${LIBS} pthread)

add_executable(estimator ${estimator_SOURCES})
target_link_libraries(estimator Emu LBFGS ${libs})
add_dependencies(estimator Emu LBFGS)

add_executable(emulator ${emulator_SOURCES})
target_link_libraries(emulator Emu LBFGS ${libs})
add_dependencies(emulator Emu LBFGS)




