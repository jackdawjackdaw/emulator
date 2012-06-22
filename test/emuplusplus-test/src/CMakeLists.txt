# cmake lists for emuplusplus-test
# created by ccs: 22.06.2012
cmake_minimum_required (VERSION 2.6)

link_directories(${CMAKE_INSTALL_PREFIX}/lib)


set(plusTest_SOURCES example)
add_executable(test_Emu++ ${plusTest_SOURCES})
target_link_libraries(test_Emu++ EmuPlusPlus ${libs} pthread)
add_dependencies(test_Emu++ EmuPlusPlus)

install(TARGETS test_Emu++ 
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib)
