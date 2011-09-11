# cmake for rbind lib

# set the include path for the headers for LBFGS and Emu
include_directories("../libEmu" "..")

IF(APPLE)

	SET(CMAKE_C_COMPILER "gcc-4.2")

	SET(CMAKE_SKIP_BUILD_RPATH  TRUE)
	SET(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE) 
	SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

	# this is cheating but at the end we export -DAPPLE
	# which will turn on the lines in the header we need that use the R framework
	SET(CMAKE_C_FLAGS "-arch x86_64 -std=gnu99 -I/Library/Frameworks/R.framework/Resources/include -I/Library/Frameworks/R.framework/Resources/include/x86_64 -I/usr/local/include -g -O2 -fPIC -DAPPLE")

	SET(EXTRA_LIBS "-undefined dynamic_lookup -single_module -multiply_defined suppress -F/Library/Frameworks/R.framework/.. -framework R -Wl,-framework -Wl,CoreFoundation")

ENDIF(APPLE)


add_library(RBIND SHARED rbind ../bin_support)
target_link_libraries(RBIND Emu LBFGS ${libs} ${EXTRA_LIBS})									
add_dependencies(RBIND Emu LBFGS)

install (TARGETS RBIND
 				RUNTIME DESTINATION bin
 				LIBRARY DESTINATION lib)
