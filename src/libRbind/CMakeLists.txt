# cmake for rbind lib

# set the include path for the headers for LBFGS and Emu
include_directories("../libEmu" "..")

IF(APPLE)
	SET(EXTRA_LIBS "-undefined dynamic_lookup -single_module -multiply_defined suppress -F/Library/Frameworks/R.framework/.. -framework R -Wl,-framework -Wl,CoreFoundation")
ENDIF(APPLE)


add_library(RBIND SHARED rbind ../bin_support)
target_link_libraries(RBIND Emu LBFGS ${libs} ${EXTRA_LIBS})									
add_dependencies(RBIND Emu LBFGS)

install (TARGETS RBIND
 				RUNTIME DESTINATION bin
 				LIBRARY DESTINATION lib)
