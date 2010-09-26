# cmake for rbind lib

# set the include path for the headers for LBFGS and Emu
include_directories("../libEmu" "..")


add_library(RBIND SHARED rbind ../bin_support)
target_link_libraries(RBIND Emu LBFGS ${libs})									
add_dependencies(RBIND Emu LBFGS)
