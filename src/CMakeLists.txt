# cmake lists
# created by ccs: 15-may 2010
# doesn't work well on os-x laptop
# there's a problem with the gsl includes
cmake_minimum_required(VERSION 2.8)


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")

get_filename_component(default_prefix ".." ABSOLUTE)
set(CMAKE_INSTALL_PREFIX $(default_prefix))

ADD_SUBDIRECTORY(libEmu)



find_package(GSL REQUIRED)

message(GSL_INCLUDE_DIR="${GSL_INCLUDE_DIR}")
message(GSL_CXX_FLAGS="${GSL_CXX_FLAGS}")

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}-std=c99 -pg -O2 -Wextra -pedantic ${GSL_CXX_FLAGS}")

INCLUDE_DIRECTORIES(${GSL_INCLUDE_DIR})
set(libs ${libs} ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})

set(emulator_SOURCES emulator-main ioread useful)
set(estimator_SOURCES main ioread useful)
set(lhood_SOURCES likelyhoody ioread useful)

link_directories(${PROJECT_BINARY_DIR}/libEmu)
link_libraries(${LIBS} pthread)

add_executable(estimator ${estimator_SOURCES})
target_link_libraries(estimator Emu LBFGS ${libs})
add_dependencies(estimator Emu LBFGS)

add_executable(emulator ${emulator_SOURCES})
target_link_libraries(emulator Emu LBFGS ${libs})
add_dependencies(emulator Emu LBFGS)

add_executable(lhood ${lhood_SOURCES})
target_link_libraries(lhood Emu LBFGS ${libs})
add_dependencies(lhood Emu LBFGS)

install(TARGETS emulator estimator RUNTIME DESTINATION bin LIBRARY DESTINATION LIB)
