# cmake lists
# created by ccs: 15-may 2010
cmake_minimum_required(VERSION 2.8)


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")

ADD_SUBDIRECTORY(libEmu)

ADD_DEFINITIONS("-std=c99 -pg -O2 -Wextra -pedantic")

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS} ${GSLCBLAS_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/libEmu)
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
