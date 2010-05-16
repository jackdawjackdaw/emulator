cmake_minimum_required (VERSION 2.8)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")
project (Emulator)

SET(${project_name}_major_version 1)
SET(${project_name}_minor_version 0)

add_subdirectory(src)
