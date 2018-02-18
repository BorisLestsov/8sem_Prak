cmake_minimum_required(VERSION 3.5)
project(Prak)

set(CMAKE_CXX_COMPILER mpicxx)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

include_directories(include)

set(SOURCE_FILES source/main.cpp )
add_executable(Prak ${SOURCE_FILES})

set_target_properties(Prak PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY_DEBUG   ${CMAKE_SOURCE_DIR}/build/bin
        RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_SOURCE_DIR}/build/bin
        )


add_custom_target(make COMMAND $(MAKE) -C ${CMAKE_SOURCE_DIR})
