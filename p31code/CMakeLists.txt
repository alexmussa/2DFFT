cmake_minimum_required(VERSION 3.8)
#project(p31)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -Werror")

include_directories(
        src
)

set(SOURCE_FILES
        src/complex.h
	src/complex.cc
        src/input_image.h
        src/input_image.cc
        )

add_executable(p31 ${SOURCE_FILES} src/p31.cc)
