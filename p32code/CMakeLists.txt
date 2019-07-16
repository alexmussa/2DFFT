cmake_minimum_required(VERSION 3.8)
#project(p32)

SET(CMAKE_C_COMPILER mpicc)
SET(CMAKE_CXX_COMPILER mpicxx)
SET(CMAKE_CXX_FLAGS${CMAKE_CXX_FLAGS};-std=c++11;)

include_directories(src)
add_executable(p32 ${SOURCE_FILES} src/p32.cc src/input_image.h src/input_image.cc src/complex.h src/complex.cc)
