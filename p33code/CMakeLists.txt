cmake_minimum_required(VERSION 3.8)
#project(p33)

FIND_PACKAGE(CUDA QUIET REQUIRED)
SET(CUDA_NVCC_FLAGS${CUDDA_NVCC_FLAGS};-std=c++11;)

include_directories(src)
cuda_add_executable(p33 src/p33.cu)

