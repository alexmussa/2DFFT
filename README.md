# 2DFFT
Three codes that implement the Forward and Inverse 2DFFT in MPI, C++ Threads, and CUDA.

p31 - C++ Threads Implementation

p32 - MPI Implementation

p33 - CUDA

These codes use cmake version 3.9.1, mpicc 1.8, gcc 4.9, and CUDA to be installed. Make sure to load gcc/4.9.0, cmake/3.9.1, openmpi/1.8, and cuda/9.1 before following the instructions below.

INSTRUCTIONS:

1.) Download all of the files in the repo and create a build folder in the directory where the files are.

2.) Execute the command ```cmake ..``` from within the build folder. 

3.) Upon succesful creation of the make file, type ```make``` to build the project.

4.) The projects will all be built, with their titles as their executable names.

Example Execution:

```./p3x forward/reverse inputfile outputfile```

```./p33 forward Tower1024.txt CUDAoutput1024.txt```
