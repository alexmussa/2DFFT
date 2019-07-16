# 2DFFT
Three codes that implement the Forward and Inverse 2DFFT in MPI, C++ Threads, and CUDA.

p31 - C++ Threads Implementation

p32 - MPI Implementation

p33 - CUDA

These codes use cmake version 3.8, mpicc, gcc 7.3, and CUDA to be installed.

INSTRUCTIONS:

1.) Download all of the files in the repo and create a build folder in the directory where the files are.

2.) Execute the command ```cmake ..``` from within the build folder. 

3.) Upon succesful creation of the make file, type ```make``` to build the project.

4.) The projects will all be built, with their titles as their executable names.

Example Execution:

```./p3x forward/reverse inputfile outputfile```

```./p33 forward Tower1024.txt CUDAoutput1024.txt```
