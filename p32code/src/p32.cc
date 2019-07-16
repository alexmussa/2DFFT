#include <mpi.h>
#include <iostream>
#include <string>
#include <array>
#include <vector>
#include <thread>
#include <cmath>
#include <chrono>
 
#include "complex.h"
#include "input_image.h"
 
const int numThreads = 8;
int im_width, im_height, RowsPerThread, count =0, numprocs = 0;
bool inverse;
 
Complex *im_data, *Wnk;

std::string outfilename;
 
using namespace std;
 
unsigned bitreversal(unsigned v, int N){
    int i = N;
    unsigned reversed = 0;
 
    for(--i; i>0;i >>= 1){
        reversed <<= 1;
        reversed |= (v & 0x1);
        v >>= 1;
    }
 
    return reversed;
}
 
void FFT1D(Complex* h, int im_width, bool inverse){
     //Daniel-Lanczos method for computing DFT.
    int Setsize, Nset, Nsetmember;
    Complex* h_tmp = new Complex [im_width];
 
    for(int i = 0; i < im_width; i++){
        *(h_tmp + i) = *(h + i);
    }

    for(int i = 0; i < im_width; i++){
        int reversedIdx = bitreversal(i, im_width);
        *(h + i) = *(h_tmp + reversedIdx);
    }

    for(Setsize = 2; Setsize <= im_width; Setsize = Setsize * 2){
        for(Nset = 0; Nset < (im_width/Setsize); Nset++){
            for(Nsetmember = 0; Nsetmember < (Setsize/2); Nsetmember++){
                Complex h_temp_fft = h[Setsize * Nset + Nsetmember];

                *(h + Setsize * Nset + Nsetmember)               = h_temp_fft + Wnk[Nsetmember * im_width / Setsize] * h[Setsize * Nset + Nsetmember + (Setsize / 2)];
                *(h + Setsize * Nset + Nsetmember + (Setsize/2)) = h_temp_fft - Wnk[Nsetmember * im_width / Setsize] * h[Setsize * Nset + Nsetmember + (Setsize / 2)];
            }
        }
    }
    
    if(inverse == true){
        for(int i = 0; i < im_width; i++){
            *(h + i) = (*(h + i)) * float(1/float(im_width));
        }
    }
}
 
void Transpose(Complex *h, int im_width, int im_height){
    Complex* h_tmp = new Complex[im_width*im_height];
    for(int row = 0; row < im_height; row++){
        for(int col = 0; col < im_width; col++){
            *(h_tmp + col + row*im_width) = *(h + col + row*im_width);
        }
    }
    for(int row = 0; row < im_height; row++){
        for(int col = 0; col < im_width; col++){
            *(h + col*im_width + row) = *(h_tmp + row*im_width + col);
        }
    }
}
 
void* threads(int threadIdx){
    for(int row = 0; row < RowsPerThread; row++){
        FFT1D((im_data + (threadIdx*RowsPerThread + row)*im_width), im_width, inverse);
    }
}

int main(int argc, char** argv){
    auto t1 = chrono::high_resolution_clock::now();

    MPI_Init(&argc, &argv);
    MPI_Status status;
    int rank, size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    numprocs = size;

    if(string(argv[1]) == "reverse"){
        inverse = true;
    }else if(string(argv[1]) == "forward"){
        inverse = false;
    }else{
        printf("Error: expecting argument 1 to be 'forward' or 'reverse'. Instead, received %s. Exitting...", argv[1]);
        exit(1);
    }

    string filename = string(argv[2]);
    outfilename = string(argv[3]);
    
    InputImage image = InputImage(filename.c_str());

    im_width = image.get_width();
    im_height = image.get_height();

    im_data = new Complex [im_width * im_height];
    im_data = image.get_image_data();

    RowsPerThread = im_height / size;

    Wnk = new Complex [im_width];
    
    if(inverse == false){
        for(int i = 0; i < im_width; i++){
            *(Wnk + i) = Complex(cos(2*M_PI*i/im_width), -1*sin(2*M_PI*i/im_width));
        }
    }else{
        for(int i = 0; i < im_width; i++){
            *(Wnk + i) = Complex(cos(2*M_PI*i/im_width), sin(2*M_PI*i/im_width));
        }
    }

    threads(rank);

    MPI_Barrier(MPI_COMM_WORLD);

    if(rank == 0){
        for(int i = 1; i<numprocs; i++){
            float buff_r[im_width*RowsPerThread];
            float buff_i[im_width*RowsPerThread];

            MPI_Recv(&buff_r,im_width*RowsPerThread,MPI_FLOAT,i,0, MPI_COMM_WORLD, &status);
            MPI_Recv(&buff_i,im_width*RowsPerThread,MPI_FLOAT,i,0, MPI_COMM_WORLD, &status);
 
            for(int j = 0; j < im_width*RowsPerThread; j++){
                *(im_data + i*RowsPerThread*im_width + j) = Complex(buff_r[j],buff_i[j]);
            }
        }
    }else{
        //Decompose complex result from each rank into array of floats for real and imaginary to send, then send to rank 0S.
        float buff_r[im_width*RowsPerThread];
        float buff_i[im_width*RowsPerThread];

        for(int i = 0; i < im_width*RowsPerThread; i++){
            *(buff_r + i) =  im_data[rank*im_width*RowsPerThread+i].real;
            *(buff_i + i) =  im_data[rank*im_width*RowsPerThread+i].imag;
        }

        MPI_Send(&buff_r,im_width*RowsPerThread,MPI_FLOAT,0,0,MPI_COMM_WORLD);
        MPI_Send(&buff_i,im_width*RowsPerThread,MPI_FLOAT,0,0,MPI_COMM_WORLD);
    }

    if(rank == 0){
        //Perform transpose on rank 0 and continue.
        Transpose(im_data,im_width,im_height);

        //Decompose complex into array of floats for real and imaginary to send, then send.
        float* im_buff_r = new float [im_width*im_height];
        float* im_buff_i = new float [im_width*im_height];

        for(int i = 0; i < im_width*im_height; i++){
            *(im_buff_r + i) = (*(im_data + i)).real;
            *(im_buff_i + i) = (*(im_data + i)).imag;
        }

        for(int i = 1; i < numprocs; i++){
            MPI_Send(&(*(im_buff_r)),im_width*im_height,MPI_FLOAT,i,0,MPI_COMM_WORLD);
            MPI_Send(&(*(im_buff_i)),im_width*im_height,MPI_FLOAT,i,0,MPI_COMM_WORLD);
        }
    }else{
        float* im_buff_r = new float [im_width*im_height];
        float* im_buff_i = new float [im_width*im_height];

        MPI_Recv(&(*(im_buff_r)),im_width*im_height,MPI_FLOAT,0,0, MPI_COMM_WORLD, &status);
        MPI_Recv(&(*(im_buff_i)),im_width*im_height,MPI_FLOAT,0,0, MPI_COMM_WORLD, &status);

        for(int i = 0; i < im_width*im_height; i++){
            *(im_data+i) = Complex(*(im_buff_r+i),*(im_buff_i+i)); 
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    threads(rank);

    if(rank == 0){
        for(int i = 1; i<numprocs; i++){
            float* buff_r = new float [im_width*RowsPerThread];
            float* buff_i = new float [im_width*RowsPerThread];

            MPI_Recv(&(*(buff_r)),im_width*RowsPerThread,MPI_FLOAT,i,0, MPI_COMM_WORLD, &status);
            MPI_Recv(&(*(buff_i)),im_width*RowsPerThread,MPI_FLOAT,i,0, MPI_COMM_WORLD, &status);
            
            for(int j = 0; j < im_width*RowsPerThread; j++){
                *(im_data + i*RowsPerThread*im_width + j) = Complex(buff_r[j],buff_i[j]);
            }
        }
    }else{
        float* buff_r = new float [im_width*RowsPerThread];
        float* buff_i = new float [im_width*RowsPerThread];
        for(int i = 0; i < im_width*RowsPerThread; i++){
            *(buff_r + i) =  im_data[rank*im_width*RowsPerThread+i].real;
            *(buff_i + i) =  im_data[rank*im_width*RowsPerThread+i].imag;
        }
        MPI_Send(&(*(buff_r)),im_width*RowsPerThread,MPI_FLOAT,0,0,MPI_COMM_WORLD);
        MPI_Send(&(*(buff_i)),im_width*RowsPerThread,MPI_FLOAT,0,0,MPI_COMM_WORLD);
    }

    if(rank == 0){
        Transpose(im_data,im_width,im_height);
        InputImage image = InputImage((string(argv[2]).c_str()));
        
        if(inverse == false){
            image.save_image_data(outfilename.c_str(),im_data,im_width,im_height);
        }else{
            image.save_image_data_real(outfilename.c_str(),im_data,im_width,im_height);
        }
    }

    MPI_Finalize();

    auto t2 = chrono::high_resolution_clock::now();
    chrono::duration<double> duration = (t2-t1);
    cout << duration.count() << '\n';
}