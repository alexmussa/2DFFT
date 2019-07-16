//Allie Alexander I. Mussa - ECE 6122 - Final Project - Fall '18
//2D-DFT - Normal - CUDA Implementation (p33)
//Last Edit - 12/13/18 09:30 AM

#include <iostream>
#include <cmath>
#include <chrono>
#include "complex.h"
#include "input_image.h"
#define T_P_B 1024

int im_width = 0;
int im_height = 0;
int RowsPerThread = 0;
bool inverse;
 
Complex *im_data, *im_data_tps, *im_transformed;
Complex *Wnk;

using namespace std;

__global__ void DFT1D(Complex* data, Complex* Wnk,Complex* data_fft, const int im_height, const int im_width, const bool inverse){

    int index = threadIdx.x + blockIdx.x * blockDim.x;
    
    if(index < im_width*im_height){
        data_fft[index] = 0;
        int factor = 1;

        if(im_width > blockDim.x){
            factor = im_width/blockDim.x;
            for(int i = 0; i < im_width; i++){
                *(data_fft + index) = *(data_fft + index) + (*(data + i + (blockIdx.x/factor)*im_width + int(threadIdx.x/im_width)*im_width) * *(Wnk + i + (index%im_width)*im_width));
            }
        }else{
            for(int i = 0; i < im_width; i++){
                *(data_fft + index) = *(data_fft + index) + (*(data + i + (blockIdx.x)*T_P_B + int(threadIdx.x/im_width)*im_width) * *(Wnk + i + (index%im_width)*im_width));
            }
        }

        

        if(inverse == true){
            *(data_fft + index) = *(data_fft + index) * float(1/float(im_width)); 
        }
    }
}

__global__ void transpose(Complex* data, Complex* data_transpose, int im_width, int im_height){
    int index = threadIdx.x + blockIdx.x * blockDim.x;

    if(index < im_width*im_height){
        *(data_transpose + index) = *(data + int(index/im_width) + (index%im_width)*im_width);
    }
}

int main(int argc, char **argv){
    auto t1 = chrono::high_resolution_clock::now();
    if(string(argv[1]) != string("forward") && string(argv[1]) != string("reverse")){
        printf("Error: expecting argument 1 to be 'forward' or 'reverse'. Instead, received %s. Exitting...", argv[1]);
        exit(1);
    }
 
    if(string(argv[1]) == "reverse"){
        inverse = true;
    }else if(string(argv[1]) == "forward"){
        inverse = false;
    }
 
    string filename = string(argv[2]);
    string outfilename = string(argv[3]);
    
    InputImage image(filename.c_str());

    im_width = image.get_width();
    im_height = image.get_height();
 
    im_data = new Complex [im_width * im_height];
    im_data = image.get_image_data();

    Complex *d_Wnk, *d_im_data, *d_im_data_ft;
    Complex *im_data_rec = new Complex [im_width*im_height];

    Wnk = new Complex [im_width*im_width];
    
    if(inverse == false){
        for(int n = 0; n < im_width; n++){
            for(int k = 0; k < im_width; k++){
                *(Wnk + n*im_width + k) = Complex(cos(2*M_PI*k*n/im_width), -1*sin(2*M_PI*k*n/im_width));
            }
        }
    }else{
        for(int n = 0; n < im_width; n++){
            for(int k = 0; k < im_width; k++){
                *(Wnk + n*im_width + k) = Complex(cos(2*M_PI*k*n/im_width), sin(2*M_PI*k*n/im_width));
            }
        }
    }

    int im_size = im_width * im_height;

    cudaMalloc((void**)&d_Wnk, im_width*im_width*sizeof(Complex));
    cudaMalloc((void**)&d_im_data, im_width*im_height*sizeof(Complex));
    cudaMalloc((void**)&d_im_data_ft, im_width*im_height*sizeof(Complex));

    cudaMemcpy(d_Wnk, Wnk, im_width*im_width*sizeof(Complex), cudaMemcpyHostToDevice);
    cudaMemcpy(d_im_data, im_data, im_width*im_height*sizeof(Complex), cudaMemcpyHostToDevice);

    DFT1D<<<(im_size + T_P_B -1)/T_P_B,T_P_B>>>(d_im_data,d_Wnk,d_im_data_ft,im_height,im_width,inverse);
    cudaDeviceSynchronize();
    transpose<<<(im_size + T_P_B -1)/T_P_B,T_P_B>>>(d_im_data_ft,d_im_data, im_width, im_height);
    cudaDeviceSynchronize();
    DFT1D<<<(im_size + T_P_B -1)/T_P_B,T_P_B>>>(d_im_data,d_Wnk,d_im_data_ft,im_height,im_width,inverse);
    cudaDeviceSynchronize();
    transpose<<<(im_size + T_P_B -1)/T_P_B,T_P_B>>>(d_im_data_ft,d_im_data, im_width, im_height);

    cudaMemcpy(im_data_rec, d_im_data, im_width*im_height*sizeof(Complex), cudaMemcpyDeviceToHost);
    
    if(inverse == false){
        image.save_image_data(outfilename.c_str(),im_data_rec,im_width,im_height);
    }else{
        image.save_image_data_real(outfilename.c_str(),im_data_rec,im_width,im_height);
    }
    
    cudaFree(d_Wnk);
    cudaFree(d_im_data);
    cudaFree(d_im_data_ft);

    auto t2 = chrono::high_resolution_clock::now();
    chrono::duration<double> duration = (t2-t1);
    cout << duration.count() << '\n';
}