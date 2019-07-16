#include <iostream>
#include <string>
#include <array>
#include <vector>
#include <thread>
#include <cmath>
#include <chrono>
 
#include "complex.h"
#include "input_image.h"
 
int numThreads = 8; 
int im_width = 0;
int im_height = 0;
int RowsPerThread = 0;
bool inverse;
 
Complex *im_data, *im_data_tps, *im_transformed;
Complex *Wnk;
 
using namespace std;
 
void DFT1D(Complex* h, int im_width, bool inverse){
    Complex h_tmp[im_width], sum;
    for(int i = 0; i < im_width; i++){
        *(h_tmp + i) = *(h + i);
    }

    
    for(int i = 0; i < im_width; i++){
        sum = 0;
        for(int j = 0; j < im_width; j++){
            sum = sum + h_tmp[j]*Wnk[i*im_width + j];
        }
        *(h+i) = sum;
    }

    if(inverse == true){
        for(int i = 0; i < im_width; i++){
            *(h + i) = *(h + i) * float(1/float(im_width));
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
        DFT1D((im_data + (threadIdx*RowsPerThread + row)*im_width), im_width, inverse);
    }
}

int main(int argc, char** argv){
    
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
 
    RowsPerThread = im_height / numThreads;

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

    vector<thread> vecOfThreads;
    for(int threadIdx = 0; threadIdx < numThreads; threadIdx++){
        vecOfThreads.push_back(thread(threads,threadIdx));
    }
    for(auto &t : vecOfThreads){
        t.join();
    }

    Transpose(im_data,im_width,im_height);

    vector<thread> vecOfThreads2;
    for(int threadIdx = 0; threadIdx < numThreads; threadIdx++){
        vecOfThreads2.push_back(thread(threads,threadIdx));
    }
    for(auto &t : vecOfThreads2){
        t.join();
    }

    Transpose(im_data,im_width,im_height);

    if(inverse == false){
            image.save_image_data(outfilename.c_str(),im_data,im_width,im_height);
    }else{
        image.save_image_data_real(outfilename.c_str(),im_data,im_width,im_height);
    }

    auto t2 = chrono::high_resolution_clock::now();
    chrono::duration<double> duration = (t2-t1);
    cout << duration.count() << '\n';

}