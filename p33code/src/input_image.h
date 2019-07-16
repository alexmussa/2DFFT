//
// Created by brian on 11/20/18.
//

#pragma once
#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>

#include "complex.h"

class Complex;

class InputImage {
public:

    InputImage(const char* filename);
    int get_width() const;
    int get_height() const;

    //returns a pointer to the image data.  Note the return is a 1D
    //array which represents a 2D image.  The data for row 1 is
    //immediately following the data for row 0 in the 1D array
    Complex* get_image_data() const;

    //use this to save output from forward DFT
    void save_image_data(const char* filename, Complex* d, int w, int h);
    //use this to save output from reverse DFT
    void save_image_data_real(const char* filename, Complex* d, int w, int h);

private:
    int w;
    int h;
    Complex* data;
};

InputImage::InputImage(const char* filename) {
    std::ifstream ifs(filename);
    if(!ifs) {
        std::cout << "Can't open image file " << filename << std::endl;
        exit(1);
    }

    ifs >> w >> h;

    data = new Complex[w * h];
    for(int r = 0; r < h; ++r) {
        for(int c = 0; c < w; ++c) {
            std::string temp;
            ifs >> temp;
            float real;
            float imag;
            int count=0, count_im=0;
            bool flag = false;
            
            std::string token, temp_r, temp_i;
            
            if(temp[0] == char('(')){
                //complex
                for(int i = 1; i < temp.length(); i++){
                    if(temp[i] != char(',') && flag == false){
                        temp_r[i-1] = temp[i];
                        count++;

                    }else if(temp[i] == char(',') && flag == false){
                        flag = true;
                        real = stof(temp.substr(1,count));

                    }else if(temp[i] != char(')') && flag == true){
                        temp_i[i - (count + 2)] = temp[i];
                        count_im++;

                    }else if(temp[i] == char(')') && flag == true){
                        imag = stof(temp.substr(count+2,count_im));
                    }
                }
                data[r * w + c] = Complex(real,imag);
            }else{
                //real
                real = stof(temp);
                data[r * w + c] = Complex(real);
            }
        }
    }
}  

int InputImage::get_width() const {
    return w;
}

int InputImage::get_height() const {
    return h;
}

Complex* InputImage::get_image_data() const {
    return data;
}

void InputImage::save_image_data(const char *filename, Complex *d, int w, int h) {
    std::ofstream ofs(filename);
    if(!ofs) {
        std::cout << "Can't create output image " << filename << std::endl;
        return;
    }

    ofs << w << " " << h << std::endl;

    for(int r = 0; r < h; ++r) {
        for(int c = 0; c < w; ++c) {
            ofs << d[r * w + c] << " ";
        }
        ofs << std::endl;
    }
}

void InputImage::save_image_data_real(const char* filename, Complex* d, int w, int h) {
    std::ofstream ofs(filename);
    if(!ofs) {
        std::cout << "Can't create output image " << filename << std::endl;
        return;
    }

    ofs << w << " " << h << std::endl;

    for (int r = 0; r < h; ++r) {
        for (int c = 0; c < w; ++c) {
            ofs << d[r * w + c].real << " ";
        }
        ofs << std::endl;
    }
}
