//
// Created by brian on 11/20/18.
//
#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __device__ __host__
#else
#define CUDA_CALLABLE_MEMBER
#endif

#pragma once

#include <iostream>
#include <cmath>

class Complex {
    public:
        Complex();
        Complex(float r, float i);
        Complex(float r);
        Complex operator+(const Complex& b) const;
        Complex operator-(const Complex& b) const;
        Complex operator*(const Complex& b) const;

        Complex mag() const;
        Complex angle() const;
        Complex conj() const;

        float real;
        float imag;
};

const float PI = 3.14159265358979f;

CUDA_CALLABLE_MEMBER Complex::Complex() : real(0.0f), imag(0.0f) {}

CUDA_CALLABLE_MEMBER Complex::Complex(float r) : real(r), imag(0.0f) {}

CUDA_CALLABLE_MEMBER Complex::Complex(float r, float i) : real(r), imag(i) {}

CUDA_CALLABLE_MEMBER Complex Complex::operator+(const Complex &b) const {
    return Complex(real + b.real, imag + b.imag);
}

CUDA_CALLABLE_MEMBER Complex Complex::operator-(const Complex &b) const {
    return Complex(real - b.real, imag - b.imag);
}

CUDA_CALLABLE_MEMBER Complex Complex::operator*(const Complex &b) const {
    return Complex(real*b.real - imag*b.imag , real*b.imag + imag*b.real);
}

CUDA_CALLABLE_MEMBER Complex Complex::mag() const {
    return Complex(sqrt(pow(real,2) + pow(imag,2)));
}

CUDA_CALLABLE_MEMBER Complex Complex::angle() const {
    return Complex(atan2(imag, real) * (180/PI));
}

CUDA_CALLABLE_MEMBER Complex Complex::conj() const {
    return Complex(real,-imag);
}

std::ostream& operator<< (std::ostream& os, const Complex& rhs) {
    Complex c(rhs);
    if(fabsf(rhs.imag) < 1e-10) c.imag = 0.0f;
    if(fabsf(rhs.real) < 1e-10) c.real = 0.0f;

    if(c.imag == 0) {
        os << c.real;
    }
    else {
        os << "(" << c.real << "," << c.imag << ")";
    }
    return os;
}