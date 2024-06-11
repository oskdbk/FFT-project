// fft_utils.cpp
#include "gfft.h"  // Include the header file for function declarations and type definitions
#include "../utils.h" // Include the header file for utility functions
#include <iostream>     // For std::cout, std::cin, etc.
#include <fstream>      // For file handling
#include <sstream>      // For string stream operations
#include <cmath>        // For mathematical functions like std::cos, std::sin, etc.
#include <string>       // For std::string
#include <vector>       // For std::vector
#include <complex>      // For std::complex
#include <tuple>        // For std::tuple

// Use the standard namespace
using namespace std;

//-----------------------------------------------------------------------------

std::vector<std::vector<complex_num>> TransposeFFT(std::vector<std::vector<complex_num>> &P){
    size_t p = P.size();
    size_t q = P[0].size();
    std::vector<std::vector<complex_num>> T(q, std::vector<complex_num>(p));
    for (size_t i = 0; i < p; ++i) {
        for (size_t j = 0; j < q; ++j) {
            T[j][i] = P[i][j]; 
        }
    }
    return T;
}

std::vector<std::vector<complex_num>> TwistFFT(std::vector<std::vector<complex_num>> &P){
    size_t p = P.size();
    size_t q = P[0].size();
    size_t n = p * q;
    std::complex<double> omega;
    double angle;
    std::vector<std::vector<complex_num>> Q(p, std::vector<complex_num>(q));
    for (size_t i = 0; i < p; ++i) {
        
        for (size_t j = 0; j < q; ++j) {
            angle = 2 * M_PI * i * j / n;
            std::complex<double> omega(std::cos(angle), -std::sin(angle));
            Q[i][j] = P[i][j] * omega; 
        }
    }
    return Q;
}

std::vector<std::vector<complex_num>> PackFFT(vector<complex_num> &P, int p, int q){
    std::vector<std::vector<complex_num>> A(p, std::vector<complex_num>(q));
    for(int j = 0; j < p; j++){
        for(int k = 0; k < q; k++){
            A[j][k] = P[j * q + k];
        }
    }
    return A;
}

vector<complex_num> UnPackFFT(std::vector<std::vector<complex_num>> &P, int p, int q){

    std::vector<std::complex<double>> output(p * q, 0);
    for(int j = 0; j < p; j++){
        for(int k = 0; k < q; k++){
            output[j * q + k] = P[j][k];
        }
    }
    return output;
}


vector<complex_num> GeneralFFT(vector<complex_num> P, bool f){
    size_t n = P.size();
    if (n == 1){
        return vector<complex_num>{P[0]};
    }
    int p, q;
    std::tie(p, q) = Decompose(n);
    if(f == true){
        PRT1(P, "P");
    }
    
    if(p == 1){
        return StandardFFT(P);
    }
    std::vector<std::vector<complex_num>> A(p, std::vector<complex_num>(q));
    A = PackFFT(P, p, q);
    
    if(f == true){
        PRT2(A, "A");
    }


    // transpose A
    std::vector<std::vector<complex_num>> A1(q, std::vector<complex_num>(p));
    A1 = TransposeFFT(A);
    if(f == true){
        PRT2(A1, "A1");
    }

    // FFT on each row of A1
    std::vector<std::vector<complex_num>> B(q, std::vector<complex_num>(p));
    for(int k = 0; k < q; k++){
        B[k] = GeneralFFT(A1[k], false);
    }
    if(f == true){
        PRT2(B, "B");
    }
   
    // Twiddle factor multiplication
    std::vector<std::vector<complex_num>> C(q, std::vector<complex_num>(p));
    C = TwistFFT(B);
    if(f == true){
        PRT2(C, "C");
    }
    

    // transpose C
    std::vector<std::vector<complex_num>> D(p, std::vector<complex_num>(q));
    D = TransposeFFT(C);
    if(f == true){
        PRT2(D, "D");
    }
    
    // FFT on each row of D
    std::vector<std::vector<complex_num>> E(p, std::vector<complex_num>(q));
    for(int k = 0; k < p; k++){
        E[k] = GeneralFFT(D[k], false);
    }
    if(f == true){
        PRT2(E, "E");
    }
    
    // Transpose E
    std::vector<std::vector<complex_num>> F(q, std::vector<complex_num>(p));
    F = TransposeFFT(E);
    if(f == true){
        PRT2(F, "F");
    }
    
    // Unpack F
    vector<complex_num> G(n);
    G = UnPackFFT(F, q, p);
    if(f == true){
        PRT1(G, "G");
    }
    
    return G;
}




