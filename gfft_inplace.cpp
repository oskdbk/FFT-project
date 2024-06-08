#include <iostream>
#include <fstream>
#include <sstream>
#include <climits>
#include <thread>
#include <numeric>
#include <iterator>
// #include <optional>
#include <vector>
#include <string>
#include <complex> // to use complex numbers
#include <cmath>
#include "gfft_inplace.h"
#include "utils.h"

typedef long double num;
typedef std::complex<num> complex_num;


using namespace std;

//-----------------------------------------------------------------------------


std::vector<complex_num> TransposeFFT_inplace(std::vector<complex_num> &P, int p, int q){
    std::vector<complex_num> transposed(q * p);

    for (int i = 0; i < p; ++i) {
        for (int j = 0; j < q; ++j) {
            int originalIndex = i * q + j;
            int transposedIndex = j * p + i;
            transposed[transposedIndex] = P[originalIndex];
        }
    }
    return transposed;
}


void TwistFFT_inplace(std::vector<complex_num> &P, int p, int q){
    int n = p * q;
    std::complex<long double> omega;
    long double angle;
    std::vector<std::vector<complex_num>> Q(p, std::vector<complex_num>(q));
    for (int i = 0; i < p; ++i) {
        
        for (int j = 0; j < q; ++j) {
            angle = 2 * M_PI * i * j / n;
            omega = std::complex<long double>(std::cos(angle), -std::sin(angle));
            P[i * q + j] *= omega;
        }
    }
}



vector<complex_num> GeneralFFT_inplace(vector<complex_num> &P, bool f){
    int n = P.size();
    if (n == 1){
        return vector<complex_num>{P[0]};
    }
    int p, q;
    std::tie(p, q) = Decompose(n);
    if(f == true){
        PRT1(P, " <- 0th step");
    }
    
    if(p == 1){
        return StandardFFT(P);
    }
    
    //1 step: transpose
    vector<complex_num> A(n);
    A=TransposeFFT_inplace(P, p, q);
    
    if(f == true){
        PRT1(A, " <- 1th step");
    }

    //2 step: FFT column-wise (after transpose row-wise)
    vector<complex_num> B;
    int first = 0;
    int last = p;
    vector<complex_num> temp;
    vector<complex_num> result;
    for(int k = 0; k < q; k++){
        temp.assign(A.begin() + first, A.begin() + last);
        result = GeneralFFT_inplace(temp, false);
        B.insert(B.end(), result.begin(), result.end());
        first += p;
        last += p;
    }
    if(f == true){
        PRT1(B, " <- 2th step");
    }


    A.clear();
    
    //3 step: Twiddle factors
    TwistFFT_inplace(B, q, p);
    
    if(f == true){
        PRT1(B, " <- 3th step");
    }

    //4 step: General transpose
    A = TransposeFFT_inplace(B, q, p);

    B.clear();

    if(f == true){
        PRT1(A, " <- 4th step");
    }
    first = 0;
    last = q;

    for(int k = 0; k < p; k++){
        temp.assign(A.begin() + first, A.begin() + last);
        result = GeneralFFT_inplace(temp, false);
        B.insert(B.end(), result.begin(), result.end());
        first += q;
        last += q;
    }

    A.clear();

    if(f == true){
        PRT1(B, " <- 5th step");
    }
    
    A = TransposeFFT_inplace(B, p, q);

    if(f == true){
        PRT1(A, " <- 6th step");
    }
    
    return A;
}
