
#include <vector>
#include <string>
#include <cmath>
#include <numeric>
#include <complex> 
#include <thread>
#include "gfft.h"
#include "gfft_parallel.h"
#include "gfft_inplace.h"
#include "utils.h"
using namespace std;

// ------------------parallel versions----------------- //

void FFTThread(vector<complex_num> &P, int thread_id, int num_threads, int p, int q, bool inverse = false){
    // input: matrix P of dimensions p x q
    // in all columns of P that are = thread_id mod num_threads,
    // replace the column by its FFT
    int n = P.size();
    for(int col = thread_id; col < q; col += num_threads){
        vector<complex_num> A(p);
        for(int i = 0; i < p; i++){
            A[i] = P[i * q + col];
        }
        // A = GeneralFFT_inplace(A, false);
        A = DFT(A, inverse);
        for(int i = 0; i < p; i++){
            P[i * q + col] = A[i];
        }
    }
}




void TwiddleThread(vector<complex_num> &P, int thread_id, int num_threads, int p, int q, bool inverse = false){
    // input: matrix P of dimensions p x q
    // for all elements of P that are = thread_id mod num_threads,
    // multiply the element by the twiddle factor
    int n = P.size();
    std::complex<double> omega;
    double angle;
    for(int i = thread_id; i < n; i += num_threads){ // i is the raw index of the element
        
        int row = i / q; // row index of the element
        int col = i % q; // column index of the element
        // if(i<20 && thread_id == 1){
        //     cout<<i<<" ("<<row<<","<<col<<"): " << P[i];
        // }
        angle = 2 * M_PI * row * col / n;
        if (inverse) angle = -angle;
        omega = std::complex<double>(std::cos(angle), -std::sin(angle));
        complex<double> temp = P[i] * omega;
        P[i] = temp;
        // if(i<20 && thread_id == 1){
        //     cout<<" * " << omega << " = "<< P[i]<< " ("<< temp<<")"<<endl;
        // }
    }
}


void TransposeThread(vector<complex_num> &P, vector<complex_num> &Q, int thread_id, int num_threads, int p, int q){
    // transpose the pxq matrix P into the qxp matrix Q
    // do all indices of Q that are = thread_id mod num_threads
    int n = P.size();

    // i is the raw index of the element in Q
    for (int i = thread_id; i < n; i+= num_threads)
    {   
        int row_Q = i / p;
        int col_Q = i % p;
        int row = col_Q;
        int col = row_Q;
        Q[row_Q * p + col_Q] = P[row * q + col];
    }   
}

// --------------------------------------------------------------------- //

// ------------------sequential versions for debugging----------------- //
vector<complex_num> FFT(vector<complex_num> &P, int p, int q){
    // input: matrix P of dimensions p x q
    // in all columns of P, replace the column by its FFT
    
    int n = P.size();
    vector<complex_num> P1(n);
    for(int col = 0; col < q; col++){
        vector<complex_num> A(p);
        for(int i = 0; i < p; i++){
            A[i] = P[i * q + col];
        }
        A = GeneralFFT(A);
        for(int i = 0; i < p; i++){
            P1[i * q + col] = A[i];
        }
    }
    return P1;
}

void Transpose(vector<complex_num> &P, vector<complex_num> &Q, int p, int q){
    // transpose the pxq matrix P into the qxp matrix Q
    int n = P.size();
    for(int i = 0; i < n; i++){
        int row_Q = i / p;
        int col_Q = i % p;
        int row = col_Q;
        int col = row_Q;
        Q[row_Q * p + col_Q] = P[row * q + col];
    }
}

vector<complex_num> Twiddle(vector<complex_num> &P, int p, int q){
    // input: matrix P of dimensions p x q
    // for all elements of P, multiply the element by the twiddle factor
    
    int n = P.size();
    vector<complex_num> P1(n);
    std::complex<double> omega;
    double angle;
    for(int i = 0; i < n; i++){ // i is the raw index of the element
        int row = i / q; // row index of the element
        int col = i % q; // column index of the element
        angle = 2 * M_PI * row * col / n;
        omega = std::complex<double>(std::cos(angle), -std::sin(angle));
        P1[i] = P[i] * omega;
        // if(i<20)
        //     cout << P[i] << " * " << omega << " = " << P1[i] << endl;
    }
    return P1;
}

// --------------------------------------------------------------------- //




vector<complex_num> GeneralFFT_Parallel(vector<complex_num> P, int num_threads, bool inverse){
    // compute the FFT of the vector P in parallel with num_threads threads
    int n = P.size();
    
    // find p closest to sqrt(n) 
    int p = sqrt(n);
    while(n % p != 0){
        p--;
    }
    int q = n / p;

    // see P as a matrix of dimensions p x q

    // do FFT on each column of the matrix P
    std::vector<std::thread> workers(num_threads - 1);
    for (size_t i = 0; i < num_threads - 1; ++i) {
        workers[i] = std::thread(&FFTThread, std::ref(P), i, num_threads, p, q, inverse);
    }
    FFTThread(P, num_threads - 1, num_threads, p, q, inverse);
    for (size_t i = 0; i < num_threads - 1; ++i) {
        workers[i].join();
    }
    
    vector<complex_num> Q(n);

    if (! inverse ){
        // transpose the matrix P
        std::vector<std::thread> workers1(num_threads - 1);
        for (size_t i = 0; i < num_threads - 1; ++i) {
            workers1[i] = std::thread(&TransposeThread, std::ref(P), std::ref(Q), i, num_threads, p, q);
        }
        TransposeThread(P, Q, num_threads -1, num_threads, p, q);
        for (size_t i = 0; i < num_threads - 1; ++i) {
            workers1[i].join();
        }

        // do twiddle factor multiplication on transposed matrix Q (Q is now a q x p matrix)
        std::vector<std::thread> workers2(num_threads - 1);
        for (size_t i = 0; i < num_threads - 1; ++i) {
            workers2[i] = std::thread(&TwiddleThread, std::ref(Q), i, num_threads, q, p, false);
        }
        TwiddleThread(Q, num_threads-1, num_threads, q, p, false);
        for (size_t i = 0; i < num_threads - 1; ++i) {
            workers2[i].join();
        }
    } else {
        // do twiddle factor multiplication on matrix P
        std::vector<std::thread> workers2(num_threads - 1);
        for (size_t i = 0; i < num_threads - 1; ++i) {
            workers2[i] = std::thread(&TwiddleThread, std::ref(P), i, num_threads, p, q, true);
        }
        TwiddleThread(P, num_threads-1, num_threads, p, q, true);
        for (size_t i = 0; i < num_threads - 1; ++i) {
            workers2[i].join();
        }

        // transpose the matrix P
        std::vector<std::thread> workers1(num_threads - 1);
        for (size_t i = 0; i < num_threads - 1; ++i) {
            workers1[i] = std::thread(&TransposeThread, std::ref(P), std::ref(Q), i, num_threads, p, q);
        }
        TransposeThread(P, Q, num_threads -1, num_threads, p, q);
        for (size_t i = 0; i < num_threads - 1; ++i) {
            workers1[i].join();
        }
    }


    // do FFT on each column of the matrix Q

    // parallel version
    for (size_t i = 0; i < num_threads - 1; ++i) {
        workers[i] = std::thread(&FFTThread, std::ref(Q), i, num_threads, q, p, inverse);
    }
    FFTThread(Q, num_threads - 1, num_threads, q, p, inverse);
    for (size_t i = 0; i < num_threads - 1; ++i) {
        workers[i].join();
    }

    return Q;
}

