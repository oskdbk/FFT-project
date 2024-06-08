
#include <vector>
#include <string>
#include <cmath>
#include <numeric>
#include <complex> 
#include "General_FFT_algo.cpp"
typedef long double num;
typedef std::complex<num> complex_num;
using namespace std;
void FFTThread(vector<complex_num> &P, int thread_id, int num_threads, int p, int q){
    // input: matrix P of dimensions p x q
    // in all columns of P that are = thread_id mod num_threads,
    // replace the column by its FFT
    int n = P.size();
    for(int col = thread_id; col < q; col += num_threads){
        vector<complex_num> A(p);
        for(int row = 0; row < p; row++){
            A[row] = P[row * q + col];
        }
        A = GeneralFFT(A);
        for(int row = 0; row < p; row++){
            P[row * q + col] = A[row];
        }
    }
}
void TwiddleThread(vector<complex_num> &P, int thread_id, int num_threads, int p, int q){
    // input: matrix P of dimensions p x q
    // for all elements of P that are = thread_id mod num_threads,
    // multiply the element by the twiddle factor
    int n = P.size();
    std::complex<long double> omega;
    long double angle;
    for(int i = thread_id; i < n; i += num_threads){ // i is the raw index of the element
        int row = i % p; // row index of the element
        int col = i / p; // column index of the element
        angle = 2 * M_PI * row * col / n;
        omega = std::complex<long double>(std::cos(angle), -std::sin(angle));
        P[i] = P[i] * omega;
    }
}

void TransposeThread(vector<complex_num> &P, vector<complex_num> &Q, int thread_id, int num_threads, int p, int q){
    // transpose the pxq matrix P into the qxp matrix Q
    // do the indexes from begin to end of Q
    int n = P.size();
    for (size_t i = thread_id; i < n; i+= num_threads)
    {
        int row = i % q;
        int col = i / q;
        Q[i] = P[row * p + col];
    }   
}
vector<complex_num> FFT_Parallel(vector<complex_num> &P, int num_threads = 1){
    int n = P.size();
    while(n % num_threads != 0){
        num_threads--;
    }
    int p = n / num_threads; // number of rows we divide P into
    int q = num_threads; // number of columns = number of threads

    // see P as a matrix of dimensions p x q
    // do FFT on each column of the matrix P
    std::vector<std::thread> workers(num_threads - 1);
    for (size_t i = 0; i < num_threads - 1; ++i) {
        workers[i] = std::thread(&FFTThread, std::ref(P), i, num_threads, p, q);
    }
    FFTThread(P, num_threads - 1, num_threads, p, q);
    for (size_t i = 0; i < num_threads - 1; ++i) {
        workers[i].join();
    }
    
    // transpose the matrix P
    vector<complex_num> Q(n);
    std::vector<std::thread> workers1(num_threads - 1);
    for (size_t i = 0; i < num_threads - 1; ++i) {
        workers1[i] = std::thread(&TransposeThread, std::ref(P), std::ref(Q), i, num_threads, p, q);
    }
    TransposeThread(P, Q, num_threads -1, num_threads, p, q);
    for (size_t i = 0; i < num_threads - 1; ++i) {
        workers1[i].join();
    }

    // Q has dimensions q x p

    // do twiddle factor multiplication on transposed matrix Q
    std::vector<std::thread> workers2(num_threads - 1);
    for (size_t i = 0; i < num_threads - 1; ++i) {
        workers2[i] = std::thread(&TwiddleThread, std::ref(Q), i, num_threads, q, p);
    }
    TwiddleThread(Q, num_threads-1, num_threads, q, p);
    for (size_t i = 0; i < num_threads - 1; ++i) {
        workers2[i].join();
    }


    // do FFT on each column of the matrix Q
    std::vector<std::thread> workers3(num_threads - 1);
    for (size_t i = 0; i < num_threads - 1; ++i) {
        workers3[i] = std::thread(&FFTThread, std::ref(Q), i, num_threads, q, p);
    }
    FFTThread(P, num_threads - 1, num_threads, q, p);
    for (size_t i = 0; i < num_threads - 1; ++i) {
        workers3[i].join();
    }

    return Q;
}