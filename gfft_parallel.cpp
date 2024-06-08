
#include <vector>
#include <string>
#include <cmath>
#include <numeric>
#include <complex> 
#include "General_FFT_algo.cpp"
typedef long double num;
typedef std::complex<num> complex_num;
using namespace std;


// ------------------parallel versions----------------- //

void FFTThread(vector<complex_num> &P, int thread_id, int num_threads, int p, int q){
    // input: matrix P of dimensions p x q
    // in all columns of P that are = thread_id mod num_threads,
    // replace the column by its FFT
    int n = P.size();
    for(int col = thread_id; col < q; col += num_threads){
        vector<complex_num> A(p);
        for(int i = 0; i < p; i++){
            A[i] = P[i * q + col];
        }
        A = GeneralFFT(A);
        for(int i = 0; i < p; i++){
            P[i * q + col] = A[i];
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
        
        int row = i / q; // row index of the element
        int col = i % q; // column index of the element
        // if(i<20 && thread_id == 1){
        //     cout<<i<<" ("<<row<<","<<col<<"): " << P[i];
        // }
        angle = 2 * M_PI * row * col / n;
        omega = std::complex<long double>(std::cos(angle), -std::sin(angle));
        complex<long double> temp = P[i] * omega;
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
    std::complex<long double> omega;
    long double angle;
    for(int i = 0; i < n; i++){ // i is the raw index of the element
        int row = i / q; // row index of the element
        int col = i % q; // column index of the element
        angle = 2 * M_PI * row * col / n;
        omega = std::complex<long double>(std::cos(angle), -std::sin(angle));
        P1[i] = P[i] * omega;
        // if(i<20)
        //     cout << P[i] << " * " << omega << " = " << P1[i] << endl;
    }
    return P1;
}

// --------------------------------------------------------------------- //




vector<complex_num> FFT_Parallel(vector<complex_num> P, int num_threads = 1){
    // compute the FFT of the vector P in parallel with num_threads threads
    int n = P.size();
    while(n % num_threads != 0){
        num_threads--; // make sure n is divisible by num_threads
    }
    int p = n / num_threads; // number of rows we divide P into
    int q = num_threads; // number of columns = number of threads

    // TODO: different decomposition of P into p x q also possible
    // doesn't have to be num_threads columns

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

    // do twiddle factor multiplication on transposed matrix Q (Q is now a q x p matrix)
    std::vector<std::thread> workers2(num_threads - 1);
    for (size_t i = 0; i < num_threads - 1; ++i) {
        workers2[i] = std::thread(&TwiddleThread, std::ref(Q), i, num_threads, q, p);
    }
    TwiddleThread(Q, num_threads-1, num_threads, q, p);
    for (size_t i = 0; i < num_threads - 1; ++i) {
        workers2[i].join();
    }

    // do FFT on each column of the matrix Q

    // parallel version
    for (size_t i = 0; i < num_threads - 1; ++i) {
        workers[i] = std::thread(&FFTThread, std::ref(Q), i, num_threads, q, p);
    }
    FFTThread(Q, num_threads - 1, num_threads, q, p);
    for (size_t i = 0; i < num_threads - 1; ++i) {
        workers[i].join();
    }

    return Q;
}

int main(){
    // This just generates a vector of length 15 of (1, 0)
    // vector<complex_num> P(15, complex_num(1, 0));

    vector<complex_num> P = Read_CSV("Weather_data.csv");
    int n = P.size();
    // copy P to Q
    vector<complex_num> Q(P.size());
    for(int i = 0; i < P.size(); i++){
        Q[i] = P[i];
    }
    // cout<<P.size()<<endl;
    

    // measure time of both versions

    // sequential version
    clock_t start, end;
    double avg_time = 0;
    for (int i = 0; i < 30; i++)
    {
        start = clock();
        vector<complex_num> P_star = GeneralFFT(P, false);
        end = clock();
        avg_time += double(end - start);
    }
    avg_time /= 30;

    cout << "Sequential version took " <<  avg_time << " ticks" << endl;

    
    
    // try with different number of threads
    for(int i = 1; i <= 30; i++){
        if (n % i != 0){
            continue;
        }
        double avg_time = 0;
        for(int j = 0; j < 30; j++){
            start = clock();
            vector<complex_num> P_star_parallel = FFT_Parallel(Q, i);
            end = clock();
            avg_time += double(end - start);
        }
        avg_time /= 30;
        cout << "Parallel version with " << i << " threads took " << avg_time << " ticks" << endl;
    }


    vector<complex_num> P_star = GeneralFFT(P, false);
    vector<complex_num> P_star_parallel = FFT_Parallel(Q, 20);
    size_t t = P_star.size();
    for(int i = 0; i < t; i++){
        if(norm(P_star[i] - P_star_parallel[i]) > 1e-6){
            cout << "Parallel differs at index " << i << endl;
            cout << P_star[i] << " != " << P_star_parallel[i] << endl;
            break;
        } 
    }
    cout << "Parallel is correct" << endl;

    // Print results
    cout << "Result:" << endl;
    for (int i = 0; i < 9; i++)
    {
        cout << P_star[i] << " ";
    }
    cout << endl;
    return 0;
}