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
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef long double num;
typedef std::vector<long double>::const_iterator NumIter;

typedef std::complex<num> complex_num;


using namespace std;

//-----------------------------------------------------------------------------

std::tuple<int, int> Decompose(int n){
    int p, q;
    for(int i = 2; i*i <= n; i++){
        if (n % i == 0){
            p = i;
            q = n/i;
            return std::make_tuple(p, q);
        }
    }
    return std::make_tuple(1, n);
}

// for prime numbers, this primitive FFT, but can be changed to Rader's algorithm
vector<complex_num> StandardFFT(vector<complex_num> &P){
    size_t n = P.size();
    std::vector<std::complex<long double>> output(n, 0);
    std::complex<long double> omega;
    for(size_t k = 0; k < n; k++){
        omega = {1.0, 0.0};
        std::complex<long double> sum(0, 0);
        long double angle = 2 * M_PI * k / n;
        std::complex<long double> omega_l(std::cos(angle), -std::sin(angle));
        for(size_t l = 0; l < n; l++){
            sum += P[l] * omega;
            omega *= omega_l;
        }
        output[k] = sum;
    }
    return output;
}

std::vector<complex_num> TransposeFFT(std::vector<complex_num> &P, int p, int q){
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


void TwistFFT(std::vector<complex_num> &P, int p, int q){
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

// std::vector<std::vector<complex_num>> PackFFT(vector<complex_num> &P, int p, int q){
//     std::vector<std::vector<complex_num>> A(p, std::vector<complex_num>(q));
//     for(int j = 0; j < p; j++){
//         for(int k = 0; k < q; k++){
//             A[j][k] = P[j * p + k];
//         }
//     }
//     return A;
// }

// vector<complex_num> UnPackFFT(std::vector<std::vector<complex_num>> &P, int p, int q){

//     std::vector<std::complex<long double>> output(p * q, 0);
//     for(int j = 0; j < p; j++){
//         for(int k = 0; k < q; k++){
//             output[j * p + k] = P[j][k];
//         }
//     }
//     return output;
// }

void PRT1(std::vector<std::complex<long double>> P, string a = ""){
    cout<< a <<endl;
    cout<<endl;
    int p;
    p = P.size();
    for(int i=0; i<p; i++){
        cout<<P[i]<<"  ";
    }
    cout<<endl;
    cout<<endl;
}

void PRT2(std::vector<std::vector<complex_num>> P, string a = ""){
    cout<< a <<endl;
    cout<<endl;
    int p,q;
    p = P.size();
    q = P[0].size();
    for(int i=0; i<p; i++){
        for(int j=0; j<q; j++){
        cout<<P[i][j]<<" ";
        }
        cout<<endl;
        cout<<endl;
    }
    cout<<endl;
}

vector<complex_num> GeneralFFT(vector<complex_num> &P, bool f = true){
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
    A=TransposeFFT(P, p, q);
    
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
        result = GeneralFFT(temp, false);
        B.insert(B.end(), result.begin(), result.end());
        first += p;
        last += p;
    }
    if(f == true){
        PRT1(B, " <- 2th step");
    }


    A.clear();
    
    //3 step: Twiddle factors
    TwistFFT(B, q, p);
    
    if(f == true){
        PRT1(B, " <- 3th step");
    }

    //4 step: General transpose
    A = TransposeFFT(B, q, p);

    B.clear();

    if(f == true){
        PRT1(A, " <- 4th step");
    }
    first = 0;
    last = q;

    for(int k = 0; k < p; k++){
        temp.assign(A.begin() + first, A.begin() + last);
        result = GeneralFFT(temp, false);
        B.insert(B.end(), result.begin(), result.end());
        first += q;
        last += q;
    }

    A.clear();

    if(f == true){
        PRT1(B, " <- 5th step");
    }
    
    A = TransposeFFT(B, p, q);

    if(f == true){
        PRT1(A, " <- 6th step");
    }
    
    return A;
}

vector<complex_num> Read_CSV(string file_path){
    std::ifstream file(file_path);
    std::string line;
    std::vector<complex_num> column_data;
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file");
    }
    long double a = 0.0;
    complex_num omega;
    while (std::getline(file, line)) {
        // Convert each line (string) to double and store in the vector
        omega = {std::stold(line), a};
        column_data.push_back(omega);
    }
    file.close();
    return column_data;
}

int main(){
    // This just generates a vector of length 15 of (1, 0)
    // vector<complex_num> P(15, complex_num(1, 0));
    
    vector<complex_num> P = Read_CSV("Weather_data.csv");

    std::vector<std::complex<long double>> vec;

    // Populate the vector with complex numbers
    for (int i = 0; i < 8; ++i) {
        vec.push_back(std::complex<double>(i + 1, 0.0)); // Real part is i + 1, imaginary part is 0.0
    }
    P = vec;
    
    P =  std::vector<complex_num>(40320, complex_num(1.0, 0.0));
    int t = P.size();
    // for(int i = 0; i < t; i++){
    //     cout<<P[i].real()<<","<<" ";
    // }
    cout<<endl<<endl;
    
    vector<complex_num> P_star = GeneralFFT(P, false);

    

    // for(int i = 0; i < t; i++){
    //     cout<<P_star[i]<<" ";
    // }
    cout<<endl<<t;
    return 0;
}