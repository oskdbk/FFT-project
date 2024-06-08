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
    for(int k = 0; k < n; k++){
        omega = {1.0, 0.0};
        std::complex<long double> sum(0, 0);
        long double angle = 2 * M_PI * k / n;
        std::complex<long double> omega_l(std::cos(angle), -std::sin(angle));
        for(int l = 0; l < n; l++){
            sum += P[l] * omega;
            omega *= omega_l;
        }
        output[k] = sum;
    }
    return output;
}

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
    std::complex<long double> omega;
    long double angle;
    std::vector<std::vector<complex_num>> Q(p, std::vector<complex_num>(q));
    for (size_t i = 0; i < p; ++i) {
        
        for (size_t j = 0; j < q; ++j) {
            angle = 2 * M_PI * i * j / n;
            std::complex<long double> omega(std::cos(angle), -std::sin(angle));
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

    std::vector<std::complex<long double>> output(p * q, 0);
    for(int j = 0; j < p; j++){
        for(int k = 0; k < q; k++){
            output[j * q + k] = P[j][k];
        }
    }
    return output;
}

void PRT1(std::vector<std::complex<long double>> P, string a = ""){
    cout<< a << ":" << endl;
    int p;
    p = P.size();
    for(int i=0; i<p; i++){
        cout<<P[i]<<"  ";
    }
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

vector<complex_num> GeneralFFT(vector<complex_num> &P, bool f = false){
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

    vector<complex_num> P_star = GeneralFFT(P, false, true);

    
    size_t t = P_star.size();
    for(int i = 0; i < t; i++){
        cout<<P_star[i]<<" ";
    }
    cout<<endl<<t;
    return 0;
}
