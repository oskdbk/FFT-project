#include "utils.h"

using namespace std;

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

void PRT1(std::vector<std::complex<long double>> P, string a){
    cout<< a << ":" << endl;
    int p;
    p = P.size();
    for(int i=0; i<p; i++){
        cout<<P[i]<<"  ";
    }
    cout<<endl;
}

void PRT2(std::vector<std::vector<complex_num>> P, string a){
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

/*
Basic DFT
*/
vector<complex_num> DFT(vector<complex_num>x) {
    size_t N = x.size();
    vector<complex_num> X = vector<complex_num>(N);
    int k, n;
    for(k = 0; k < N; k++) {
        for(n = 0; n < N; n++) {
            double angle = 2 * M_PI * n * k / N;
            
            X[k] += x[n] * complex_num(std::cos(angle), -std::sin(angle));
        }
    }
    
    return X;
}

bool is_same_vector(vector<complex_num> A, vector<complex_num> B, string a, string b, bool verbose){
    size_t t = A.size();
    if (t != B.size()){
        if (verbose)
            cout << a << " and " << b << " have different sizes" << endl;
        return false;
    }
    for(int i = 0; i < t; i++){
        if(norm(A[i] - B[i]) > 1e-6){
            if (verbose)
                cout << a << " and " << b << " are different at index" << i;
                cout << "A[" << i << "] = " << A[i] << ", B[" << i << "] = " << B[i] << endl;
            return false;
        } 
    }
    if (verbose)
        cout << a << " and " << b << " are the same" << endl;
    return true;
}
