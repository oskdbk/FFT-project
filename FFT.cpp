#include <iostream>
#include <climits>
#include <thread>
#include <numeric>
#include <iterator>
#include <optional>
#include <vector>
#include <complex> // to use complex numbers
#include <cmath>

typedef double num;
typedef std::vector<double>::const_iterator NumIter;

typedef std::complex<num> complex_num;


using namespace std;

//-----------------------------------------------------------------------------

vector<complex_num> Radix2FFT(vector<complex_num> &P){
    size_t n = P.size();
    if (n == 1){
        return vector<complex_num>{P[0]};
    }

    vector<complex_num> U(n/2);
    vector<complex_num> V(n/2);

    for (size_t j = 0; j < n/2; j++){
        U[j] = P[2*j];
        V[j] = P[2*j + 1];
    }
    
    vector<complex_num> U_star = Radix2FFT(U);
    vector<complex_num> V_star = Radix2FFT(V);

    complex_num omega_n(cos(2*M_PI/n), sin(2*M_PI/n));
    complex_num omega(1, 0);

    vector<complex_num> P_star(n);

    for (size_t j=0; j < n/2; j++){
        P_star[j] = U_star[j] + omega*V_star[j];
        P_star[j + n/2] = U_star[j] - omega*V_star[j];
        omega *= omega_n;
    }

    return P_star;
}

vector<complex_num> InverseRadix2FFT(vector<complex_num> &P_star){
    size_t n = P_star.size();
    if (n == 1){
        return vector<complex_num>{P_star[0]};
    }

    vector<complex_num> U_star(n/2);
    vector<complex_num> V_star(n/2);

    for (size_t j = 0; j < n/2; j++){
        U_star[j] = P_star[2*j];
        V_star[j] = P_star[2*j + 1];
    }
    
    vector<complex_num> U = InverseRadix2FFT(U_star);
    vector<complex_num> V = InverseRadix2FFT(V_star);

    complex_num omega_n(cos(2*M_PI/n), -sin(2*M_PI/n));
    complex_num omega(1, 0);

    vector<complex_num> P(n);

    for (size_t j=0; j < n/2; j++){
        P[j] = (U[j] + omega*V[j])/complex(num(2), num(0));
        P[j + n/2] = (U[j] - omega*V[j])/complex(num(2), num(0));
        omega *= omega_n;
    }

    return P;
}


int main(){
    vector<complex_num> P = {complex_num(2, 0), complex_num(1, 0), complex_num(3, 0), complex_num(7, 0)};
    vector<complex_num> P_star = Radix2FFT(P);
    cout << "P star:" << endl;
    for (int i = 0; i < P_star.size(); i++){
        cout << P_star[i] << endl;
    }
    vector<complex_num> P2 = InverseRadix2FFT(P_star);
    cout << "P2:" << endl;
    for (int i = 0; i < P2.size(); i++){
        cout << P2[i] << endl;
    }
    return 0;
}