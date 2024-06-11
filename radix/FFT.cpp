#include "FFT.h"


using namespace std;

//-----------------------------------------------------------------------------

vector<complex_num> Radix2FFT(vector<complex_num> P){
    size_t n_orig = P.size();
    if (n_orig == 1){
        return vector<complex_num>{P[0]};
    }
    // pad to power of 2
    while(log2(P.size()) > (int) log2(P.size())){
        P.push_back(complex_num(0.0,0.0));
    }

    int n = P.size();
    vector<complex_num> U(n/2);
    vector<complex_num> V(n/2);

    for (size_t j = 0; j < n/2; j++){
        U[j] = P[2*j];
        V[j] = P[2*j + 1];
    }
    
    vector<complex_num> U_star = Radix2FFT(U);
    vector<complex_num> V_star = Radix2FFT(V);

    complex_num omega_n(cos(2*M_PI/n), -sin(2*M_PI/n));
    complex_num omega(1, 0);

    vector<complex_num> P_star(n);

    for (size_t j=0; j < n/2; j++){
        P_star[j] = U_star[j] + omega*V_star[j];
        P_star[j + n/2] = U_star[j] - omega*V_star[j];
        omega *= omega_n;
    }

    // remove Padding
    while (P_star.size() > n_orig)
        P_star.pop_back();
    return P_star;
}

vector<complex_num> InverseRadix2FFT(vector<complex_num> P_star){
    size_t n_orig = P_star.size();
    if (n_orig == 1){
        return vector<complex_num>{P_star[0]};
    }
    // pad to power of 2
    while(log2(P_star.size()) > (int) log2(P_star.size()))
        P_star.push_back(complex_num(0.0,0.0));

    int n = P_star.size();
    vector<complex_num> U_star(n/2);
    vector<complex_num> V_star(n/2);

    for (size_t j = 0; j < n/2; j++){
        U_star[j] = P_star[2*j];
        V_star[j] = P_star[2*j + 1];
    }
    
    vector<complex_num> U = InverseRadix2FFT(U_star);
    vector<complex_num> V = InverseRadix2FFT(V_star);

    complex_num omega_n(cos(2*M_PI/n), +sin(2*M_PI/n));
    complex_num omega(1, 0);

    vector<complex_num> P(n);

    for (size_t j=0; j < n/2; j++){
        P[j] = (U[j] + omega*V[j])/complex(num(2), num(0));
        P[j + n/2] = (U[j] - omega*V[j])/complex(num(2), num(0));
        omega *= omega_n;
    }

    // remove Padding
    while (P.size() > n_orig)
        P.pop_back();
    return P;
}
