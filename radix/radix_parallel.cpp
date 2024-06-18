#include "radix_parallel.h"

using namespace std;

void get_even_odd_thread(vector<complex_num> &U, vector<complex_num> &V, vector<complex_num> P, int begin, int end){
    for (size_t j = begin; j < end; j++){
        U[j] = P[2*j];
        V[j] = P[2*j + 1];
    }
}

void radix_thread(vector<complex_num> &P_star, vector<complex_num> P, int num_threads, bool inverse = false){
    int n = P.size();
    if (n == 1){
        P_star[0] = P[0];
        return;
    }

    vector<complex_num> U(n/2);
    vector<complex_num> V(n/2);

    for (size_t j = 0; j < n/2; j++){
        U[j] = P[2*j];
        V[j] = P[2*j + 1];
    }

    vector<complex_num> U_star(n/2);
    vector<complex_num> V_star(n/2);

    if (num_threads == 1){
        radix_thread(U_star, U, 1, inverse);
        radix_thread(V_star, V, 1, inverse);
    } else {
        thread t1(radix_thread, ref(U_star), U, num_threads/2, inverse);
        radix_thread(V_star, V,  num_threads - num_threads/2, inverse);
        t1.join();
    }
    int sign = inverse ? 1 : -1;
    complex_num omega_n(cos(2*M_PI/n), sign*sin(2*M_PI/n));
    complex_num omega(1, 0);

    for (size_t j=0; j < n/2; j++){
        double factor = inverse ? 0.5 : 1;
        P_star[j] = factor * (U_star[j] + omega*V_star[j]);
        P_star[j + n/2] = factor * (U_star[j] - omega*V_star[j]);
        omega *= omega_n;
    }
}

vector<complex_num> Radix2FFT_parallel(vector<complex_num> P, int num_threads, bool inverse){
    size_t n_orig = P.size();
    if (n_orig == 1){
        return vector<complex_num>{P[0]};
    }
    // pad to power of 2
    int p = 1;
    while (p < P.size()){
        p *= 2;
    }
    
    P.insert(P.end(), p - n_orig, complex_num(0.0, 0.0));
    int n = P.size();
    
    vector<complex_num> P_star(n);
    radix_thread(P_star, P, num_threads, inverse);

    // remove Padding
    while (P_star.size() > n_orig)
        P_star.pop_back();
    return P_star;
}