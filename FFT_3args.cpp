#include <vector>
#include <complex>
#include <cmath>
#include <iostream>
#include <thread>
#include <iostream>
#include <fstream>
#include <sstream>
#include <climits>
#include <thread>
#include <numeric>
#include <iterator>
#include <vector>
#include <string>
#include <complex>
#include <cmath>

using namespace std;
typedef complex<double> complex_num;

vector<complex_num> FFT(vector<complex_num>& input);

/*
Basic DFT
*/
vector<complex_num> DFT(vector<complex_num>&x) {
    size_t N = x.size();
    vector<complex_num> X = vector<complex_num>(N);
    int k, n;
    for(k = 0; k < N; k++) {
        for(n = 0; n < N; n++) {
            double angle = -2 * M_PI * n * k / N;
            X[k] += x[n] * polar(1.0, angle);
        }
    }
    
    return X;
}

vector<complex_num> FFT(vector<complex_num>&x, int N1, int N2, int N3) {
    // Make X0
    complex_num X0[N1][N2][N3];
    for (int i1 = 0; i1<N1; i1++){
        for (int i2 = 0; i2<N2; i2++){
            for (int i3 = 0; i3<N3; i3++){
                X0[i1][i2][i3] = x[i3 + i2*N3 + i1*N3*N2];
            }
        }
    }
    // Step 1
    complex_num X1[N1][N2][N3];
    for (int j2 = 0; j2<N1; j2++){
        for (int j3 = 0; j3<N2; j3++){
            for (int k1 = 0; k1<N3; k1++){
                complex_num sm = complex_num(0,0);
                for (int j1 = 0; j1 < N1; j1++){
                    sm += X0[j1][j2][j3]*polar(1.0, 2.0 * M_PI * j1 * k1 / N1);
                }
                X1[j2][j3][k1] = sm;
            }
        }
    }
    // Step 2
    complex_num X2[N1][N2][N3];
    for (int j3 = 0; j3<N1; j3++){
        for (int k1 = 0; k1<N2; k1++){
            for (int k2 = 0; k2<N3; k2++){
                complex_num sm = complex_num(0,0);
                for (int j2 = 0; j2 < N2; j2++){
                    sm += X1[j2][j3][k1]*polar(1.0, 2.0 * M_PI * j2 * k1 / (N1*N2))*polar(1.0, 2.0 * M_PI * j2 * k2 / N2);
                }
                X2[j3][k1][k2] = sm;
            }
        }
    }
    // Step 3
    complex_num X2[N1][N2][N3];
    for (int k1 = 0; k1<N1; k1++){
        for (int k2 = 0; k2<N2; k2++){
            for (int k3 = 0; k3<N3; k3++){
                complex_num sm = complex_num(0,0);
                for (int j3 = 0; j3 < N3; j3++){
                    sm += X1[j3][k1][k2]*polar(1.0, 2.0 * M_PI * (j3*k1+k2*N1+k3*N2) / (N1*N2*N3))*polar(1.0, 2.0 * M_PI * (j3*k3) / N3);
                }
                X2[k1][k2][k3] = sm;
            }
        }
    }

    vector<complex_num> output;
    
    return output;
}