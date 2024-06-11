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
    complex_num X3[N1][N2][N3];
    for (int k1 = 0; k1<N1; k1++){
        for (int k2 = 0; k2<N2; k2++){
            for (int k3 = 0; k3<N3; k3++){
                complex_num sm = complex_num(0,0);
                for (int j3 = 0; j3 < N3; j3++){
                    sm += X2[j3][k1][k2]*polar(1.0, 2.0 * M_PI * (j3*k1+k2*N1+k3*N2) / (N1*N2*N3))*polar(1.0, 2.0 * M_PI * (j3*k3) / N3);
                }
                X3[k1][k2][k3] = sm;
            }
        }
    }

    vector<complex_num> output(x.size());
    // Flatten - pretty sure this is wrong
    for (int i1 = 0; i1<N1; i1++){
        for (int i2 = 0; i2<N2; i2++){
            for (int i3 = 0; i3<N3; i3++){
                output[i3 + i2*N3 + i1*N3*N2] = X3[i1][i2][i3];
            }
        }
    }
    return output;
}

void printVector(const vector<complex_num>& vec) {
    for (const auto& c : vec) {
        cout << c << " ";
    }
    cout << endl;
}

int main() {
    // Test case 1:
    vector<complex_num> input1 = {1.0, 1.0, 1.0, 1.0};
    cout << endl << "Input 1: " << endl;
    vector<complex_num> output1 = DFT(input1);
    printVector(output1);
    cout << endl;
    vector<complex_num> output1ct = FFT(input1, 2, 2, 1);
    printVector(output1ct);

    // Test case 2:
    vector<complex_num> input2 = {1.0, 2.0, 3.0, 4.0};
    cout << endl << "Input 2: " << endl;
    vector<complex_num> output2 = DFT(input2);
    printVector(output2);
    cout << endl;
    vector<complex_num> output2ct = FFT(input2, 2, 2, 1);
    printVector(output2ct);

    // Test case 3:
    vector<complex_num> input3 = {0.0, 1.0, 0.0, -1.0};
    cout << endl << "Input 3: " << endl;
    vector<complex_num> output3 = DFT(input3);
    printVector(output3);
    cout << endl;
    vector<complex_num> output3ct = FFT(input3, 2, 2, 1);
    printVector(output3ct);

    // Test case 4:
    vector<complex_num> input4 = {1.0, 2.0, 3.0, 4.0};
    cout << endl << "Input 4: " << endl;
    vector<complex_num> output4 = DFT(input4);
    printVector(output4);
    cout << endl;
    vector<complex_num> output4ct = FFT(input4, 2, 2, 1);
    printVector(output4ct);

    return 0;
}