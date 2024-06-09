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

typedef double num;
typedef std::complex<num> complex_num;

using namespace std;

/*
Naive DFT function with n^2 complexity
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

/*
Cooley-Tukey algorithm
*/
vector<complex_num> CooleyTukey(vector<complex_num>& input, int N1, int N2) {
    size_t N = input.size();

    // Create empty column and row matrices
    vector<vector<complex_num>> columns(N1, vector<complex_num>(N2));
    vector<vector<complex_num>> rows(N2, vector<complex_num>(N1));

    // Load column matrix
    for (int k1 = 0; k1 < N1; k1++) {
        for (int k2 = 0; k2 < N2; k2++) {
            columns[k1][k2] = input[N1 * k2 + k1];
        }
    }

    // DFT per column
    for (int k1 = 0; k1 < N1; k1++) {
        columns[k1] = DFT(columns[k1]);
    }

    // Loading rows by multiplication of columns by twiddle factors
    for (int k1 = 0; k1 < N1; k1++) {
        for (int k2 = 0; k2 < N2; k2++) {
            rows[k2][k1] = polar(1.0, 2.0 * M_PI * k1 * k2 / N) * columns[k1][k2];
        }
    }

    // DFT per row
    for (int k2 = 0; k2 < N2; k2++) {
        rows[k2] = DFT(rows[k2]);
    }

    // Finally flatten
    vector<complex_num> output(N);
    for (int k1 = 0; k1 < N1; k1++) {
        for (int k2 = 0; k2 < N2; k2++) {
            output[N2 * k1 + k2] = rows[k2][k1];
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
    // Test case 1: Simple case
    vector<complex_num> input1 = {1.0, 1.0, 1.0, 1.0};
    cout << endl << "Input 1: " << endl;
    vector<complex_num> output1 = DFT(input1);
    printVector(output1);
    cout << endl;
    vector<complex_num> output1ct = CooleyTukey(input1, 2, 2);
    printVector(output1ct);

    // Test case 2: Another simple case
    vector<complex_num> input2 = {1.0, 2.0, 3.0, 4.0};
    cout << endl << "Input 2: " << endl;
    vector<complex_num> output2 = DFT(input2);
    printVector(output2);
    cout << endl;
    vector<complex_num> output2ct = CooleyTukey(input2, 2, 2);
    printVector(output2ct);

    // Test case 3: Sinusoidal input
    vector<complex_num> input3 = {0.0, 1.0, 0.0, -1.0};
    cout << endl << "Input 3: " << endl;
    vector<complex_num> output3 = DFT(input3);
    printVector(output3);
    cout << endl;
    vector<complex_num> output3ct = CooleyTukey(input3, 2, 2);
    printVector(output3ct);

    // Test case 4: FFT with Cooley-Tukey algorithm
    vector<complex_num> input4 = {1.0, 2.0, 3.0, 4.0};
    cout << endl << "Input 4: " << endl;
    vector<complex_num> output4 = DFT(input4);
    printVector(output4);
    cout << endl;
    vector<complex_num> output4ct = CooleyTukey(input4, 2, 2);
    printVector(output4ct);

    return 0;
}
