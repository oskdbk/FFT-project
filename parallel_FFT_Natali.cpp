#include <vector>
#include <complex>
#include <cmath>
#include <iostream>
#include <thread>

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

void computeDFTColumn(vector<vector<complex_num>>& matrix, int index) {
    matrix[index] = FFT(matrix[index]);
}

void computeTwiddle(vector<vector<complex_num>>& columns, vector<vector<complex_num>>& rows, int N1, int N2, int index) {
    int N = N1 * N2;
    for (int k2 = 0; k2 < N2; k2++) {
        rows[k2][index] = polar(1.0, -2.0 * M_PI * index * k2 / N) * columns[index][k2];
    }
}

void computeDFTRow(vector<vector<complex_num>>& matrix, int index) {
    matrix[index] = FFT(matrix[index]);
}

void flattenOutput(vector<vector<complex_num>>& rows, vector<complex_num>& output, int N1, int N2, int index) {
    for (int k2 = 0; k2 < N2; k2++) {
        output[N2 * index + k2] = rows[k2][index];
    }
}


/*
Parallel CooleyTukey for N = N1 x N2
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
    vector<thread> threads;
    for (int k1 = 0; k1 < N1; k1++) {
        threads.emplace_back(computeDFTColumn, ref(columns), k1);
    }
    for (auto& t : threads) {
        t.join();
    }
    threads.clear();

    // Loading rows by multiplication of columns by twiddle factors
    for (int k1 = 0; k1 < N1; k1++) {
        threads.emplace_back(computeTwiddle, ref(columns), ref(rows), N1, N2, k1);
    }
    for (auto& t : threads) {
        t.join();
    }
    threads.clear();

    // DFT per row
    for (int k2 = 0; k2 < N2; k2++) {
        threads.emplace_back(computeDFTRow, ref(rows), k2);
    }
    for (auto& t : threads) {
        t.join();
    }
    threads.clear();

    // Finally flatten
    vector<complex_num> output(N);
    for (int k1 = 0; k1 < N1; k1++) {
        threads.emplace_back(flattenOutput, ref(rows), ref(output), N1, N2, k1);
    }
    for (auto& t : threads) {
        t.join();
    }

    return output;
}

vector<complex_num> FFT(vector<complex_num>& input) {
    // find good N1 and N2 to divide by
    int N = input.size();
    if (N < 10){
        return DFT(input);
    }
    double n = sqrt(N);
    for (int i = floor(sqrt(N)); i > 0; i--){
        if (N % i == 0){
            return CooleyTukey(input, i, N%i);
        }
    }
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