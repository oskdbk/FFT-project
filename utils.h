#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <string>
#include <complex>
#include <cmath>
#include <iostream>
#include <fstream>
#include <tuple>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;

typedef double num;
typedef std::complex<num> complex_num;

std::vector<complex_num> Read_CSV(std::string file_path);
void PRT1(std::vector<std::complex<double>> P, std::string a);
void PRT2(std::vector<std::vector<complex_num>> P, std::string a);
std::tuple<int, int> Decompose(int n);
std::vector<complex_num> StandardFFT(std::vector<complex_num> &P, bool inv = false);
std::vector<complex_num> DFT(std::vector<complex_num>x, bool inverse = false);
bool is_same_vector(vector<complex_num> A, vector<complex_num> B, string a = "A", string b = "B", bool verbose = false);

// Template function moved here
template<typename F, typename... Args>
double measure_time(F func, int iters, Args&&... args) {
    clock_t start, end;
    double avg_time = 0;
    for (int i = 0; i < iters; i++) {
        start = clock();
        func(std::forward<Args>(args)...);
        end = clock();
        avg_time += double(end - start);
    }
    avg_time /= iters;
    return avg_time / CLOCKS_PER_SEC;
}
#endif // UTILS_H