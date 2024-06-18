#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <string>
#include <complex>
#include <cmath>
#include <iostream>
#include <fstream>
#include <tuple>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;

typedef double num;
typedef std::complex<num> complex_num;

std::vector<complex_num> Read_CSV(std::string file_path);
void Write_CSV(const std::vector<complex_num>& data, const std::string& file_path);

template <typename T>
void Write_CSV_Columns(const std::vector<std::vector<T>>& data, const std::vector<std::string>& headers, const std::string& file_path) {
    std::ofstream file(file_path);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file");
    }
    
    // Write headers
    for (const auto& header : headers) {
        file << header << ",";
    }
    file << "\n";
    
    // Write data
    size_t num_columns = data.size();
    size_t num_rows = data[0].size();
    for (size_t i = 0; i < num_rows; ++i) {
        for (size_t j = 0; j < num_columns; ++j) {
            // if T = complex_num
            if constexpr (std::is_same<T, complex_num>::value) {
                file << std::to_string(data[j][i].real()) << ",";
            } else {
                file << std::to_string(data[j][i]) << ",";
            }
        }
        file << "\n";
    }
    
    file.close();
}

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

void keep_largest_n(std::vector<complex_num> &P, int n);
#endif // UTILS_H