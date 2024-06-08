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

typedef long double num;
typedef std::complex<num> complex_num;

std::vector<complex_num> Read_CSV(std::string file_path);
void PRT1(std::vector<std::complex<long double>> P, std::string a);
void PRT2(std::vector<std::vector<complex_num>> P, std::string a);
std::tuple<int, int> Decompose(int n);
std::vector<complex_num> StandardFFT(std::vector<complex_num> &P);

#endif // UTILS_H