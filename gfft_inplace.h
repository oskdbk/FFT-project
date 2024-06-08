#ifndef GFFT_INPLACE_H
#define GFFT_INPLACE_H

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

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef long double num;
typedef std::vector<long double>::const_iterator NumIter;

typedef std::complex<num> complex_num;
using namespace std;

vector<complex_num> GeneralFFT_inplace(vector<complex_num> P, bool f = false);

#endif // GFFT_INPLACE_H