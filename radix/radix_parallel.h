#include <thread>
#include <numeric>
#include <vector>
#include <complex> // to use complex numbers
#include <cmath>
#include <iostream>

typedef double num;
typedef std::vector<double>::const_iterator NumIter;

typedef std::complex<num> complex_num;

using namespace std;

vector<complex_num> Radix2FFT_parallel(vector<complex_num> P, int num_threads, bool inverse = false);