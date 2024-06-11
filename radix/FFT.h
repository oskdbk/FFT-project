#include <iostream>
#include <climits>
#include <thread>
#include <numeric>
#include <iterator>
#include <optional>
#include <vector>
#include <complex> // to use complex numbers
#include <cmath>

typedef double num;
typedef std::vector<double>::const_iterator NumIter;

typedef std::complex<num> complex_num;
using namespace std;
vector<complex_num> Radix2FFT(vector<complex_num> P);
vector<complex_num> InverseRadix2FFT(vector<complex_num> P_star);
