#ifndef PARALLEL_NATALI_H    // Include guard to prevent multiple inclusions
#define PARALLEL_NATALI_H

#include <vector>         // For std::vector
#include <complex>        // For std::complex

// Define a complex number type alias for convenience
using complex_num = std::complex<long double>;

// Function declarations
std::vector<complex_num> FFT(std::vector<complex_num> input);

#endif // PARALLEL_NATALI_H
