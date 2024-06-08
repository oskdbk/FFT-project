#ifndef GFFT_PARALLEL_H    // Include guard to prevent multiple inclusions
#define GFFT_PARALLEL_H

#include <vector>         // For std::vector
#include <complex>        // For std::complex

// Define a complex number type alias for convenience
using complex_num = std::complex<long double>;

// Function declarations
std::vector<complex_num> GeneralFFT_Parallel(std::vector<complex_num> P, int num_threads = 1);

#endif // GFFT_PARALLEL_H
