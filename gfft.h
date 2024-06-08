
#ifndef GFFT_H
#define GFFT_H

#include <iostream>
#include <vector>
#include <string>
#include <complex>
#include <tuple>

typedef long double num;
typedef std::complex<num> complex_num;

// Function declarations

/**
 * Compute the General FFT of a vector of complex numbers.
 * @param P A vector of complex numbers.
 * @param f A flag to enable/disable detailed output for debugging.
 * @return A vector containing the FFT result.
 */
std::vector<complex_num> GeneralFFT(std::vector<complex_num> &P, bool f = false);

/**
 * Read a CSV file containing a list of numbers and convert them to a vector of complex numbers.
 * @param file_path The path to the CSV file.
 * @return A vector of complex numbers read from the file.
 */
std::vector<complex_num> Read_CSV(std::string file_path);

#endif // GFFT_H