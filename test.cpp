#include "utils.h"  // Include the header file containing measure_time
#include "gfft.h"    // Include the header file or declare the GeneralFFT function
#include "gfft_inplace.h"
#include "gfft_parallel.h"

// Sample data
std::vector<complex_num> sample_data = {
    {1.0, 2.0}, {2.0, 3.0}, {3.0, 4.0}, {4.0, 5.0}
    // Add more complex numbers if needed
};

int main() {
    // Number of iterations for timing
    // int iterations = 1000;

    // Measure the time taken by GeneralFFT function
    // double time_taken = measure_time(GeneralFFT, iterations, sample_data, false);


    std::vector<complex_num> vec = sample_data;

    std::vector<complex_num> vec_1 = GeneralFFT_inplace(vec, false);
    std::vector<complex_num> vec_2 = GeneralFFT_inplace(vec_1, false, true);

    std::vector<complex_num> pvec_1 = GeneralFFT_Parallel(vec, 8, false);
    std::vector<complex_num> pvec_2 = GeneralFFT_Parallel(pvec_1, 8, true);

    PRT1(vec_2, "Normal Vector");
    PRT1(pvec_2, "Parallel");
    return 0;
}