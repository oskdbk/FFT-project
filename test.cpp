#include "utils.h"  // Include the header file containing measure_time
#include "general/gfft.h"    // Include the header file or declare the GeneralFFT function
#include "general/gfft_parallel.h"
#include "radix/FFT.h"
#include "radix/radix_parallel.h"

// Sample data
std::vector<complex_num> sample_data = {
    {1.0, 2.0}, {2.0, 3.0}, {3.0, 4.0}, {4.0, 5.0}, {5.0, 6.0}, {6.0, 7.0}, {7.0, 8.0}, {8.0, 9.0}
};

int main() {
    // Number of iterations for timing
    // int iterations = 1000;

    // Measure the time taken by GeneralFFT function
    // double time_taken = measure_time(GeneralFFT, iterations, sample_data, false);

    std::vector<complex_num> vec = sample_data;

    std::vector<complex_num> vec_1 = Radix2FFT(vec);
    PRT1(vec_1, "Radix2FFT");
    return 0;
    std::vector<complex_num> vec_2 = InverseRadix2FFT(vec_1);
    std::vector<complex_num> pvec_1 = Radix2FFT_parallel(vec, 8, false);
    std::vector<complex_num> pvec_2 = Radix2FFT_parallel(pvec_1, 8, true);

    // test on 10 random vectors of different lengths 2^p
    for (int i = 1; i < 10; i++) {

        int n = pow(2, i);
        cout << "Testing on vector of length " << n << endl;
        
        // generate random vector of length n
        vector<complex_num> vec(n);
        for (int j = 0; j < n; j++) {
            vec[j] = complex_num(rand() % 100, rand() % 100);
        }

        std::vector<complex_num> vec_1 = DFT(vec);
        std::vector<complex_num> vec_2 = InverseRadix2FFT(vec_1);
        std::vector<complex_num> pvec_1 = Radix2FFT_parallel(vec, 8, false);
        std::vector<complex_num> pvec_2 = Radix2FFT_parallel(pvec_1, 8, true);
        if (!is_same_vector(vec_1, pvec_1, "Radix2FFT", "Radix2FFT_parallel", false)){
            std::cout << "Failed on vector of length " << n << " for Radix2FFT" << std::endl;
            break;
        }
        if (!is_same_vector(vec_2, pvec_2, "InverseRadix2FFT", "InverseRadix2FFT_parallel", false)){
            std::cout << "Failed on vector of length " << n << " for InverseRadix2FFT" << std::endl;
            break;
        }
    }
    return 0;
}