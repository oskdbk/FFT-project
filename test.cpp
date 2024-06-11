#include "utils.h"  // Include the header file containing measure_time
#include "gfft.h"    // Include the header file or declare the GeneralFFT function
#include "gfft_inplace.h"

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


    std::vector<complex_num> vec(15, complex_num(1.0, 0.0));

    std::vector<complex_num> vec_1 = GeneralFFT_inplace(vec);
    std::vector<complex_num> vec_2 = GeneralFFT_inplace(vec_1, false, true);

    for(const auto& elem : vec) {
        std::cout << elem << std::endl;
    }

    std::cout << endl;

    for(const auto& elem : vec_1) {
        std::cout << elem << std::endl;
    }

    std::cout << endl;

    for(const auto& elem : vec_2) {
        std::cout << elem << std::endl;
    }

    // Output the average time taken
    return 0;
}