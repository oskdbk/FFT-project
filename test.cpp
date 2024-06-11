#include "utils.h"  // Include the header file containing measure_time
#include "gfft.h"    // Include the header file or declare the GeneralFFT function

// Sample data
std::vector<complex_num> sample_data = {
    {1.0, 2.0}, {2.0, 3.0}, {3.0, 4.0}, {4.0, 5.0}
    // Add more complex numbers if needed
};

int main() {
    // Number of iterations for timing
    int iterations = 1000;

    // Measure the time taken by GeneralFFT function
    double time_taken = measure_time(GeneralFFT, iterations, sample_data, false);

    // Output the average time taken
    std::cout << "Average time taken by GeneralFFT: " << time_taken << " seconds" << std::endl;

    return 0;
}