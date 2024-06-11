#include <iostream>
#include <vector>
#include <fstream>
#include <complex>
#include <cmath>
#include <fftw3.h>
#include <matplotlibcpp.h>
#include "svg.h" // Ensure you have a library or implementation for parsing SVG and extracting points.

namespace plt = matplotlibcpp;

const int num_samples = 2000;
const double sampling_frequency = 2000.0;
const double function_period = 2 * M_PI;

std::complex<double> f(double t) {
    return std::cos(t / 2) * std::sin(8 * t) * std::exp(std::complex<double>(0, t));
}

int main() {
    std::string ImagePath = "SampleImages/earth.svg";
    std::string OutputFileName = "output.txt";
    std::vector<std::complex<double>> samples;

    std::vector<double> x_p, y_p;
    for (const auto& sample : samples) {
        x_p.push_back(sample.real());
        y_p.push_back(sample.imag());
    }

    // Plot the complex numbers
    plt::scatter(x_p, y_p);
    plt::xlabel("Real");
    plt::ylabel("Imaginary");
    plt::save("Samples_in_complex_plane.png");
    plt::show();

    // FFT calculation
    fftw_complex *in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_samples);
    fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_samples);
    fftw_plan p = fftw_plan_dft_1d(num_samples, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    for (int i = 0; i < num_samples; ++i) {
        in[i][0] = samples[i].real();
        in[i][1] = samples[i].imag();
    }

    fftw_execute(p);

    std::vector<double> fourier_frequencies(num_samples);
    for (int i = 0; i < num_samples; ++i) {
        fourier_frequencies[i] = i * sampling_frequency / num_samples;
    }

    std::vector<double> magnitudes(num_samples);
    double max_abs = 0.0;
    for (int i = 0; i < num_samples; ++i) {
        double abs_val = std::sqrt(out[i][0] * out[i][0] + out[i][1] * out[i][1]);
        magnitudes[i] = abs_val;
        if (abs_val > max_abs) {
            max_abs = abs_val;
        }
    }

    // Normalize magnitudes
    for (int i = 0; i < num_samples; ++i) {
        magnitudes[i] /= max_abs;
    }

    // Plot Fourier Transform
    plt::plot(fourier_frequencies, magnitudes);
    plt::save("Fourier_Transform.png");
    plt::show();

    // Write to output file
    std::ofstream outfile(OutputFileName);
    for (int i = 0; i < num_samples; ++i) {
        double phase = std::atan2(out[i][1], out[i][0]);
        outfile << fourier_frequencies[i] << ";" << magnitudes[i] << ";" << phase << "\n";
    }
    outfile.close();

    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);

    return 0;
}