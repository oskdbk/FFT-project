#include <opencv2/opencv.hpp>
#include <iostream>
#include <complex>
#include <SFML/Graphics.hpp>
#include <cmath>
#include "../gfft_parallel.h"

using namespace cv;
using namespace std;

typedef double num;
typedef std::complex<num> complex_num;

vector<complex_num> DFT(vector<complex_num>& x) {
    size_t N = x.size();
    vector<complex_num> X(N);
    for (size_t k = 0; k < N; k++) {
        for (size_t n = 0; n < N; n++) {
            double angle = -2 * M_PI * n * k / N;
            X[k] += x[n] * polar(1.0, angle);
        }
    }
    return X;
}

vector<complex_num> contour_to_complex_vector(const std::vector<Point>& contour) {
    size_t N = contour.size();
    std::vector<complex_num> complex_vector(N);
    for (size_t i = 0; i < N; i++) {
        auto point = contour[i];
        complex_vector[i] = complex_num(point.x, point.y);
    }
    return complex_vector;
}

complex_num parametric(float t, const vector<complex_num>& fft) {
    complex_num res = 0;
    size_t N = fft.size();
    for (size_t i = 0; i < N; i++) {
        double angle = -2 * i * M_PI * t / N;
        res += fft[i] * polar(1.0, angle);
    }
    return res;
}

int main() {
    // Load the image
    Mat image = imread("pi_try.png");
    if (image.empty()) {
        cout << "Could not open or find the image" << endl;
        return -1;
    }

    // Convert to grayscale and binary
    Mat gray, binary;
    cvtColor(image, gray, COLOR_BGR2GRAY);
    threshold(gray, binary, 100, 255, THRESH_BINARY);

    // Find contours
    vector<vector<Point>> contours;
    vector<Vec4i> hierarchy;
    findContours(binary, contours, hierarchy, RETR_TREE, CHAIN_APPROX_SIMPLE);

    // Check if any contours were found
    if (contours.empty()) {
        cout << "No contours found" << endl;
        return -1;
    }

    // Store generated points for drawing
    sf::VertexArray curve(sf::LineStrip);
    auto contour = contours[0];
    vector<complex_num> contour_complex = contour_to_complex_vector(contour);

    // Calculate Fourier coefficients
    vector<complex_num> fft = GeneralFFT_Parallel(contour_complex);

    // Generate points
    for (float t = 0; t <= 1; t += 0.01) {
        complex_num pnt = parametric(t, fft);
        float x = pnt.real() / 1000;  // Scale x
        float y = pnt.imag() / 1000;  // Scale y
        curve.append(sf::Vertex(sf::Vector2f(x, y), sf::Color::Red));
        cout << "Point: (" << x << ", " << y << ")" << endl;
    }

    // Display the result
    int WIDTH = 800;
    int HEIGHT = 600;

    sf::RenderWindow window(sf::VideoMode(WIDTH, HEIGHT), "Parametric Curve");
    sf::View view = window.getView();
    view.setCenter(0.0f, 0.0f);
    view.setSize(WIDTH, HEIGHT);
    window.setView(view);

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed) {
                window.close();
            }
        }
        window.clear(sf::Color::White);
        window.draw(curve);
        window.display();
    }

    return 0;
}
