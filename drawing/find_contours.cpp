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

vector<complex_num> contour_to_complex_vector(const std::vector<Point>& contour) {
    size_t N = contour.size();
    std::vector<complex_num> complex_vector(N);

    for (size_t i = 0; i < N; i++){
        auto point = contour[i];
        complex_vector[i] = complex_num(point.x, point.y);
    }
    return complex_vector;
}

/*
    Takes a contour and returns a complex vector: each Point(x, y) => complex_num(x, y)
*/
vector<complex_num> contour_to_complex_vector_x(const std::vector<Point>& contour) {
    size_t N = contour.size();
    std::vector<complex_num> complex_vector(N);

    for (size_t i = 0; i < N; i++){
        auto point = contour[i];
        complex_vector[i] = complex_num(point.x, 0);
    }
    return complex_vector;
}

vector<complex_num> contour_to_complex_vector_y(const std::vector<Point>& contour) {
    size_t N = contour.size();
    std::vector<complex_num> complex_vector(N);

    for (size_t i = 0; i < N; i++){
        auto point = contour[i];
        complex_vector[i] = complex_num(point.y, 0);
    }
    return complex_vector;
}

// To test drawing: parametric equations for a circle

// float parametricX(float t) {
//     return cos(t) + sin(t);
// }

// float parametricY(float t) {
//     return sin(t*5);
// }

/*
    Parametric equation using fft
 */
complex_num parametric(complex_num init_res, double rotation, float t, vector<complex_num> &fft) {
    complex_num res = init_res;
    size_t N = fft.size();

    for (int i = 0; i < N; i++){
        double amp = sqrt(fft[i].imag()*fft[i].imag() + fft[i].real()*fft[i].real());
        double freq = i;
        double phase = atan2(fft[i].imag(), fft[i].real());
        res += complex_num(amp*cos(freq*t + phase + rotation), amp*sin(freq*t + phase + rotation));
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

    // Make image gray-scale and use binary threshold
    Mat gray, binary;
    cvtColor(image, gray, COLOR_BGR2GRAY);
    threshold(gray, binary, 100, 255, THRESH_BINARY);

    // Find contours
    vector<vector<Point>> contours;
    vector<Vec4i> hierarchy;
    findContours(binary, contours, hierarchy, RETR_TREE, CHAIN_APPROX_SIMPLE);

    /*
    FFT
    */
    // Store generated points for drawing
    sf::VertexArray curve(sf::LineStrip);
    
    int WIDTH = 800;
    int HEIGHT = 600;
    float PI = 3.14159265358979323846f;

    for (int c = 0; c < 1; c++){
        auto contour = contours[c];
        // vector<complex_num> contour_complex = contour_to_complex_vector(contour);

        vector<complex_num> contour_complex_x = contour_to_complex_vector_x(contour);
        vector<complex_num> contour_complex_y = contour_to_complex_vector_y(contour);
        vector<complex_num> fft_x = GeneralFFT_Parallel(contour_complex_x);
        vector<complex_num> fft_y = GeneralFFT_Parallel(contour_complex_y);

        // // debugging
        // for (int i = 0; i < contour.size(); i++){
        //     float x = contour[i].x;//pnt.real()/500;  // Scale x
        //     float y = contour[i].y;//pnt.imag()/500;  // Scale y
        //     std::cout << x << "  " << y << std::endl;
        //     curve.append(sf::Vertex(sf::Vector2f(x, y), sf::Color::Red));
        // }

        // Calculate fourier coefficients and generate points
        for (float t = 0; t < 2*M_PI; t += 2*M_PI/fft_x.size()) {
            double x = parametric(complex_num(WIDTH/2.0 + 100, 100), 0, t, fft_x).real();
            double y = parametric(complex_num(100.0, HEIGHT/2.0+100), M_PI/2, t, fft_y).imag();
            x = x/300;
            y = y/300;
            std::cout << x << "  " << y << std::endl;
            curve.append(sf::Vertex(sf::Vector2f(x, y), sf::Color::Red));
        }
    }

    // // Draw contours
    // Mat contourOutput = Mat::zeros(image.size(), CV_8UC3);
    // for (size_t i = 0; i < contours.size(); i++) {
    //     Scalar color = Scalar(0, 255, 0); // Green color for contours
    //     drawContours(contourOutput, contours, (int)i, color, 2, LINE_8, hierarchy, 0);
    // }

    // for (const auto& contour : contours) {
    //     for (const auto& point : contour) {
    //         int x = point.x;
    //         int y = point.y;
    //         cout << "X: " << x << " Y: " << y << endl;
    //     }
    // }

    // Display the result
    // imshow("Original Image", image);
    // imshow("Contours", contourOutput);

    sf::RenderWindow window(sf::VideoMode(WIDTH, HEIGHT), "Parametric Curve");

    sf::View view = window.getView();
    view.setCenter(0.0f, 0.0f);
    view.setSize(WIDTH, HEIGHT); // flip the y axis
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