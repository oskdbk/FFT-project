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


vector<complex_num> DFT(vector<complex_num>&x) {
    size_t N = x.size();
    vector<complex_num> X = vector<complex_num>(N);
    int k, n;
    for(k = 0; k < N; k++) {
        for(n = 0; n < N; n++) {
            double angle = -2 * M_PI * n * k/ N;
            X[k] += x[n] * polar(1.0, angle);
        }
    }
    
    return X;
}

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
complex_num parametric(float t, vector<complex_num> &fft) {
    complex_num res = 0;
    size_t N = fft.size();
    // real
    // complex_num a0 = fft[0]/complex_num(N, 0);
    // vector<double> an(N-1);
    // vector<double> bn(N-1);

    // for (int i = 1; i < N; i++){
    //     an[i-1] = -2*(fft[i]/complex_num(N, 0)).real();
    //     bn[i-1] = -2*(fft[i]/complex_num(N, 0)).imag();
    // }

    // res += a0.real() / 2.0;

    // for (int n = 1; n < N; n++){
    //     res += an[n-1]*cos(n*t) + bn[n-1]*sin(n*t);
    // }

    for (int i = 0; i < fft.size(); i++){
        double angle = -2*i*M_PI*t/fft.size();
        complex_num cn = fft[i] * (complex_num(1, 0)*cos(angle) + complex_num(1, 1)*sin(angle));
        // res += cn*(complex_num(1, 0)*cos(angle) + complex_num(1, 1)*sin(angle));
        //res += fft[i]*polar(1.0, (i*2*M_PI*t)/fft.size());
        // res += fft[i]*polar(1.0, -(i*2*M_PI*t));
        double amp = (fft[i].imag()*fft[i].imag() + fft[i].real()*fft[i].real());
        double freq = i;
        double phase = atan2(fft[i].imag(), fft[i].real());
        res += complex_num(amp*cos(freq*t + phase + M_PI/2), amp*sin(freq*t + phase + M_PI/2));
        // res += amp*polar(1.0, (i*2*M_PI*t/fft.size()));
        //res += fft[i]*polar(1.0, -(i*2*M_PI*t));
        //res += fft[i]*polar(1.0, -(i*2*M_PI*t));
    }
    cout << "resultttttttttttttttttt" << res << endl;
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

        // Calculate fourier coefficients
        // vector<complex_num> fft = GeneralFFT_Parallel(contour_complex);
        // Generate points
        for (float t = 0; t <= 2*M_PI; t += 2*M_PI/fft.size()) {
            // double x = parametric(t, fft_x);
            // double y = parametric(t, fft_y);
            complex_num pnt = parametric(t, fft);
            float x = pnt.real()/20000000;  // Scale x
            float y = pnt.imag()/20000000;  // Scale y
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

    int WIDTH = 800;
    int HEIGHT = 600;
    float PI = 3.14159265358979323846f;

    sf::RenderWindow window(sf::VideoMode(WIDTH, HEIGHT), "Parametric Curve");

    sf::View view = window.getView();
    view.setCenter(0.0f, 0.0f);
    view.setSize(WIDTH, HEIGHT); // flip the y axis
    window.setView(view);

    // // Vertex array
    // sf::VertexArray curve(sf::LineStrip);

    // // Generate points
    // for (float t = 0; t <= 2 * PI; t += 0.01f) {
    //     float x = parametricX(t) * 200;  // Scale x
    //     float y = parametricY(t) * 200;  // Scale y
    //     curve.append(sf::Vertex(sf::Vector2f(x, y), sf::Color::Red));
    // }

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