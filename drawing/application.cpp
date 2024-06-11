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

/*
    Given contours, combine them and return just one vector
*/
vector<cv::Point_<int>> combine_contours(const std::vector<vector<cv::Point_<int>>> &contours){
    int total_size = 0;
    cout << total_size;
    for (auto contour: contours){
        total_size += contour.size();
    }
    vector<cv::Point_<int>>res;
    res.reserve(total_size);
    for (const auto& contour : contours) {
        res.insert(res.end(), contour.begin(), contour.end());
    }
    return res;
}

/*
    Takes a contour and returns a complex vector: each Point(x, y) => complex_num(x, y)
*/
vector<complex_num> contour_to_complex_vector_x(const std::vector<cv::Point>& contour) {
    size_t N = contour.size();
    std::vector<complex_num> complex_vector(N);

    for (size_t i = 0; i < N; i++){
        auto point = contour[i];
        complex_vector[i] = complex_num(point.x, 0);
    }
    return complex_vector;
}

vector<complex_num> contour_to_complex_vector_y(const std::vector<cv::Point>& contour) {
    size_t N = contour.size();
    std::vector<complex_num> complex_vector(N);

    for (size_t i = 0; i < N; i++){
        auto point = contour[i];
        complex_vector[i] = complex_num(point.y, 0);
    }
    return complex_vector;
}

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
    Mat image = imread("pictures/gojo.png");
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

    cout << "Number of contours is " << contours.size() << endl;

    /*
    FFT
    */
    // sf::VertexArray curve(sf::LineStrip);
    sf::VertexArray curve(sf::Points);
    
    int WIDTH = 800;
    int HEIGHT = 600;
    float PI = 3.14159265358979323846f;

    auto contour = combine_contours(contours);
    // Calculate Fourier transform
    vector<complex_num> contour_complex_x = contour_to_complex_vector_x(contour);
    vector<complex_num> contour_complex_y = contour_to_complex_vector_y(contour);
    vector<complex_num> fft_x = GeneralFFT_Parallel(contour_complex_x);
    vector<complex_num> fft_y = GeneralFFT_Parallel(contour_complex_y);

    // Generate points
    vector<complex_num> points;
    double xmin = INFINITY;
    double xmax = -INFINITY;
    double ymin = INFINITY;
    double ymax = -INFINITY;
    for (float t = 0; t < 2*M_PI; t += 2*M_PI/fft_x.size()) {
        double x = parametric(complex_num(WIDTH/2.0 + 100, 100), 0, t, fft_x).real();
        double y = parametric(complex_num(100.0, HEIGHT/2.0+100), M_PI/2, t, fft_y).imag();
        if (x < xmin) xmin = x;
        if (x > xmax) xmax = x;
        if (y < ymin) ymin = y;
        if (y > ymax) ymax = y;
        points.push_back(complex_num(x, y));
    }
    double scale_x = WIDTH / (xmax - xmin);
    double scale_y = HEIGHT / (ymax - ymin);
    double scale = std::min(scale_x, scale_y);

    for (size_t i = 0; i < points.size(); ++i) {
        double x = (points[i].real() - xmin) * scale/4;
        double y = (points[i].imag() - ymin) * scale/4;
        curve.append(sf::Vertex(sf::Vector2f(x, y), sf::Color::Red));
    }

    // Rendering
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