#include <opencv2/opencv.hpp>
#include <iostream>
#include <complex>
#include <SFML/Graphics.hpp>
#include <cmath>
#include "../general/gfft_parallel.h"
#include "../utils.h"

using namespace cv;
using namespace std;

typedef double num;
typedef std::complex<num> complex_num;

/*
    Given contours, combine them and return just one vector
*/
vector<cv::Point_<int>> combine_contours(const std::vector<vector<cv::Point_<int>>> &contours){
    int total_size = 0;
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
vector<complex_num> contour_to_complex_vector(const std::vector<cv::Point>& contour) {
    size_t N = contour.size();
    std::vector<complex_num> complex_vector(N);

    for (size_t i = 0; i < N; i++){
        auto point = contour[i];
        complex_vector[i] = complex_num(point.x, point.y);
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

int main(int argc, char* argv[]) {
    bool usePoints = false;
    int num_threads = 1; // Default number of threads for FFT
    int top = 0; // Default value for top

    // Parse command line arguments
    if (argc < 1 || argc > 7) {
        std::cerr << "Usage: " << argv[0] << " [-points] [-n <num_threads>] [-top <value>]\n";
        return 1;
    } else {
        for (int i = 1; i < argc; ++i) {
            std::string arg = argv[i];
            if (arg == "-points") {
                usePoints = true;
            } else if (arg == "-n") {
                if (i + 1 < argc) {
                    num_threads = std::atoi(argv[++i]);
                    if (num_threads <= 0) {
                        std::cerr << "Error: Number of threads must be a positive integer.\n";
                        return 1;
                    }
                } else {
                    std::cerr << "Error: Missing value for -n option.\n";
                    return 1;
                }
            } else if (arg == "-top") {
                if (i + 1 < argc) {
                    top = std::atoi(argv[++i]);
                    if (top < 0) {
                        std::cerr << "Error: Value for -top must be a non-negative integer.\n";
                        return 1;
                    }
                } else {
                    std::cerr << "Error: Missing value for -top option.\n";
                    return 1;
                }
            } else {
                std::cerr << "Unknown argument: " << arg << "\n";
                return 1;
            }
        }
    }

    /*
    FFT
    */

    sf::VertexArray curve;

    if (usePoints){
        curve = sf::VertexArray(sf::Points);
    } else {
        curve = sf::VertexArray(sf::LineStrip);
    }
    
    int WIDTH = 800;
    int HEIGHT = 600;
    float PI = 3.14159265358979323846f;

    std::vector<cv::Point_<int>> contour;

    sf::RenderWindow window2(sf::VideoMode(800, 600), "Draw Curve");

    // Vector to store the coordinates
    std::vector<sf::Vector2f> coordinates;

    while (window2.isOpen()) {
        sf::Event event;
        while (window2.pollEvent(event)) {
            // Close window: exit
            if (event.type == sf::Event::Closed) {
                window2.close();
            }

            // Mouse button pressed: start drawing
            if (event.type == sf::Event::MouseButtonPressed) {
                if (event.mouseButton.button == sf::Mouse::Left) {
                    // Clear previous coordinates if needed
                    coordinates.clear();
                }
            }

            // Mouse moved: capture coordinates
            if (event.type == sf::Event::MouseMoved) {
                if (sf::Mouse::isButtonPressed(sf::Mouse::Left)) {
                    coordinates.push_back(sf::Vector2f(event.mouseMove.x, event.mouseMove.y));
                }
            }
        }

        // Clear the window
        window2.clear(sf::Color::White);

        // Draw the curve
        if (coordinates.size() > 1) {
            sf::VertexArray lines(sf::LineStrip, coordinates.size());
            for (size_t i = 0; i < coordinates.size(); ++i) {
                lines[i].position = coordinates[i];
                lines[i].color = sf::Color::Black;
            }
            window2.draw(lines);
        }

        // Update the window
        window2.display();
    }

    contour.resize(0);
    for (int i=0; i<coordinates.size(); i++){
        cv::Point_<int> point(coordinates[i].x, coordinates[i].y);
        contour.push_back(point);
    }



    cout << "Number of points is: " << contour.size() << endl;
    // Calculate Fourier transform
    vector<complex_num> contour_complex = contour_to_complex_vector(contour);
    vector<complex_num> fft = GeneralFFT_Parallel(contour_complex, num_threads);

    if (top > 0){
        keep_largest_n(fft, top);
    }

    // Generate points

    vector<complex_num> points;
    double xmin = INFINITY;
    double xmax = -INFINITY;
    double ymin = INFINITY;
    double ymax = -INFINITY;
    for (float t = 0; t < 2*M_PI; t += 4*M_PI/fft.size()) {
        auto pointt = parametric(complex_num(WIDTH/2.0 + 100, 100), 0, t, fft);
        double x = pointt.real();
        double y = pointt.imag();
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