#include <opencv2/opencv.hpp>
#include <iostream>
#include <complex>
#include <SFML/Graphics.hpp>
#include <cmath>
#include <thread>
#include "../general/gfft_parallel.h"

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
    Takes a contour and returns a complex vector: each Point(x, y) => complex_num(x, 0) or complex_num(y, 0)
*/
void contour_to_complex_vector_x_helper(const vector<cv::Point>& contour, vector<complex_num>& result, size_t start, size_t end) {
    for (size_t i = start; i < end; i++) {
        auto point = contour[i];
        result[i] = complex_num(point.x, 0);
    }
}

vector<complex_num> contour_to_complex_vector_x(const vector<cv::Point>& contour, size_t num_threads = 8) {
    size_t N = contour.size();
    vector<complex_num> result(N);

    size_t chunk_size = (N + num_threads - 1) / num_threads;
    vector<thread> threads;

    for (size_t i = 0; i < num_threads; ++i) {
        size_t start = i * chunk_size;
        size_t end = (i == num_threads - 1) ? N : (i + 1) * chunk_size;

        threads.emplace_back([&contour, &result, start, end]() {
            contour_to_complex_vector_x_helper(contour, result, start, end);
        });
    }

    for (auto& t : threads) {
        t.join();
    }

    return result;
}

void contour_to_complex_vector_y_helper(const vector<cv::Point>& contour, vector<complex_num>& result, size_t start, size_t end) {
    for (size_t i = start; i < end; i++) {
        auto point = contour[i];
        result[i] = complex_num(point.y, 0);
    }
}

vector<complex_num> contour_to_complex_vector_y(const vector<cv::Point>& contour, size_t num_threads = 8) {
    size_t N = contour.size();
    vector<complex_num> result(N);

    size_t chunk_size = (N + num_threads - 1) / num_threads;
    vector<thread> threads;

    for (size_t i = 0; i < num_threads; ++i) {
        size_t start = i * chunk_size;
        size_t end = (i == num_threads - 1) ? N : (i + 1) * chunk_size;

        threads.emplace_back([&contour, &result, start, end]() {
            contour_to_complex_vector_y_helper(contour, result, start, end);
        });
    }

    for (auto& t : threads) {
        t.join();
    }

    return result;
}

/*
    Parametric equation using fft
 */
void compute_partial_results(const vector<complex_num>& fft, double rotation, float t, vector<complex_num>& partial_results, size_t start, size_t end) {
    size_t N = fft.size();
    for (size_t j = start; j < end; ++j) {
        double amp = sqrt(fft[j].imag()*fft[j].imag() + fft[j].real()*fft[j].real());
        double freq = j;
        double phase = atan2(fft[j].imag(), fft[j].real());
        partial_results[j] = complex_num(amp*cos(freq*t + phase + rotation), amp*sin(freq*t + phase + rotation));
    }
}
complex_num parametric(complex_num init_res, double rotation, float t, vector<complex_num>& fft, int num_threads = 8) {
    complex_num res = init_res;
    size_t N = fft.size();
    
    vector<complex_num> partial_results(N);
    size_t chunk_size = (N + num_threads - 1) / num_threads;

    vector<thread> threads;

    for (size_t i = 0; i < num_threads; ++i) {
        size_t start = i * chunk_size;
        size_t end = (i == num_threads - 1) ? N : (i + 1) * chunk_size;

        threads.emplace_back([&fft, &partial_results, rotation, t, start, end]() {
            compute_partial_results(fft, rotation, t, partial_results, start, end);
        });
    }

    for (auto& t : threads) {
        t.join();
    }

    for (size_t i = 0; i < N; ++i) {
        res += partial_results[i];
    }

    return res;
}

/*
    Scaling vertices and adding to final result 
*/
void construct_vertices(const vector<complex_num>& points, double xmin, double ymin, double scale, sf::VertexArray& curve, size_t start, size_t end) {
    for (size_t i = start; i < end; ++i) {
        double x = (points[i].real() - xmin) * scale / 4;
        double y = (points[i].imag() - ymin) * scale / 4;
        curve[i] = sf::Vertex(sf::Vector2f(x, y), sf::Color::Red);
    }
}

void parallel_construct_vertices(const vector<complex_num>& points, double xmin, double ymin, double scale, sf::VertexArray& curve, int num_threads = 8) {
    size_t num_points = points.size();
    size_t chunk_size = num_points / num_threads;

    vector<thread> threads;

    for (size_t i = 0; i < num_threads; ++i) {
        size_t start = i * chunk_size;
        size_t end = (i == num_threads - 1) ? num_points : (i + 1) * chunk_size;

        threads.emplace_back([&points, xmin, ymin, scale, &curve, start, end]() {
            construct_vertices(points, xmin, ymin, scale, curve, start, end);
        });
    }

    for (auto& t : threads) {
        t.join();
    }
}

int main(int argc, char* argv[]) {
    std::string filename;
    bool usePoints = false;
    int num_threads = 8; // Default number of threads

    // Parse command line arguments
    if (argc < 2 || argc > 5) {
        std::cerr << "Usage: " << argv[0] << " <filename> [-points] [-n <num_threads>]\n";
        return 1;
    } else {
        filename = argv[1];
        
        for (int i = 2; i < argc; ++i) {
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
            } else {
                std::cerr << "Unknown argument: " << arg << "\n";
                return 1;
            }
        }
    }

    // Load the image
    Mat image = imread("pictures/" + filename);
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
    int WIDTH = 800;
    int HEIGHT = 600;
    float PI = 3.14159265358979323846f;

    auto contour = combine_contours(contours);
    // Calculate Fourier transform
    vector<complex_num> contour_complex_x = contour_to_complex_vector_x(contour, num_threads);
    vector<complex_num> contour_complex_y = contour_to_complex_vector_y(contour, num_threads);
    vector<complex_num> fft_x = GeneralFFT_Parallel(contour_complex_x, num_threads);
    vector<complex_num> fft_y = GeneralFFT_Parallel(contour_complex_y, num_threads);

    // Generate points
    vector<complex_num> points;
    double xmin = INFINITY;
    double xmax = -INFINITY;
    double ymin = INFINITY;
    double ymax = -INFINITY;
    for (float t = 0; t < 2*M_PI; t += 2*M_PI/fft_x.size()) {
        double x = parametric(complex_num(WIDTH/2.0 + 100, 100), 0, t, fft_x, num_threads).real();
        double y = parametric(complex_num(100.0, HEIGHT/2.0+100), M_PI/2, t, fft_y, num_threads).imag();
        if (x < xmin) xmin = x;
        if (x > xmax) xmax = x;
        if (y < ymin) ymin = y;
        if (y > ymax) ymax = y;
        points.push_back(complex_num(x, y));
    }
    double scale_x = WIDTH / (xmax - xmin);
    double scale_y = HEIGHT / (ymax - ymin);
    double scale = std::min(scale_x, scale_y);

    sf::VertexArray curve;

    if (usePoints){
        curve = sf::VertexArray(sf::Points, points.size());
    } else {
        curve = sf::VertexArray(sf::LineStrip, points.size());
    }
    parallel_construct_vertices(points, xmin, ymin, scale, curve, num_threads);

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