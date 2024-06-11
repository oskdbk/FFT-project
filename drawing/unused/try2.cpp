#include <iostream>
#include <random>
#include <math.h>
#include <algorithm>
#include <stack>

#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"


class Vector {
public:
	explicit Vector(double x = 0, double y = 0, double z = 0) {
		data[0] = x;
		data[1] = y;
		data[2] = z;
	}
	double norm2() const {
		return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
	}
	double norm() const {
		return sqrt(norm2());
	}
	void normalize() {
		double n = norm();
		data[0] /= n;
		data[1] /= n;
		data[2] /= n;
	}
	double operator[](int i) const { return data[i]; };
	double& operator[](int i) { return data[i]; };
	double data[3];
};

Vector operator+(const Vector& a, const Vector& b) {
	return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b) {
	return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(const double a, const Vector& b) {
	return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector& a, const double b) {
	return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator/(const Vector& a, const double b) {
	return Vector(a[0] / b, a[1] / b, a[2] / b);
}
double dot(const Vector& a, const Vector& b) {
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector& a, const Vector& b) {
	return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

double norm (double a, double b, double c){
    return sqrt(a*a + b*b + c*c);
}
double norm (double a, double b){
    return sqrt(a*a + b*b);
}


struct Point {
    int x, y;
};

void floodFill(const std::vector<unsigned char> &edges, std::vector<bool> &visited, std::vector<Point> &contour, int x, int y, int W, int H) {
    std::stack<Point> stack;
    stack.push({x, y});
    while (!stack.empty()) {
        Point p = stack.top();
        stack.pop();
        if (p.x < 0 || p.x >= W || p.y < 0 || p.y >= H) continue;
        if (visited[p.y * W + p.x] || edges[p.y * W + p.x] == 0) continue;

        visited[p.y * W + p.x] = true;
        contour.push_back(p);

        stack.push({p.x + 1, p.y});
        stack.push({p.x - 1, p.y});
        stack.push({p.x, p.y + 1});
        stack.push({p.x, p.y - 1});
        stack.push({p.x + 1, p.y + 1});
        stack.push({p.x - 1, p.y + 1});
        stack.push({p.x + 1, p.y - 1});
        stack.push({p.x - 1, p.y - 1});
    }
}

std::vector<std::vector<Point>> findContours(const std::vector<unsigned char> &edges, int W, int H) {
    std::vector<std::vector<Point>> contours;
    std::vector<bool> visited(W * H, false);

    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            if (edges[y * W + x] == 255 && !visited[y * W + x]) {
                std::vector<Point> contour;
                floodFill(edges, visited, contour, x, y, W, H);
                if (!contour.empty()) {
                    contours.push_back(contour);
                }
            }
        }
    }

    return contours;
}

void drawLine(std::vector<unsigned char> &image, Point p1, Point p2, int W, int H) {
    int dx = abs(p2.x - p1.x), sx = p1.x < p2.x ? 1 : -1;
    int dy = -abs(p2.y - p1.y), sy = p1.y < p2.y ? 1 : -1;
    int err = dx + dy, e2;

    while (true) {
        if (p1.x >= 0 && p1.x < W && p1.y >= 0 && p1.y < H) {
            image[p1.y * W + p1.x] = 255;
        }
        if (p1.x == p2.x && p1.y == p2.y) break;
        e2 = 2 * err;
        if (e2 >= dy) { err += dy; p1.x += sx; }
        if (e2 <= dx) { err += dx; p1.y += sy; }
    }
}

int main() {

	int W, H, C;
	
	// Load input image
	unsigned char *image_input_raw = stbi_load("pi.png",
                                 &W,
                                 &H,
                                 &C,
                                 STBI_rgb);

	std::vector<double> input_image(W*H*3);
	for (int i=0; i<W*H*3; i++)
		input_image[i] = image_input_raw[i];

    
    // Make it black and white
    std::vector<double>image_grey(W*H);
    for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			image_grey[i*W + j] = norm(input_image[(i*W+j)*3+0], input_image[(i*W+j)*3+1], input_image[(i*W+j)*3+2]);
		}
	}
    // Sobel kernels
    int Gx[3][3] = {
        {-1, 0, 1},
        {-2, 0, 2},
        {-1, 0, 1}
    };

    int Gy[3][3] = {
        {-1, -2, -1},
        {0, 0, 0},
        {1, 2, 1}
    };

    // Compute gradients
    std::vector<double> grad_x(W*H, 0);
    std::vector<double> grad_y(W*H, 0);
    std::vector<double> grad_norm(W*H, 0);
    std::vector<double> grad_dir(W*H, 0);

    for (int i = 1; i < H - 1; i++) {
        for (int j = 1; j < W - 1; j++) {
            double gx = 0;
            double gy = 0;
            for (int ki = -1; ki <= 1; ki++) {
                for (int kj = -1; kj <= 1; kj++) {
                    gx += Gx[ki + 1][kj + 1] * image_grey[(i + ki) * W + (j + kj)];
                    gy += Gy[ki + 1][kj + 1] * image_grey[(i + ki) * W + (j + kj)];
                }
            }
            grad_x[i * W + j] = gx;
            grad_y[i * W + j] = gy;
            grad_norm[i * W + j] = std::sqrt(gx * gx + gy * gy);
            grad_dir[i * W + j] = std::atan2(gy, gx);
        }
    }

    // Non-maximum suppression
    std::vector<unsigned char> edges(W * H, 0);
    for (int i = 1; i < H - 1; i++) {
        for (int j = 1; j < W - 1; j++) {
            double angle = grad_dir[i * W + j] * 180.0 / M_PI;
            angle = angle < 0 ? angle + 180 : angle;
            
            double q = 255, r = 255;

            // Check pixels in the gradient direction
            if ((0 <= angle && angle < 22.5) || (157.5 <= angle && angle <= 180)) {
                q = grad_norm[i * W + (j + 1)];
                r = grad_norm[i * W + (j - 1)];
            } else if (22.5 <= angle && angle < 67.5) {
                q = grad_norm[(i + 1) * W + (j - 1)];
                r = grad_norm[(i - 1) * W + (j + 1)];
            } else if (67.5 <= angle && angle < 112.5) {
                q = grad_norm[(i + 1) * W + j];
                r = grad_norm[(i - 1) * W + j];
            } else if (112.5 <= angle && angle < 157.5) {
                q = grad_norm[(i - 1) * W + (j - 1)];
                r = grad_norm[(i + 1) * W + (j + 1)];
            }

            if (grad_norm[i * W + j] >= q && grad_norm[i * W + j] >= r) {
                edges[i * W + j] = static_cast<unsigned char>(grad_norm[i * W + j]);
            } else {
                edges[i * W + j] = 0;
            }
        }
    }

    // Threshold edges
    double threshold = 10.0;
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            if (edges[i * W + j] >= threshold) {
                edges[i * W + j] = 255;
            } else {
                edges[i * W + j] = 0;
            }
        }
    }

    // Find the nearest neighbor given a point
    // Algotithm:
    // Flag already visited
    // Find the closest one which is 1 in the radius
	
	
    // Find contours
    std::vector<std::vector<Point>> contours = findContours(edges, W, H);

    // Print contours size
    std:: cout << contours.size();

    /// Create a blank image to draw contours
    std::vector<unsigned char> contour_image(W * H, 0);

    // Draw contours
    for (const auto &contour : contours) {
        for (size_t i = 0; i < contour.size() - 1; ++i) {
            drawLine(contour_image, contour[i], contour[i + 1], W, H);
        }
        // Optionally, close the contour by drawing a line from the last point to the first point
        // drawLine(contour_image, contour.back(), contour.front(), W, H);
    }

    // Save the contour image
    stbi_write_png("contours.png", W, H, 1, contour_image.data(), W);

    // Free image memory
    stbi_image_free(image_input_raw);
	
	std::vector<unsigned char> image_result(W*H * 3, 0);
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {

			image_result[(i*W + j) * 3 + 0] = edges[(i*W+j)];
			image_result[(i*W + j) * 3 + 1] = edges[(i*W+j)];
			image_result[(i*W + j) * 3 + 2] = edges[(i*W+j)];
		}
	}
	stbi_write_png("image.png", W, H, 3, &image_result[0], 0);

	return 0;
}