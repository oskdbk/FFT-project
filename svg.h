#ifndef SVG_H
#define SVG_H

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <complex>
#include <regex>

namespace svg {

std::vector<std::complex<double>> getpoints(const std::string& filepath, int num_samples) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file");
    }

    std::string line;
    std::stringstream buffer;
    buffer << file.rdbuf();
    std::string content = buffer.str();
    
    // Regular expression to extract path data
    std::regex path_regex(R"(d\s*=\s*\"([^\"]+)\")");
    std::smatch match;
    std::string path_data;

    if (std::regex_search(content, match, path_regex) && match.size() > 1) {
        path_data = match.str(1);
    } else {
        throw std::runtime_error("No path data found in SVG");
    }

    // Regular expression to extract the commands and points from path data
    std::regex coord_regex(R"(([MLHVCSQTAZ])\s*([^MLHVCSQTAZ]*)?)");
    std::smatch coord_match;
    std::vector<std::complex<double>> points;

    auto it = path_data.begin();
    while (std::regex_search(it, path_data.end(), coord_match, coord_regex)) {
        char command = coord_match.str(1)[0];
        std::string coords = coord_match.str(2);
        
        std::istringstream coord_stream(coords);
        double x, y;

        switch (command) {
            case 'M':
            case 'L':
                while (coord_stream >> x >> y) {
                    points.emplace_back(x, y);
                }
                break;
            case 'Z':
                if (!points.empty()) {
                    points.push_back(points.front());
                }
                break;
            // Add handling for other SVG path commands as necessary
            default:
                break;
        }

        it = coord_match.suffix().first;
    }

    if (points.size() > num_samples) {
        points.resize(num_samples);
    }

    return points;
}

} // namespace svg

#endif // SVG_H
