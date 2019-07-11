/*
 * C++ Mandelbrot.
 *
 * Usage: mandelbrot_cpp <image_width> <image_height> <max_iterations> <repetitions (1+)> <center x> <center y> <section height> <gradient filename> <output filename>
 * Example: mandelbrot_cpp 320 200 20 1 -0.5 0.0 2.0 blue.gradient mandelbrot.raw
 *
 * Gradient file example:
 *   0.0: 0.0, 0.0, 0.0
 *   0.5: 0.0, 0.0, 1.0
 *   1.0: 1.0, 1.0, 1.0
 *
 * Tobias Brückner, 2019
 */

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <regex>
#include <string>
#include <tuple>
#include <vector>

enum class Error {
    AllocMemory = 1,
    EvalArgs,
    LoadGradient,
    SaveImage,
    GetTime
};

struct PixelColor {
    unsigned char r, g, b;
};

struct GradientColor {
    double pos;
    double r, g, b;
};

struct Gradient {
    std::vector<GradientColor> colors;
};


void die(Error error)
{
    std::cerr << "Error: " << static_cast<int>(error) << std::endl;
    std::exit(static_cast<int>(error));
}

// Compare two double values for "enough" equality.
bool equal_enough(double a, double b)
{
    constexpr double epsilon = std::numeric_limits<double>::epsilon();

    a = std::abs(a);
    b = std::abs(b);

    return std::abs(a - b) <= std::max(a, b) * epsilon;
}

GradientColor* gradient_get_color_at_position(Gradient& gradient, const double pos)
{
    for (GradientColor& col : gradient.colors)
        if (equal_enough(col.pos, pos))
            return &col;

    return nullptr;
}

Gradient load_gradient(const std::string& filename)
{
    Gradient gradient;
    gradient.colors.push_back(GradientColor{0.0, 0.0, 0.0, 0.0});
    gradient.colors.push_back(GradientColor{1.0, 1.0, 1.0, 1.0});

    std::ifstream in{filename};
    std::string line;

    if (!in.is_open())
        die(Error::LoadGradient);

    const std::regex re{R"(([0-9]*\.?[0-9]+):\s*([0-9]*\.?[0-9]+),\s*([0-9]*\.?[0-9]+),\s*([0-9]*\.?[0-9]+))"};
    std::smatch m;

    while (std::getline(in, line)) {
        if (std::regex_match(line, m, re)) {
            if (m.size() == 4+1) {
                double pos = std::stod(m[1]);
                GradientColor* col;

                if ((col = gradient_get_color_at_position(gradient, pos))) {
                    col->r = std::stod(m[2]);
                    col->g = std::stod(m[3]);
                    col->b = std::stod(m[4]);
                } else {
                    GradientColor col{pos, std::stod(m[2]), std::stod(m[3]), std::stod(m[4])};
                    gradient.colors.push_back(col);
                }
            }
        }
    }

    std::sort(gradient.colors.begin(), gradient.colors.end(), [](const auto& col1, const auto& col2) { return col1.pos < col2.pos; });
    return gradient;
}

void color_from_gradient_range(const GradientColor& left, const GradientColor& right, const double pos, PixelColor& pixel_color)
{
    const double pos2 = (pos - left.pos) / (right.pos - left.pos);

    pixel_color.r = static_cast<unsigned char>(255.0 * (((right.r - left.r) * pos2) + left.r));
    pixel_color.g = static_cast<unsigned char>(255.0 * (((right.g - left.g) * pos2) + left.g));
    pixel_color.b = static_cast<unsigned char>(255.0 * (((right.b - left.b) * pos2) + left.b));
}

void color_from_gradient(const Gradient& gradient, double pos, PixelColor& pixel_color)
{
    if (pos < 0.0)
        pos = 0.0;

    if (pos > 1.0)
        pos = 1.0;

    for (std::size_t i = 1; i < gradient.colors.size(); ++i) {
        const GradientColor& left = gradient.colors[i - 1];
        const GradientColor& right = gradient.colors[i];

        if (pos >= left.pos && pos <= right.pos) {
            color_from_gradient_range(left, right, pos, pixel_color);
            break;
        }
    }
}

void mandelbrot_calc(const int image_width, const int image_height, const int max_iterations, const double center_x, const double center_y, const double height,
                     std::vector<int>& histogram, std::vector<int>& iterations_per_pixel, std::vector<double>& smoothed_distances_to_next_iteration_per_pixel)
{
    const double width = height * (static_cast<double>(image_width) / static_cast<double>(image_height));

    const double x_left   = center_x - width / 2.0;
 // const double x_right  = center_x + width / 2.0;
    const double y_top    = center_y + height / 2.0;
 // const double y_bottom = center_y - height / 2.0;

    constexpr double bailout = 20.0;
    constexpr double bailout_squared = bailout * bailout;
    const double log_log_bailout = std::log(std::log(bailout));
    const double log_2 = std::log(2.0);

    double final_magnitude = 0.0;

    std::fill(histogram.begin(), histogram.end(), 0);

    for (int pixel_y = 0; pixel_y < image_height; ++pixel_y) {
        for (int pixel_x = 0; pixel_x < image_width; ++pixel_x) {
            const double x0 = x_left + width * (static_cast<double>(pixel_x) / static_cast<double>(image_width));
            const double y0 = y_top - height * (static_cast<double>(pixel_y) / static_cast<double>(image_height));

            double x = 0.0;
            double y = 0.0;

            // iteration, will be from 1 to max_iterations once the loop is done
            int iter = 0;

            while (iter < max_iterations) {
                const double x_squared = x*x;
                const double y_squared = y*y;

                if (x_squared + y_squared >= bailout_squared) {
                    final_magnitude = std::sqrt(x_squared + y_squared);
                    break;
                }

                const double xtemp = x_squared - y_squared + x0;
                y = 2.0*x*y + y0;
                x = xtemp;

                ++iter;
            }

            const int pixel = pixel_y * image_width + pixel_x;

            if (iter < max_iterations) {
                smoothed_distances_to_next_iteration_per_pixel[static_cast<std::size_t>(pixel)] = 1.0 - std::min(1.0, (std::log(std::log(final_magnitude)) - log_log_bailout) / log_2);
                ++histogram[static_cast<std::size_t>(iter)];  // no need to count histogram[max_iterations]
            }

            iterations_per_pixel[static_cast<std::size_t>(pixel)] = iter;  // 1 .. max_iterations
        }
    }
}

void mandelbrot_colorize(const int max_iterations, const Gradient& gradient,
                         std::vector<PixelColor>& image_data, const std::vector<int>& histogram, const std::vector<int>& iterations_per_pixel, const std::vector<double>& smoothed_distances_to_next_iteration_per_pixel, std::vector<double>& normalized_colors)
{
    // Sum all iterations, not counting the last one at position histogram[max_iterations] (which
    // are points in the Mandelbrot Set).
    const double total_iterations = std::accumulate(std::next(histogram.cbegin()), std::prev(histogram.cend()), 0);

    // Normalize the colors (0.0 .. 1.0) based on how often they are used in the image, not counting
    // histogram[max_iterations] (which are points in the Mandelbrot Set).
    int running_total = 0;

    for (std::size_t i = 1; i < static_cast<std::size_t>(max_iterations); ++i) {
        running_total += histogram[i];
        normalized_colors[i] = running_total / total_iterations;
    }

    auto iter = iterations_per_pixel.cbegin();  // in range of 1 .. max_iterations
    auto smoothed_distance_to_next_iteration = smoothed_distances_to_next_iteration_per_pixel.cbegin();  // in range of 0 .. <1.0

    for (auto& pixel : image_data) {
        if (*iter == max_iterations) {
            // pixels with max. iterations (aka. inside the Mandelbrot Set) are always black
            pixel.r = 0;
            pixel.g = 0;
            pixel.b = 0;
        } else {
            // we use the color of the previous iteration in order to cover the full gradient range
            const double color_of_previous_iter = normalized_colors[static_cast<std::size_t>(*iter - 1)];
            const double color_of_current_iter  = normalized_colors[static_cast<std::size_t>(*iter)];
            const double pos_in_gradient = color_of_previous_iter + *smoothed_distance_to_next_iteration * (color_of_current_iter - color_of_previous_iter);

            color_from_gradient(gradient, pos_in_gradient, pixel);
        }

        ++iter;
        ++smoothed_distance_to_next_iteration;
    }
}

bool save_image(const std::string& filename, const std::vector<PixelColor>& image_data)
{
    std::ofstream out{filename};

    if (!out.is_open())
        return false;

    for (auto p : image_data)
        out << p.r << p.g << p.b;

    return true;
}

double mean(const std::vector<double>& values)
{
    return std::accumulate(values.begin(), values.end(), 0.0) / values.size();
}

double median(const std::vector<double>& values)
{
    std::vector<double> sorted_values{values};
    std::sort(sorted_values.begin(), sorted_values.end());

    if (sorted_values.size() % 2)
        return sorted_values[(sorted_values.size() - 1) / 2];
    else
        return (sorted_values[sorted_values.size() / 2 - 1] + sorted_values[sorted_values.size() / 2]) / 2.0;
}

void show_summary(const std::vector<double>& durations)
{
    if (durations.size() == 1) {
        std::cout << durations[0] << " s" << std::endl;
    } else {
        std::vector<double> sorted_values{durations};
        std::sort(sorted_values.begin(), sorted_values.end());

        std::cout << "mean: " << mean(sorted_values) << " s, median: " << median(sorted_values) << " s (repetitions=" << sorted_values.size() << ") ["
                  << sorted_values[0];

        for (auto p = std::next(sorted_values.cbegin()); p != sorted_values.cend(); ++p)
            std::cout << ", " << *p;

        std::cout << "]" << std::endl;
    }
}

int eval_int_arg(const char* s, int min, int max)
{
    int value;

    try {
        value = std::stoi(s);

        if (value < min || value > max)
            die(Error::EvalArgs);

        return value;
    } catch (...) {
        die(Error::EvalArgs);
    }

    return value;
}

double eval_double_arg(const char* s, double min, double max)
{
    double value;

    try {
        value = std::stod(s);

        if (value < min || value > max || std::isnan(value) || std::isinf(value))
            die(Error::EvalArgs);

        return value;
    } catch (...) {
        die(Error::EvalArgs);
    }

    return value;
}

std::tuple<int, int, int, int, double, double, double, std::string, std::string> eval_args(const int argc, char const* argv[])
{
    if (argc < 10)
        die(Error::EvalArgs);

    auto image_width    = eval_int_arg(argv[1], 1, 100000);
    auto image_height   = eval_int_arg(argv[2], 1, 100000);
    auto max_iterations = eval_int_arg(argv[3], 1, 1000000000);
    auto repetitions    = eval_int_arg(argv[4], 1, 1000000);
    auto center_x       = eval_double_arg(argv[5], -100.0, 100.0);
    auto center_y       = eval_double_arg(argv[6], -100.0, 100.0);
    auto height         = eval_double_arg(argv[7], -100.0, 100.0);
    auto colors         = std::string{argv[8]};
    auto filename       = std::string{argv[9]};

    return std::make_tuple(image_width, image_height, max_iterations, repetitions, center_x, center_y, height, colors, filename);
}

void go(const int image_width, const int image_height, const int max_iterations, const double center_x, const double center_y, const double height, const Gradient& gradient, std::vector<PixelColor>& image_data, std::vector<double>& durations, const int repetitions)
{
    // histogram & normalized_colors: for simplicity we only use indices [1] .. [max_iterations], [0] is unused
    std::vector<int> histogram(static_cast<std::size_t>(max_iterations + 1));
    std::vector<double> normalized_colors(static_cast<std::size_t>(max_iterations + 1));
    std::vector<int> iterations_per_pixel(static_cast<std::size_t>(image_width * image_height));
    std::vector<double> smoothed_distances_to_next_iteration_per_pixel(static_cast<std::size_t>(image_width * image_height));

    for (int i = 0; i < repetitions; ++i) {
        const auto t1 = std::chrono::high_resolution_clock::now();
        mandelbrot_calc(image_width, image_height, max_iterations, center_x, center_y, height, histogram, iterations_per_pixel, smoothed_distances_to_next_iteration_per_pixel);
        mandelbrot_colorize(max_iterations, gradient, image_data, histogram, iterations_per_pixel, smoothed_distances_to_next_iteration_per_pixel, normalized_colors);
        const auto t2 = std::chrono::high_resolution_clock::now();

        durations.push_back(std::chrono::duration<double>{t2 - t1}.count());
    }
}

int main(int argc, char const* argv[])
{
    auto [image_width, image_height, max_iterations, repetitions, center_x, center_y, height, gradient_filename, filename] = eval_args(argc, argv);
    auto gradient = load_gradient(gradient_filename);

    std::vector<PixelColor> image_data(static_cast<unsigned long>(image_width * image_height));
    std::vector<double> durations;

    go(image_width, image_height, max_iterations, center_x, center_y, height, gradient, image_data, durations, repetitions);

    if (!save_image(filename, image_data))
        die(Error::SaveImage);

    show_summary(durations);
}
