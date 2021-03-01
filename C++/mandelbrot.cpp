/*
 * C++ Mandelbrot.
 *
 * Usage: ./build/C++/mandelbrot_cpp <image width> <image height> <max iterations> <repetitions (1+)> <center x> <center y> <section height> <gradient filename> <output filename>
 * Example: ./build/C++/mandelbrot_cpp 800 600 200 1 -0.8 0.0 2.2 gradients/benchmark.gradient mandelbrot.raw
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
#include <sstream>
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

// Compare two float values for "enough" equality.
bool equal_enough(float a, float b) noexcept
{
    constexpr float epsilon = std::numeric_limits<float>::epsilon();

    a = std::abs(a);
    b = std::abs(b);

    return std::abs(a - b) <= std::max(a, b) * epsilon;
}

struct PixelColor {
    unsigned char r, g, b;
};

struct GradientColor {
    float pos;
    float r, g, b;
    bool operator==(const float p) const noexcept { return equal_enough(p, pos); }
    bool operator<(const GradientColor& col) const noexcept { return pos < col.pos; }
};

struct Gradient {
    std::vector<GradientColor> colors;
};

void die(Error error)
{
    std::cerr << "Error: " << static_cast<int>(error) << std::endl;
    std::exit(static_cast<int>(error));
}

Gradient load_gradient(const std::string& filename)
{
    Gradient gradient;
    gradient.colors.push_back(GradientColor{0.0f, 0.0f, 0.0f, 0.0f});
    gradient.colors.push_back(GradientColor{1.0f, 1.0f, 1.0f, 1.0f});

    std::ifstream in{filename};
    std::string line;

    if (!in.is_open())
        die(Error::LoadGradient);

    const std::regex re{R"(([0-9]*\.?[0-9]+):\s*([0-9]*\.?[0-9]+),\s*([0-9]*\.?[0-9]+),\s*([0-9]*\.?[0-9]+))"};
    std::smatch m;

    while (std::getline(in, line)) {
        if (std::regex_match(line, m, re)) {
            if (m.size() == 4+1) {
                float pos = std::stof(m[1]);
                auto col = std::find(gradient.colors.begin(), gradient.colors.end(), pos);

                if (col != gradient.colors.end())
                    *col = GradientColor{pos, std::stof(m[2]), std::stof(m[3]), std::stof(m[4])};
                else
                    gradient.colors.push_back({pos, std::stof(m[2]), std::stof(m[3]), std::stof(m[4])});
            }
        }
    }

    std::sort(gradient.colors.begin(), gradient.colors.end());
    return gradient;
}

void color_from_gradient_range(const GradientColor& left, const GradientColor& right, const float pos, PixelColor& pixel_color) noexcept
{
    const float pos2 = (pos - left.pos) / (right.pos - left.pos);

    pixel_color.r = static_cast<unsigned char>(255.0f * (((right.r - left.r) * pos2) + left.r));
    pixel_color.g = static_cast<unsigned char>(255.0f * (((right.g - left.g) * pos2) + left.g));
    pixel_color.b = static_cast<unsigned char>(255.0f * (((right.b - left.b) * pos2) + left.b));
}

void color_from_gradient(const Gradient& gradient, const float pos, PixelColor& pixel_color) noexcept
{
    const auto end = gradient.colors.cend();

    const auto it = std::adjacent_find(gradient.colors.cbegin(), end,
        [&](const GradientColor& left, const GradientColor& right) { return left.pos <= pos && pos <= right.pos; });

    if (it != end)
        color_from_gradient_range(*it, *(it+1), pos, pixel_color);
}

void mandelbrot_calc(const int image_width, const int image_height, const int max_iterations, const double center_x, const double center_y, const double height,
                     std::vector<int>& histogram, std::vector<int>& iterations_per_pixel, std::vector<float>& smoothed_distances_to_next_iteration_per_pixel) noexcept
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
        const double y0 = y_top - height * (static_cast<double>(pixel_y) / static_cast<double>(image_height));

        for (int pixel_x = 0; pixel_x < image_width; ++pixel_x) {
            const double x0 = x_left + width * (static_cast<double>(pixel_x) / static_cast<double>(image_width));

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
                smoothed_distances_to_next_iteration_per_pixel[static_cast<std::size_t>(pixel)] = 1.0f - std::min(1.0f, static_cast<float>((std::log(std::log(final_magnitude)) - log_log_bailout) / log_2));
                ++histogram[static_cast<std::size_t>(iter)];  // no need to count histogram[max_iterations]
            }

            iterations_per_pixel[static_cast<std::size_t>(pixel)] = iter;  // 1 .. max_iterations
        }
    }
}

void mandelbrot_colorize(const int max_iterations, const Gradient& gradient,
                         std::vector<PixelColor>& image_data, const std::vector<int>& histogram, const std::vector<int>& iterations_per_pixel, const std::vector<float>& smoothed_distances_to_next_iteration_per_pixel, std::vector<float>& normalized_colors) noexcept
{
    // Sum all iterations, not counting the last one at position histogram[max_iterations] (which
    // are points in the Mandelbrot Set).
    const float total_iterations = static_cast<float>(std::accumulate(std::next(histogram.cbegin()), std::prev(histogram.cend()), 0));

    // Normalize the colors (0.0 .. 1.0) based on how often they are used in the image, not counting
    // histogram[max_iterations] (which are points in the Mandelbrot Set).
    int running_total = 0;

    for (std::size_t i = 1; i < static_cast<std::size_t>(max_iterations); ++i) {
        running_total += histogram[i];
        normalized_colors[i] = static_cast<float>(running_total) / total_iterations;
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
            const float color_of_previous_iter = normalized_colors[static_cast<std::size_t>(*iter - 1)];
            const float color_of_current_iter  = normalized_colors[static_cast<std::size_t>(*iter)];
            const float pos_in_gradient = color_of_previous_iter + *smoothed_distance_to_next_iteration * (color_of_current_iter - color_of_previous_iter);

            color_from_gradient(gradient, pos_in_gradient, pixel);
        }

        ++iter;
        ++smoothed_distance_to_next_iteration;
    }
}

bool save_image(const std::string& filename, const std::vector<PixelColor>& image_data)
{
    std::ofstream out{filename, std::ofstream::binary};

    if (!out.is_open())
        return false;

    for (auto p : image_data)
        out << p.r << p.g << p.b;

    return true;
}

double mean(const std::vector<double>& values)
{
    return std::accumulate(values.begin(), values.end(), 0.0) / static_cast<double>(values.size());
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

template <typename T>
T eval_arg(const char* s, T min, T max)
{
    T value;
    std::istringstream in{s};

    if (!(in >> value) || value < min || value > max)
        die(Error::EvalArgs);

    return value;
}

std::tuple<int, int, int, int, double, double, double, std::string, std::string> eval_args(const int argc, char const* argv[])
{
    if (argc != 10)
        die(Error::EvalArgs);

    auto image_width    = eval_arg(argv[1], 1, 100000);
    auto image_height   = eval_arg(argv[2], 1, 100000);
    auto max_iterations = eval_arg(argv[3], 1, 1000000000);
    auto repetitions    = eval_arg(argv[4], 1, 1000000);
    auto center_x       = eval_arg(argv[5], -100.0, 100.0);
    auto center_y       = eval_arg(argv[6], -100.0, 100.0);
    auto height         = eval_arg(argv[7], -100.0, 100.0);
    auto colors         = std::string{argv[8]};
    auto filename       = std::string{argv[9]};

    return std::make_tuple(image_width, image_height, max_iterations, repetitions, center_x, center_y, height, colors, filename);
}

void go(const int image_width, const int image_height, const int max_iterations, const double center_x, const double center_y, const double height, const Gradient& gradient, std::vector<PixelColor>& image_data, std::vector<double>& durations, const int repetitions) noexcept
{
    // histogram & normalized_colors: for simplicity we only use indices [1] .. [max_iterations], [0] is unused
    std::vector<int> histogram(static_cast<std::size_t>(max_iterations + 1));
    std::vector<float> normalized_colors(static_cast<std::size_t>(max_iterations + 1));
    std::vector<int> iterations_per_pixel(static_cast<std::size_t>(image_width * image_height));
    std::vector<float> smoothed_distances_to_next_iteration_per_pixel(static_cast<std::size_t>(image_width * image_height));

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
