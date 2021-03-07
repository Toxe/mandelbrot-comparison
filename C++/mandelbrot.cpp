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
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

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

struct CalculationResult {
    int iter;
    float distance_to_next_iteration;
};

Gradient load_gradient(const std::string& filename)
{
    Gradient gradient;
    gradient.colors.push_back(GradientColor{0.0f, 0.0f, 0.0f, 0.0f});
    gradient.colors.push_back(GradientColor{1.0f, 1.0f, 1.0f, 1.0f});

    std::ifstream in{filename};
    std::string line;

    if (!in.is_open())
        throw std::runtime_error("unable to open gradient file");

    const std::regex re{R"(([0-9]*\.?[0-9]+):\s*([0-9]*\.?[0-9]+),\s*([0-9]*\.?[0-9]+),\s*([0-9]*\.?[0-9]+)$)"};
    std::smatch m;

    while (std::getline(in, line)) {
        if (std::regex_match(line, m, re)) {
            float pos = std::stof(m[1]);
            auto col = std::find(gradient.colors.begin(), gradient.colors.end(), pos);

            if (col != gradient.colors.end())
                *col = GradientColor{pos, std::stof(m[2]), std::stof(m[3]), std::stof(m[4])};
            else
                gradient.colors.push_back({pos, std::stof(m[2]), std::stof(m[3]), std::stof(m[4])});
        }
    }

    std::sort(gradient.colors.begin(), gradient.colors.end());
    return gradient;
}

void color_from_gradient_range(const GradientColor& left, const GradientColor& right, const float pos, PixelColor& pixel_color) noexcept
{
    const float relative_pos_between_colors = (pos - left.pos) / (right.pos - left.pos);
    pixel_color.r = static_cast<unsigned char>(255.0f * std::lerp(left.r, right.r, relative_pos_between_colors));
    pixel_color.g = static_cast<unsigned char>(255.0f * std::lerp(left.g, right.g, relative_pos_between_colors));
    pixel_color.b = static_cast<unsigned char>(255.0f * std::lerp(left.b, right.b, relative_pos_between_colors));
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
                     std::vector<int>& iterations_histogram, std::vector<CalculationResult>& results_per_point) noexcept
{
    const double width = height * (static_cast<double>(image_width) / static_cast<double>(image_height));

    const double x_left   = center_x - width / 2.0;
    const double x_right  = center_x + width / 2.0;
    const double y_top    = center_y + height / 2.0;
    const double y_bottom = center_y - height / 2.0;

    constexpr double bailout = 20.0;
    constexpr double bailout_squared = bailout * bailout;
    const double log_log_bailout = std::log(std::log(bailout));
    const double log_2 = std::log(2.0);

    double final_magnitude = 0.0;

    std::fill(iterations_histogram.begin(), iterations_histogram.end(), 0);

    int pixel = 0;

    for (int pixel_y = 0; pixel_y < image_height; ++pixel_y) {
        const double y0 = std::lerp(y_top, y_bottom, static_cast<double>(pixel_y) / static_cast<double>(image_height));

        for (int pixel_x = 0; pixel_x < image_width; ++pixel_x) {
            const double x0 = std::lerp(x_left, x_right, static_cast<double>(pixel_x) / static_cast<double>(image_width));

            double x = 0.0;
            double y = 0.0;

            // iteration, will be from 1 .. max_iterations once the loop is done
            int iter = 0;

            while (iter < max_iterations) {
                const double x_squared = x * x;
                const double y_squared = y * y;

                if (x_squared + y_squared >= bailout_squared) {
                    final_magnitude = std::sqrt(x_squared + y_squared);
                    break;
                }

                const double xtemp = x_squared - y_squared + x0;
                y = 2.0 * x * y + y0;
                x = xtemp;

                ++iter;
            }

            if (iter < max_iterations) {
                ++iterations_histogram[static_cast<std::size_t>(iter)]; // iter: 1 .. max_iterations-1, no need to count iterations_histogram[max_iterations]
                results_per_point[static_cast<std::size_t>(pixel)] = CalculationResult{iter, 1.0f - std::min(1.0f, static_cast<float>((std::log(std::log(final_magnitude)) - log_log_bailout) / log_2))};
            } else {
                results_per_point[static_cast<std::size_t>(pixel)] = CalculationResult{iter, 0.0};
            }

            ++pixel;
        }
    }
}

std::vector<float> equalize_histogram(const std::vector<int>& iterations_histogram, const int max_iterations)
{
    // Calculate the CDF (Cumulative Distribution Function) by accumulating all iteration counts.
    // Element [0] is unused and iterations_histogram[max_iterations] should be zero (as we do not count
    // the iterations of the points inside the Mandelbrot Set).
    std::vector<int> cdf(iterations_histogram.size());
    std::partial_sum(iterations_histogram.cbegin(), iterations_histogram.cend(), cdf.begin());

    // Get the minimum value in the CDF that is bigger than zero and the sum of all iteration counts
    // from iterations_histogram (which is the last value of the CDF).
    const auto cdf_min = std::find_if(cdf.cbegin(), cdf.cend(), [](auto n) { return n > 0; });
    const auto total_iterations = cdf[cdf.size() - 1];

    // normalize all values from the CDF that are bigger than zero to a range of 0.0 .. max_iterations
    const auto f = static_cast<float>(max_iterations) / static_cast<float>(total_iterations - *cdf_min);
    std::vector<float> equalized_iterations(iterations_histogram.size());
    std::transform(cdf.cbegin(), cdf.cend(), equalized_iterations.begin(),
                   [=](const auto& c) { return c > 0 ? f * static_cast<float>(c - *cdf_min) : 0.0; });

    return equalized_iterations;
}

void mandelbrot_colorize(const int max_iterations, const Gradient& gradient,
                         std::vector<PixelColor>& image_data, const std::vector<int>& iterations_histogram, const std::vector<CalculationResult>& results_per_point) noexcept
{
    const auto equalized_iterations = equalize_histogram(iterations_histogram, max_iterations);
    int pixel = 0;

    for (auto& results : results_per_point) {
        if (results.iter == max_iterations) {
            // points inside the Mandelbrot Set are always painted black
            image_data[static_cast<std::size_t>(pixel)] = PixelColor{0, 0, 0};
        } else {
            // The equalized iteration value (in the range of 0 .. max_iterations) represents the
            // position of the pixel color in the color gradiant and needs to be mapped to 0.0 .. 1.0.
            // To achieve smooth coloring we need to edge the equalized iteration towards the next
            // iteration, determined by the distance between the two iterations.
            const auto iter_curr = equalized_iterations[static_cast<std::size_t>(results.iter)];
            const auto iter_next = equalized_iterations[static_cast<std::size_t>(results.iter + 1)];

            const auto smoothed_iteration = std::lerp(iter_curr, iter_next, results.distance_to_next_iteration);
            const auto pos_in_gradient = smoothed_iteration / static_cast<float>(max_iterations);

            color_from_gradient(gradient, pos_in_gradient, image_data[static_cast<std::size_t>(pixel)]);
        }

        ++pixel;
    }
}

void save_image(const std::string& filename, const std::vector<PixelColor>& image_data)
{
    std::ofstream out{filename, std::ofstream::binary};

    if (!out.is_open())
        throw std::runtime_error("unable to open output file");

    for (auto p : image_data)
        out << p.r << p.g << p.b;
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
        std::cout << durations[0] << " s" << '\n';
    } else {
        std::vector<double> sorted_values{durations};
        std::sort(sorted_values.begin(), sorted_values.end());

        std::cout << "mean: " << mean(sorted_values) << " s, median: " << median(sorted_values) << " s (repetitions=" << sorted_values.size() << ") ["
                  << sorted_values[0];

        for (auto p = std::next(sorted_values.cbegin()); p != sorted_values.cend(); ++p)
            std::cout << ", " << *p;

        std::cout << "]" << '\n';
    }
}

template <typename T>
T eval_arg(const char* s, T min, T max)
{
    T value;
    std::istringstream in{s};

    if (!(in >> value) || value < min || value > max)
        throw std::runtime_error("invalid value");

    return value;
}

auto eval_args(const int argc, char const* argv[])
{
    if (argc != 10)
        throw std::runtime_error("invalid number of arguments");

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

auto go(const int image_width, const int image_height, const int max_iterations, const double center_x, const double center_y, const double height, const Gradient& gradient, const int repetitions) noexcept
{
    // iterations_histogram: for simplicity we only use indices [1] .. [max_iterations], [0] is unused
    std::vector<int> iterations_histogram(static_cast<std::size_t>(max_iterations + 1));

    // For every point store a tuple consisting of the final iteration and (for escaped points)
    // the distance to the next iteration (as value of 0.0 .. 1.0).
    std::vector<CalculationResult> results_per_point(static_cast<std::size_t>(image_width * image_height));

    std::vector<PixelColor> image_data(static_cast<std::size_t>(image_width * image_height));
    std::vector<double> durations;

    for (int i = 0; i < repetitions; ++i) {
        const auto t1 = std::chrono::high_resolution_clock::now();
        mandelbrot_calc(image_width, image_height, max_iterations, center_x, center_y, height, iterations_histogram, results_per_point);
        mandelbrot_colorize(max_iterations, gradient, image_data, iterations_histogram, results_per_point);
        const auto t2 = std::chrono::high_resolution_clock::now();

        durations.push_back(std::chrono::duration<double>{t2 - t1}.count());
    }

    return std::make_tuple(image_data, durations);
}

int main(int argc, char const* argv[])
{
    auto [image_width, image_height, max_iterations, repetitions, center_x, center_y, height, gradient_filename, filename] = eval_args(argc, argv);
    auto gradient = load_gradient(gradient_filename);

    auto [image_data, durations] = go(image_width, image_height, max_iterations, center_x, center_y, height, gradient, repetitions);

    save_image(filename, image_data);
    show_summary(durations);
}
