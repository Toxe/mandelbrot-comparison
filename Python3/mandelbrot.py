# coding=utf8

# Python3 Mandelbrot.
#
# Usage: python3 mandelbrot.py <image_width> <image_height> <max_iterations> <repetitions (1+)> <center x> <center y> <section height> <gradient filename> <output filename>
# Example: python3 mandelbrot.py 320 200 20 1 -0.5 0.0 2.0 blue.gradient mandelbrot.raw
#
# Gradient file example:
#   0.0: 0.0, 0.0, 0.0
#   0.5: 0.0, 0.0, 1.0
#   1.0: 1.0, 1.0, 1.0
#
# Tobias Brückner, 2021


import re
import sys
import math
import statistics

from time import time
from math import log, sqrt


REGEXP_GRADIENT_LINE = re.compile(r'([0-9]*\.?[0-9]+):\s*([0-9]*\.?[0-9]+),\s*([0-9]*\.?[0-9]+),\s*([0-9]*\.?[0-9]+)')


class GradientColor:
    def __init__(self, pos, r, g, b):
        self.pos = float(pos)
        self.update_colors(r, g, b)

    def update_colors(self, r, g, b):
        self.r = float(r)
        self.g = float(g)
        self.b = float(b)

    def __lt__(self, other):
        return self.pos < other.pos


class Gradient:
    def __init__(self):
        self.colors = []


# Compare two float values for "enough" equality.
float_epsilon = sys.float_info.epsilon

def equal_enough(a, b):
    a = abs(a)
    b = abs(b)
    return abs(a - b) <= max(a, b) * float_epsilon


def gradient_get_color_at_position(gradient, pos):
    return next(filter(lambda col: equal_enough(col.pos, pos), gradient.colors), None)


def load_gradient(filename):
    gradient = Gradient()
    gradient.colors.append(GradientColor(0.0, 0.0, 0.0, 0.0))
    gradient.colors.append(GradientColor(1.0, 1.0, 1.0, 1.0))

    with open(filename, "r") as f:
        for line in f:
            match = REGEXP_GRADIENT_LINE.match(line)
            if match:
                data = match.groups()
                col = gradient_get_color_at_position(gradient, float(data[0]))
                if col:
                    col.update_colors(data[1], data[2], data[3])
                else:
                    gradient.colors.append(GradientColor(data[0], data[1], data[2], data[3]))
    gradient.colors.sort()
    return gradient


def color_from_gradient_range(left_color, right_color, pos):
    pos2 = (pos - left_color.pos) / (right_color.pos - left_color.pos)
    r = (right_color.r - left_color.r) * pos2 + left_color.r
    g = (right_color.g - left_color.g) * pos2 + left_color.g
    b = (right_color.b - left_color.b) * pos2 + left_color.b
    return r, g, b


def color_from_gradient(gradient, pos):
    for colors in zip(gradient.colors[:-1], gradient.colors[1:]):
        if pos >= colors[0].pos and pos <= colors[1].pos:
            return color_from_gradient_range(colors[0], colors[1], pos)
    return None


def mandelbrot_calc(image_width, image_height, max_iterations, center_x, center_y, height, histogram, iterations_per_pixel, smoothed_distances_to_next_iteration_per_pixel):
    width = height * (float(image_width) / float(image_height))
    x_left = center_x - width / 2.0
    # x_right = center_x + width / 2.0
    y_top = center_y + height / 2.0
    # y_bottom = center_y - height / 2.0

    bailout = 20.0
    bailout_squared = bailout * bailout
    log_log_bailout = log(log(bailout))
    log_2 = log(2.0)

    histogram[:] = [0] * len(histogram)

    for pixel_y in range(image_height):
        y0 = y_top - height * (float(pixel_y) / float(image_height))

        for pixel_x in range(image_width):
            x0 = x_left + width * (float(pixel_x) / float(image_width))

            x = 0.0
            y = 0.0

            x_squared = 0.0
            y_squared = 0.0

            # iteration, will be from 1 to max_iterations once the loop is done
            iter = 0

            for iter in range(max_iterations + 1):
                x_squared = x*x
                y_squared = y*y

                if x_squared + y_squared >= bailout_squared:
                    break

                y = 2.0*x*y + y0
                x = x_squared - y_squared + x0

            if iter < max_iterations:
                final_magnitude = sqrt(x_squared + y_squared)
                smoothed_distances_to_next_iteration_per_pixel[pixel_y * image_width + pixel_x] = 1.0 - min(1.0, (log(log(final_magnitude)) - log_log_bailout) / log_2)
                histogram[iter] += 1  # no need to count histogram[max_iterations]

            iterations_per_pixel[pixel_y * image_width + pixel_x] = iter  # 1 .. max_iterations


def mandelbrot_colorize(image_width, image_height, max_iterations, gradient, image_data, histogram, iterations_per_pixel, smoothed_distances_to_next_iteration_per_pixel, normalized_colors):
    # Sum all iterations, not counting the last one at position histogram[max_iterations] (which
    # are points in the Mandelbrot Set).
    total_iterations = float(sum(histogram[1:-1]))

    # Normalize the colors (0.0 .. 1.0) based on how often they are used in the image, not counting
    # histogram[max_iterations] (which are points in the Mandelbrot Set).
    running_total = 0

    for i in range(1, max_iterations):
        running_total += histogram[i]
        normalized_colors[i] = float(running_total) / total_iterations

    for pixel in range(image_width * image_height):
        iter = iterations_per_pixel[pixel]  # 1 .. max_iterations

        if iter == max_iterations:
            # pixels with max. iterations (aka. inside the Mandelbrot Set) are always black
            image_data[3 * pixel + 0] = 0
            image_data[3 * pixel + 1] = 0
            image_data[3 * pixel + 2] = 0
        else:
            # we use the color of the previous iteration in order to cover the full gradient range
            color_of_previous_iter = normalized_colors[iter - 1]
            color_of_current_iter  = normalized_colors[iter]
            smoothed_distance_to_next_iteration = smoothed_distances_to_next_iteration_per_pixel[pixel]  # 0 .. <1.0
            pos_in_gradient = color_of_previous_iter + smoothed_distance_to_next_iteration * (color_of_current_iter - color_of_previous_iter)

            r, g, b = color_from_gradient(gradient, pos_in_gradient)

            image_data[3 * pixel + 0] = int(255.0 * r)
            image_data[3 * pixel + 1] = int(255.0 * g)
            image_data[3 * pixel + 2] = int(255.0 * b)


def save_image(filename, image_data):
    with open(filename, "wb") as f:
        f.write(image_data)


def show_summary(durations):
    if len(durations) == 1:
        print(f"{durations[0]} s")
    else:
        print(f"mean: {statistics.mean(durations)} s, median: {statistics.median(durations)} s (repetitions={len(durations)})", sorted(durations))


def eval_int_arg(s, min, max):
    value = int(s)
    if value < min or value > max:
        raise ValueError(f"invalid value {s}")
    return value


def eval_float_arg(s, min, max):
    value = float(s)
    if value < min or value > max or math.isnan(value) or math.isinf(value):
        raise ValueError(f"invalid value {s}")
    return value


def eval_args():
    if len(sys.argv) < 10:
        raise RuntimeError("invalid number of arguments")

    image_width    = eval_int_arg(sys.argv[1], 1, 100000)
    image_height   = eval_int_arg(sys.argv[2], 1, 100000)
    max_iterations = eval_int_arg(sys.argv[3], 1, 1000000000)
    repetitions    = eval_int_arg(sys.argv[4], 1, 1000000)
    center_x       = eval_float_arg(sys.argv[5], -100.0, 100.0)
    center_y       = eval_float_arg(sys.argv[6], -100.0, 100.0)
    height         = eval_float_arg(sys.argv[7], -100.0, 100.0)
    colors         = sys.argv[8]
    filename       = sys.argv[9]

    return image_width, image_height, max_iterations, repetitions, center_x, center_y, height, colors, filename


def go(image_width, image_height, max_iterations, center_x, center_y, height, gradient, image_data, durations, repetitions):
    # histogram & normalized_colors: for simplicity we only use indices [1] .. [max_iterations], [0] is unused
    histogram = [0] * (max_iterations + 1)
    normalized_colors = [0.0] * (max_iterations + 1)
    iterations_per_pixel = [0] * (image_width * image_height)
    smoothed_distances_to_next_iteration_per_pixel = [0.0] * (image_width * image_height)

    for _ in range(repetitions):
        t1 = time()
        mandelbrot_calc(image_width, image_height, max_iterations, center_x, center_y, height, histogram, iterations_per_pixel, smoothed_distances_to_next_iteration_per_pixel)
        mandelbrot_colorize(image_width, image_height, max_iterations, gradient, image_data, histogram, iterations_per_pixel, smoothed_distances_to_next_iteration_per_pixel, normalized_colors)
        t2 = time()
        durations.append(t2 - t1)


def main():
    image_width, image_height, max_iterations, repetitions, center_x, center_y, height, gradient_filename, filename = eval_args()
    gradient = load_gradient(gradient_filename)
    image_data = bytearray([0] * (image_width * image_height * 3))
    durations = []

    go(image_width, image_height, max_iterations, center_x, center_y, height, gradient, image_data, durations, repetitions)

    save_image(filename, image_data)
    show_summary(durations)


if __name__ == '__main__':
    main()
