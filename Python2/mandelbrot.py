# coding=utf8

# Python2 Mandelbrot.
#
# Usage: python2.7 mandelbrot.py <image_width> <image_height> <max_iterations> <repetitions (1+)> <center x> <center y> <section height> <gradient filename> <output filename>
# Example: python2.7 mandelbrot.py 320 200 20 1 -0.5 0.0 2.0 blue.gradient mandelbrot.raw
#
# Gradient file example:
#   0.0: 0.0, 0.0, 0.0
#   0.5: 0.0, 0.0, 1.0
#   1.0: 1.0, 1.0, 1.0
#
# Tobias Brückner, 2016

from __future__ import print_function

import re
import sys
import math

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

    def __cmp__(self, other):
        if equal_enough(self.pos, other.pos):
            return 0
        elif self.pos < other.pos:
            return -1
        else:
            return 1


class Gradient:
    def __init__(self):
        self.colors = []


def lerp(a, b, t):
    return (1.0 - t) * a + t * b


def cumsum(hist):
    cdf = [0] * len(hist)
    total = 0
    for i, n in enumerate(hist):
        total += n
        cdf[i] = total
    return cdf


# Compare two float values for "enough" equality.
float_epsilon = sys.float_info.epsilon

def equal_enough(a, b):
    a = abs(a)
    b = abs(b)
    return abs(a - b) <= max(a, b) * float_epsilon


def gradient_get_color_at_position(gradient, pos):
    return next(iter(filter(lambda col: equal_enough(col.pos, pos), gradient.colors)), None)


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
    relative_pos_between_colors = (pos - left_color.pos) / (right_color.pos - left_color.pos)
    r = lerp(left_color.r, right_color.r, relative_pos_between_colors)
    g = lerp(left_color.g, right_color.g, relative_pos_between_colors)
    b = lerp(left_color.b, right_color.b, relative_pos_between_colors)
    return r, g, b


def color_from_gradient(gradient, pos):
    for colors in zip(gradient.colors[:-1], gradient.colors[1:]):
        if pos >= colors[0].pos and pos <= colors[1].pos:
            return color_from_gradient_range(colors[0], colors[1], pos)
    return None


def mandelbrot_calc(image_width, image_height, max_iterations, center_x, center_y, height, iterations_per_pixel, distances_to_next_iteration_per_pixel):
    width = height * image_width / image_height
    x_left = center_x - width / 2.0
    # x_right = center_x + width / 2.0
    y_top = center_y + height / 2.0
    # y_bottom = center_y - height / 2.0

    bailout = 20.0
    bailout_squared = bailout * bailout
    log_log_bailout = log(log(bailout))
    log_2 = log(2.0)

    # for simplicity we only use indices [1] .. [max_iterations], [0] is unused
    iterations_histogram = [0] * (max_iterations + 1)

    pixel = 0

    for pixel_y in xrange(image_height):
        y0 = y_top - height * pixel_y / image_height

        for pixel_x in xrange(image_width):
            x0 = x_left + width * pixel_x / image_width

            x = 0.0
            y = 0.0

            x_squared = 0.0
            y_squared = 0.0

            # iteration, will be from 1 .. max_iterations once the loop is done
            iter = 0

            for iter in xrange(max_iterations + 1):
                x_squared = x * x
                y_squared = y * y

                if x_squared + y_squared >= bailout_squared:
                    break

                y = 2.0 * x * y + y0
                x = x_squared - y_squared + x0

            if iter < max_iterations:
                final_magnitude = sqrt(x_squared + y_squared)
                distances_to_next_iteration_per_pixel[pixel] = 1.0 - min(1.0, (log(log(final_magnitude)) - log_log_bailout) / log_2)
                iterations_histogram[iter] += 1  # no need to count histogram[max_iterations]

            iterations_per_pixel[pixel] = iter
            pixel += 1
    return iterations_histogram


def equalize_histogram(iterations_histogram, max_iterations):
    # Calculate the CDF (Cumulative Distribution Function) by accumulating all iteration counts.
    # Element [0] is unused and iterations_histogram[max_iterations] should be zero (as we do not count
    # the iterations of the points inside the Mandelbrot Set).
    cdf = cumsum(iterations_histogram)

    # Get the minimum value in the CDF that is bigger than zero and the sum of all iteration counts
    # from iterations_histogram (which is the last value of the CDF).
    cdf_min = next(iter(filter(lambda x: x > 0, cdf)))
    total_iterations = cdf[-1]

    # normalize all values from the CDF that are bigger than zero to a range of 0.0 .. max_iterations
    f = max_iterations / float(total_iterations - cdf_min)
    return list(map(lambda c: (c - cdf_min) * f if c > 0 else 0, cdf))


def mandelbrot_colorize(image_width, image_height, max_iterations, gradient, image_data, iterations_histogram, iterations_per_pixel, distances_to_next_iteration_per_pixel):
    equalized_iterations = equalize_histogram(iterations_histogram, max_iterations)

    for pixel in xrange(image_width * image_height):
        iter = iterations_per_pixel[pixel]  # 1 .. max_iterations

        if iter == max_iterations:
            # points inside the Mandelbrot Set are always painted black
            image_data[3 * pixel + 0] = 0
            image_data[3 * pixel + 1] = 0
            image_data[3 * pixel + 2] = 0
        else:
            # The equalized iteration value (in the range of 0 .. max_iterations) represents the
            # position of the pixel color in the color gradiant and needs to be mapped to 0.0 .. 1.0.
            # To achieve smooth coloring we need to edge the equalized iteration towards the next
            # iteration, determined by the distance between the two iterations.
            iter_curr = equalized_iterations[iter]
            iter_next = equalized_iterations[iter + 1]

            smoothed_iteration = lerp(iter_curr, iter_next, distances_to_next_iteration_per_pixel[pixel])
            pos_in_gradient = smoothed_iteration / max_iterations

            r, g, b = color_from_gradient(gradient, pos_in_gradient)

            image_data[3 * pixel + 0] = int(255.0 * r)
            image_data[3 * pixel + 1] = int(255.0 * g)
            image_data[3 * pixel + 2] = int(255.0 * b)


def save_image(filename, image_data):
    with open(filename, "wb") as f:
        f.write(image_data)


def mean(values):
    return math.fsum(values) / len(values)


def median(values):
    values = sorted(values)
    count = len(values)
    if count % 2 == 1:
        return values[int((count-1) / 2)]
    else:
        return (values[int(count/2 - 1)] + values[int(count/2)]) / 2.0


def show_summary(durations):
    if len(durations) == 1:
        print("%f s" % durations[0])
    else:
        print("mean: %f s, median: %f s (repetitions=%d) %s" % (mean(durations), median(durations), len(durations), str(sorted(durations))))


def eval_int_arg(s, min, max):
    value = int(s)
    if value < min or value > max:
        raise ValueError("invalid value %s" % s)
    return value


def eval_float_arg(s, min, max):
    value = float(s)
    if value < min or value > max or math.isnan(value) or math.isinf(value):
        raise ValueError("invalid value %s" % s)
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
    iterations_per_pixel = [0] * (image_width * image_height)
    distances_to_next_iteration_per_pixel = [0.0] * (image_width * image_height)

    for _ in xrange(repetitions):
        t1 = time()
        iterations_histogram = mandelbrot_calc(image_width, image_height, max_iterations, center_x, center_y, height, iterations_per_pixel, distances_to_next_iteration_per_pixel)
        mandelbrot_colorize(image_width, image_height, max_iterations, gradient, image_data, iterations_histogram, iterations_per_pixel, distances_to_next_iteration_per_pixel)
        t2 = time()
        durations.append(t2 - t1)


def main():
    image_width, image_height, max_iterations, repetitions, center_x, center_y, height, gradient_filename, filename = eval_args()
    gradient = load_gradient(gradient_filename)
    image_data = bytearray(image_width * image_height * 3)
    durations = []

    go(image_width, image_height, max_iterations, center_x, center_y, height, gradient, image_data, durations, repetitions)

    save_image(filename, image_data)
    show_summary(durations)


if __name__ == "__main__":
    main()
