# coding=utf8

# Python Mandelbrot.
#
# Usage: python mandelbrot.py <image_width> <image_height> <max_iterations> <repetitions (1+)> <center x> <center y> <section height> <gradient filename> <output filename>
# Example: python mandelbrot.py 320 200 20 1 -0.5 0.0 2.0 blue.gradient mandelbrot.raw
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


class ExitCode:
    ALLOC_MEMORY = 1
    EVAL_ARGS = 2
    LOAD_GRADIENT = 3
    SAVE_IMAGE = 4
    GETTIME = 5


class GradientColor:
    def __init__(self, pos, r, g, b):
        self.pos = pos
        self.r = r
        self.g = g
        self.b = b
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


# Compare two float values for "enough" equality.
float_epsilon = sys.float_info.epsilon

def equal_enough(a, b):
    a = abs(a)
    b = abs(b)
    return abs(a - b) <= max(a, b) * float_epsilon


def gradient_get_color_at_position(gradient, pos):
    for col in gradient.colors:
        if equal_enough(col.pos, pos):
            return col
    return None


def load_gradient(filename):
    try:
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
                        col.r = float(data[1])
                        col.g = float(data[2])
                        col.b = float(data[3])
                    else:
                        gradient.colors.append(GradientColor(float(data[0]), float(data[1]), float(data[2]), float(data[3])))
        gradient.colors.sort()
        return gradient
    except IOError as e:
        return None
    return None


def color_from_gradient_range(left_color, right_color, pos):
    pos2 = (pos - left_color.pos) / (right_color.pos - left_color.pos)
    r = (right_color.r - left_color.r) * pos2 + left_color.r
    g = (right_color.g - left_color.g) * pos2 + left_color.g
    b = (right_color.b - left_color.b) * pos2 + left_color.b
    return r, g, b


def color_from_gradient(gradient, pos):
    if pos < 0.0:
        pos = 0.0
    if pos > 1.0:
        pos = 1.0

    left_color = gradient.colors[0]

    for right_color in gradient.colors[1:]:
        if pos >= left_color.pos and pos <= right_color.pos:
            return color_from_gradient_range(left_color, right_color, pos)
        left_color = right_color

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

    for i in xrange(max_iterations + 1):
        histogram[i] = 0

    for pixel_y in xrange(image_height):
        for pixel_x in xrange(image_width):
            x0 = x_left + width * (float(pixel_x) / float(image_width))
            y0 = y_top - height * (float(pixel_y) / float(image_height))

            x = 0.0
            y = 0.0

            # iteration, will be from 1 to max_iterations once the loop is done
            iter = 0

            while iter < max_iterations:
                x_squared = x*x
                y_squared = y*y

                if x_squared + y_squared >= bailout_squared:
                    break

                xtemp = x_squared - y_squared + x0
                y = 2.0*x*y + y0
                x = xtemp
                iter += 1

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

    for i in xrange(1, max_iterations):
        running_total += histogram[i]
        normalized_colors[i] = float(running_total) / total_iterations

    for pixel_y in xrange(image_height):
        for pixel_x in xrange(image_width):
            pixel = pixel_y * image_width + pixel_x
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
    try:
        with open(filename, "wb") as f:
            f.write(image_data)
    except IOError as e:
        return False
    return True


def mean(values):
    return math.fsum(values) / len(values)


def median(values):
    values = sorted(values)
    count = len(values)
    if count % 2 == 1:
        return values[(count-1) / 2]
    else:
        return (values[count/2 - 1] + values[count/2]) / 2.0


def show_summary(durations):
    if len(durations) == 1:
        print("%f s" % durations[0])
    else:
        print("mean: %f s, median: %f s (repetitions=%d) %s" % (mean(durations), median(durations), len(durations), str(sorted(durations))))


def die(error):
    print("Error: %d" % error, file=sys.stderr)
    sys.exit(error)


def eval_int_arg(s, min, max):
    try:
        value = int(s)
        if value < min or value > max:
            die(ExitCode.EVAL_ARGS)
        return value
    except ValueError as e:
        die(ExitCode.EVAL_ARGS)


def eval_float_arg(s, min, max):
    try:
        value = float(s)
        if value < min or value > max or math.isnan(value) or math.isinf(value):
            die(ExitCode.EVAL_ARGS)
        return value
    except ValueError as e:
        die(ExitCode.EVAL_ARGS)


def eval_args():
    if len(sys.argv) < 10:
        die(ExitCode.EVAL_ARGS)

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

    for i in xrange(repetitions):
        t1 = time()
        mandelbrot_calc(image_width, image_height, max_iterations, center_x, center_y, height, histogram, iterations_per_pixel, smoothed_distances_to_next_iteration_per_pixel)
        mandelbrot_colorize(image_width, image_height, max_iterations, gradient, image_data, histogram, iterations_per_pixel, smoothed_distances_to_next_iteration_per_pixel, normalized_colors)
        t2 = time()
        durations.append(t2 - t1)


if __name__ == '__main__':
    image_width, image_height, max_iterations, repetitions, center_x, center_y, height, gradient_filename, filename = eval_args()
    gradient = load_gradient(gradient_filename)

    if gradient == None:
        die(ExitCode.LOAD_GRADIENT)

    image_data = bytearray([0] * (image_width * image_height * 3))
    durations = []

    go(image_width, image_height, max_iterations, center_x, center_y, height, gradient, image_data, durations, repetitions)

    if not save_image(filename, image_data):
        die(ExitCode.SAVE_IMAGE)
    show_summary(durations)
