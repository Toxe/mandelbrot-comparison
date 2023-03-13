/*
 * C Mandelbrot.
 *
 * Usage: ./build/C/mandelbrot_c <image width> <image height> <max iterations> <repetitions (1+)> <center x> <center y> <section height> <gradient filename> <output filename>
 * Example: ./build/C/mandelbrot_c 800 600 200 1 -0.8 0.0 2.2 gradients/benchmark.gradient mandelbrot.raw
 *
 * Tobias Brückner, 2016
 */

#if _WIN32 || _WIN64
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <errno.h>
#include <limits.h>

#if _WIN32 || _WIN64
#define WIN32_LEAN_AND_MEAN
#include <time.h>
#include <Windows.h>
#else
#include <sys/time.h>
#endif

typedef struct
{
    unsigned char r, g, b;
} pixel_color_t;

typedef struct
{
    float pos;
    float r, g, b;
} gradient_color_t;

typedef struct
{
    gradient_color_t *colors;
    int num_colors;
} gradient_t;

typedef struct
{
    int iter;
    float distance_to_next_iteration;
} calculation_result_t;


void die(char *reason)
{
    fprintf(stderr, "Runtime error: %s\n", reason);
    exit(1);
}

float lerp_float(const float a, const float b, const float t)
{
    return (1.0f - t) * a + t * b;
}

double lerp_double(const double a, const double b, const double t)
{
    return (1.0 - t) * a + t * b;
}

void cumsum(const int *iterations_histogram, int *cdf, int size)
{
    int total = 0;

    for (int i = 0; i < size; ++i) {
        total += iterations_histogram[i];
        cdf[i] = total;
    }
}

// Compare two float values for "enough" equality.
int equal_enough_float(float a, float b)
{
    a = fabsf(a);
    b = fabsf(b);

    return fabsf(a - b) <= fmaxf(a, b) * FLT_EPSILON;
}

// Compare two double values for "enough" equality.
int equal_enough_double(double a, double b)
{
    a = fabs(a);
    b = fabs(b);

    return fabs(a - b) <= fmax(a, b) * DBL_EPSILON;
}

gradient_color_t *gradient_get_color_at_position(gradient_t *gradient, float pos)
{
    for (int i = 0; i < gradient->num_colors; ++i)
        if (equal_enough_float(gradient->colors[i].pos, pos))
            return &gradient->colors[i];

    return NULL;
}

int cmp_color_pos_func(const void *p1, const void *p2)
{
    const gradient_color_t *a = (const gradient_color_t *) p1;
    const gradient_color_t *b = (const gradient_color_t *) p2;

    if (equal_enough_float(a->pos, b->pos))
        return 0;

    return (a->pos < b->pos) ? -1 : 1;
}

gradient_t *load_gradient(char *filename)
{
    FILE *fp;
    gradient_t *gradient;
    char buf[256];

    if (!(gradient = (gradient_t *) malloc(sizeof(gradient_t))))
        die("alloc memory");

    gradient->num_colors = 2;

    if (!(gradient->colors = (gradient_color_t *) malloc((size_t) gradient->num_colors * sizeof(gradient_color_t))))
        die("alloc memory");

    gradient->colors[0].pos = 0.0;
    gradient->colors[0].r = 0.0;
    gradient->colors[0].g = 0.0;
    gradient->colors[0].b = 0.0;

    gradient->colors[1].pos = 1.0;
    gradient->colors[1].r = 1.0;
    gradient->colors[1].g = 1.0;
    gradient->colors[1].b = 1.0;

    if (!(fp = fopen(filename, "r")))
        return NULL;

    while (fgets(buf, sizeof(buf), fp)) {
        float pos, r, g, b;
        gradient_color_t *col;

        if (sscanf(buf, "%f: %f, %f, %f", &pos, &r, &g, &b) != 4)
            continue;

        if (!(col = gradient_get_color_at_position(gradient, pos))) {
            gradient->num_colors++;
            gradient->colors = (gradient_color_t *) realloc(gradient->colors, (size_t) gradient->num_colors * sizeof(gradient_color_t));
            col = &gradient->colors[gradient->num_colors - 1];
        }

        col->pos = pos;
        col->r = r;
        col->g = g;
        col->b = b;
    }

    qsort(gradient->colors, (size_t) gradient->num_colors, sizeof(gradient_color_t), cmp_color_pos_func);

    fclose(fp);
    return gradient;
}

void free_gradient(gradient_t *gradient)
{
    free(gradient->colors);
    free(gradient);
}

void color_from_gradient_range(const gradient_color_t *left, const gradient_color_t *right, float pos, pixel_color_t *color)
{
    const float relative_pos_between_colors = (pos - left->pos) / (right->pos - left->pos);
    color->r = (unsigned char) (255.0f * lerp_float(left->r, right->r, relative_pos_between_colors));
    color->g = (unsigned char) (255.0f * lerp_float(left->g, right->g, relative_pos_between_colors));
    color->b = (unsigned char) (255.0f * lerp_float(left->b, right->b, relative_pos_between_colors));
}

int color_from_gradient(const gradient_t *gradient, float pos, pixel_color_t *color)
{
    gradient_color_t *left = &gradient->colors[0];

    for (int i = 1; i < gradient->num_colors; ++i) {
        gradient_color_t *right = &gradient->colors[i];

        if (pos >= left->pos && pos <= right->pos) {
            color_from_gradient_range(left, right, pos, color);
            return 0;
        }

        left = right;
    }

    return -1;
}

void mandelbrot_calc(int image_width, int image_height, int max_iterations, double center_x, double center_y, double height,
                     int *iterations_histogram, calculation_result_t *results_per_point)
{
    const double width = height * ((double) image_width / (double) image_height);

    const double x_left   = center_x - width / 2.0;
    const double x_right  = center_x + width / 2.0;
    const double y_top    = center_y + height / 2.0;
    const double y_bottom = center_y - height / 2.0;

    const double bailout = 20.0;
    const double bailout_squared = bailout * bailout;
    const double log_log_bailout = log(log(bailout));
    const double log_2 = log(2.0);

    double final_magnitude = 0.0;

    memset(iterations_histogram, 0, (size_t) (max_iterations + 1) * sizeof(int));

    int pixel = 0;

    for (int pixel_y = 0; pixel_y < image_height; ++pixel_y) {
        const double y0 = lerp_double(y_top, y_bottom, (double) pixel_y / (double) image_height);

        for (int pixel_x = 0; pixel_x < image_width; ++pixel_x) {
            const double x0 = lerp_double(x_left, x_right, (double) pixel_x / (double) image_width);

            double x = 0.0;
            double y = 0.0;

            // iteration, will be from 1 .. max_iterations once the loop is done
            int iter = 0;

            while (iter < max_iterations) {
                const double x_squared = x * x;
                const double y_squared = y * y;

                if (x_squared + y_squared >= bailout_squared) {
                    final_magnitude = sqrt(x_squared + y_squared);
                    break;
                }

                y = 2.0 * x * y + y0;
                x = x_squared - y_squared + x0;

                ++iter;
            }

            if (iter < max_iterations) {
                ++iterations_histogram[iter]; // iter: 1 .. max_iterations-1, no need to count iterations_histogram[max_iterations]
                results_per_point[pixel].iter = iter;
                results_per_point[pixel].distance_to_next_iteration = 1.0f - fminf(1.0f, (float) ((log(log(final_magnitude)) - log_log_bailout) / log_2));
            } else {
                results_per_point[pixel].iter = iter;
                results_per_point[pixel].distance_to_next_iteration = 0.0;
            }

            ++pixel;
        }
    }
}

float *equalize_histogram(const int *iterations_histogram, const int max_iterations)
{
    const int size = max_iterations + 1;

    int *cdf;
    float *equalized_iterations;

    if (!(cdf = (int *) malloc((size_t) size * sizeof(int))))
        die("alloc memory");

    if (!(equalized_iterations = (float *) malloc((size_t) size * sizeof(float))))
        die("alloc memory");

    // Calculate the CDF (Cumulative Distribution Function) by accumulating all iteration counts.
    // Element [0] is unused and iterations_histogram[max_iterations] should be zero (as we do not count
    // the iterations of the points inside the Mandelbrot Set).
    cumsum(iterations_histogram, cdf, size);

    // Get the minimum value in the CDF that is bigger than zero and the sum of all iteration counts
    // from iterations_histogram (which is the last value of the CDF).
    int cdf_min = 0;

    for (int i = 0; i < size; ++i) {
        if (cdf[i] > 0) {
            cdf_min = cdf[i];
            break;
        }
    }

    const int total_iterations = cdf[size - 1];

    // normalize all values from the CDF that are bigger than zero to a range of 0.0 .. max_iterations
    const float f = (float) max_iterations / (float) (total_iterations - cdf_min);

    for (int i = 0; i < size; ++i)
        equalized_iterations[i] = cdf[i] > 0 ? f * (float) (cdf[i] - cdf_min) : 0.0f;

    free(cdf);

    return equalized_iterations;
}

void mandelbrot_colorize(const int image_width, const int image_height, const int max_iterations, const gradient_t *gradient,
                         pixel_color_t *image_data, const int *iterations_histogram, const calculation_result_t *results_per_point)
{
    float *equalized_iterations = equalize_histogram(iterations_histogram, max_iterations);

    for (int pixel = 0; pixel < image_width * image_height; ++pixel) {
        const calculation_result_t *results = &results_per_point[pixel];

        if (results->iter == max_iterations) {
            // points inside the Mandelbrot Set are always painted black
            image_data[pixel].r = 0;
            image_data[pixel].g = 0;
            image_data[pixel].b = 0;
        } else {
            // The equalized iteration value (in the range of 0 .. max_iterations) represents the
            // position of the pixel color in the color gradiant and needs to be mapped to 0.0 .. 1.0.
            // To achieve smooth coloring we need to edge the equalized iteration towards the next
            // iteration, determined by the distance between the two iterations.
            const float iter_curr = equalized_iterations[results->iter];
            const float iter_next = equalized_iterations[results->iter + 1];

            const float smoothed_iteration = lerp_float(iter_curr, iter_next, results->distance_to_next_iteration);
            const float pos_in_gradient = smoothed_iteration / (float) max_iterations;

            color_from_gradient(gradient, pos_in_gradient, &image_data[pixel]);
        }
    }

    free(equalized_iterations);
}

int save_image(const char *filename, const pixel_color_t *image_data, int width, int height)
{
    FILE *fp;

    if (!(fp = fopen(filename, "wb")))
        return -1;

    fwrite(image_data, sizeof(pixel_color_t), (size_t) (width * height), fp);
    fclose(fp);

    return 0;
}

int cmp_doubles_func(const void *p1, const void *p2)
{
    double a = *((const double *) p1);
    double b = *((const double *) p2);

    if (equal_enough_double(a, b))
        return 0;

    return (a < b) ? -1 : 1;
}

double mean(const double *values, int num_values)
{
    double sum = 0.0;

    for (int i = 0; i < num_values; ++i)
        sum += values[i];

    return sum / num_values;
}

double median(const double *values, int num_values)
{
    double *sorted_values;

    if (!(sorted_values = (double *) malloc((size_t) num_values * sizeof(double))))
        die("alloc memory");

    memcpy(sorted_values, values, (size_t) num_values * sizeof(double));
    qsort(sorted_values, (size_t) num_values, sizeof(double), cmp_doubles_func);

    double d;

    if (num_values % 2 == 1)
        d = sorted_values[(num_values - 1) / 2];
    else
        d = (sorted_values[num_values/2 - 1] + sorted_values[num_values/2]) / 2.0;

    free(sorted_values);
    return d;
}

char *durations2string(const double *values, int num_values)
{
    double *sorted_values;
    char *buf;
    char tmp[32];
    int buf_len = 256;
    int tmp_len;
    int pos;

    if (!(sorted_values = (double *) malloc((size_t) num_values * sizeof(double))))
        die("alloc memory");

    memcpy(sorted_values, values, (size_t) num_values * sizeof(double));
    qsort(sorted_values, (size_t) num_values, sizeof(double), cmp_doubles_func);

    if (!(buf = (char *) malloc((size_t) buf_len * sizeof(char))))
        die("alloc memory");

    sprintf(buf, "[");
    pos = (int) strlen(buf);

    for (int i = 0; i < num_values; ++i) {
        sprintf(tmp, "%s%f", (i > 0) ? ", " : "", sorted_values[i]);
        tmp_len = (int) strlen(tmp);

        if ((pos + tmp_len + 1) > buf_len) {
            buf_len += 256;

            if (!(buf = (char *) realloc(buf, (size_t) buf_len)))
                die("alloc memory");
        }

        memcpy(buf + pos, tmp, (size_t) tmp_len + 1);
        pos += tmp_len;
    }

    sprintf(buf + pos, "]");
    free(sorted_values);

    return buf;
}

void show_summary(const double *durations, int repetitions)
{
    if (repetitions == 1) {
        printf("%f s\n", durations[0]);
    } else {
        char *s = durations2string(durations, repetitions);
        printf("mean: %f s, median: %f s (repetitions=%d) %s\n", mean(durations, repetitions), median(durations, repetitions), repetitions, s);
        free(s);
    }
}

double gettime()
{
#if _WIN32 || _WIN64
    LARGE_INTEGER t, freq;

    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&t);

    t.QuadPart *= 1000000;
    t.QuadPart /= freq.QuadPart;

    return (double) t.QuadPart / 1000000.0;
#else
    struct timeval tv;

    if (gettimeofday(&tv, NULL) != 0)
        die("gettimeofday");

    return (double) tv.tv_sec + (double) tv.tv_usec / 1000000.0;
#endif
}

int eval_int_arg(const char *s, int min, int max)
{
    int value;

    errno = 0;
    value = strtol(s, NULL, 0);

    if (errno == ERANGE || value < min || value > max)
        die("invalid value");

    return value;
}

double eval_double_arg(const char *s, double min, double max)
{
    double value;

    errno = 0;
    value = strtod(s, NULL);

    if (errno == ERANGE || value < min || value > max || isnan(value) || isinf(value))
        die("invalid value");

    return value;
}

void eval_args(int argc, char **argv, int *image_width, int *image_height, int *max_iterations, int *repetitions,
               double *center_x, double *center_y, double *height, char **colors, char **filename)
{
    if (argc != 10)
        die("invalid number of arguments");

    *image_width    = eval_int_arg(argv[1], 1, 100000);
    *image_height   = eval_int_arg(argv[2], 1, 100000);
    *max_iterations = eval_int_arg(argv[3], 1, 1000000000);
    *repetitions    = eval_int_arg(argv[4], 1, 1000000);
    *center_x       = eval_double_arg(argv[5], -100.0, 100.0);
    *center_y       = eval_double_arg(argv[6], -100.0, 100.0);
    *height         = eval_double_arg(argv[7], -100.0, 100.0);
    *colors         = argv[8];
    *filename       = argv[9];
}

void go(int image_width, int image_height, int max_iterations, double center_x, double center_y, double height, gradient_t *gradient, pixel_color_t *image_data, double *durations, int repetitions)
{
    int *iterations_histogram;
    calculation_result_t *results_per_point;

    // iterations_histogram: for simplicity we only use indices [1] .. [max_iterations], [0] is unused
    if (!(iterations_histogram = (int *) malloc((size_t) (max_iterations + 1) * sizeof(int))))
        die("alloc memory");

    // For every point store a tuple consisting of the final iteration and (for escaped points)
    // the distance to the next iteration (as value of 0.0 .. 1.0).
    if (!(results_per_point = (calculation_result_t *) malloc((size_t) (image_width * image_height) * sizeof(calculation_result_t))))
        die("alloc memory");

    for (int i = 0; i < repetitions; ++i) {
        double t1, t2;

        t1 = gettime();
        mandelbrot_calc(image_width, image_height, max_iterations, center_x, center_y, height, iterations_histogram, results_per_point);
        mandelbrot_colorize(image_width, image_height, max_iterations, gradient, image_data, iterations_histogram, results_per_point);
        t2 = gettime();

        durations[i] = t2 - t1;
    }

    free(results_per_point);
    free(iterations_histogram);
}

int main(int argc, char **argv)
{
    int image_width, image_height;
    int max_iterations, repetitions;
    double center_x, center_y, height;
    double *durations;
    char *gradient_filename, *filename;
    pixel_color_t *image_data;
    gradient_t *gradient;

    eval_args(argc, argv, &image_width, &image_height, &max_iterations, &repetitions, &center_x, &center_y, &height, &gradient_filename, &filename);

    if (!(gradient = load_gradient(gradient_filename)))
        die("unable to load gradient");

    if (!(image_data = (pixel_color_t *) malloc((size_t) (image_width * image_height) * sizeof(pixel_color_t))))
        die("alloc memory");

    if (!(durations = (double *) malloc((size_t) repetitions * sizeof(double))))
        die("alloc memory");

    go(image_width, image_height, max_iterations, center_x, center_y, height, gradient, image_data, durations, repetitions);

    if (save_image(filename, image_data, image_width, image_height) < 0)
        die("unable to save output file");

    show_summary(durations, repetitions);

    free(image_data);
    free(durations);
    free_gradient(gradient);

    return 0;
}
