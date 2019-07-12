/*
 * C Mandelbrot.
 *
 * Usage: mandelbrot_c <image_width> <image_height> <max_iterations> <repetitions (1+)> <center x> <center y> <section height> <gradient filename> <output filename>
 * Example: mandelbrot_c 320 200 20 1 -0.5 0.0 2.0 blue.gradient mandelbrot.raw
 *
 * Gradient file example:
 *   0.0: 0.0, 0.0, 0.0
 *   0.5: 0.0, 0.0, 1.0
 *   1.0: 1.0, 1.0, 1.0
 *
 * Tobias Brückner, 2016
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <errno.h>
#include <limits.h>
#include <sys/time.h>

typedef enum exitcode_t {
    ERROR_ALLOC_MEMORY = 1,
    ERROR_EVAL_ARGS,
    ERROR_LOAD_GRADIENT,
    ERROR_SAVE_IMAGE,
    ERROR_GETTIME
} exitcode_t;

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


void die(exitcode_t error)
{
    fprintf(stderr, "Error: %d\n", error);
    exit(error);
}

// Compare two float values for "enough" equality.
int equal_enough(float a, float b)
{
    a = fabsf(a);
    b = fabsf(b);

    return fabsf(a - b) <= fmaxf(a, b) * FLT_EPSILON;
}

gradient_color_t *gradient_get_color_at_position(gradient_t *gradient, float pos)
{
    for (int i = 0; i < gradient->num_colors; ++i)
        if (equal_enough(gradient->colors[i].pos, pos))
            return &gradient->colors[i];

    return NULL;
}

int cmp_color_pos_func(const void *p1, const void *p2)
{
    const gradient_color_t *a = (const gradient_color_t *) p1;
    const gradient_color_t *b = (const gradient_color_t *) p2;

    if (equal_enough(a->pos, b->pos))
        return 0;

    return (a->pos < b->pos) ? -1 : 1;
}

gradient_t *load_gradient(char *filename)
{
    FILE *fp;
    gradient_t *gradient;
    char buf[256];

    if (!(gradient = malloc(sizeof(gradient_t))))
        die(ERROR_ALLOC_MEMORY);

    gradient->num_colors = 2;

    if (!(gradient->colors = malloc(gradient->num_colors * sizeof(gradient_color_t))))
        die(ERROR_ALLOC_MEMORY);

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
            gradient->colors = realloc(gradient->colors, gradient->num_colors * sizeof(gradient_color_t));
            col = &gradient->colors[gradient->num_colors - 1];
        }

        col->pos = pos;
        col->r = r;
        col->g = g;
        col->b = b;
    }

    qsort(gradient->colors, gradient->num_colors, sizeof(gradient_color_t), cmp_color_pos_func);

    fclose(fp);
    return gradient;
}

void free_gradient(gradient_t *gradient)
{
    free(gradient->colors);
    free(gradient);
}

void color_from_gradient_range(gradient_color_t *left, gradient_color_t *right, float pos, float *r, float *g, float *b)
{
    float pos2 = (pos - left->pos) / (right->pos - left->pos);

    *r = ((right->r - left->r) * pos2) + left->r;
    *g = ((right->g - left->g) * pos2) + left->g;
    *b = ((right->b - left->b) * pos2) + left->b;
}

int color_from_gradient(gradient_t *gradient, float pos, float *r, float *g, float *b)
{
    gradient_color_t *left, *right;

    left = &gradient->colors[0];

    for (int i = 1; i < gradient->num_colors; ++i) {
        right = &gradient->colors[i];

        if (pos >= left->pos && pos <= right->pos) {
            color_from_gradient_range(left, right, pos, r, g, b);
            return 0;
        }

        left = right;
    }

    return -1;
}

void mandelbrot_calc(int image_width, int image_height, int max_iterations, double center_x, double center_y, double height,
                     int *histogram, int *iterations_per_pixel, float *smoothed_distances_to_next_iteration_per_pixel)
{
    const double width = height * ((double) image_width / (double) image_height);

    const double x_left   = center_x - width / 2.0;
 // const double x_right  = center_x + width / 2.0;
    const double y_top    = center_y + height / 2.0;
 // const double y_bottom = center_y - height / 2.0;

    const double bailout = 20.0;
    const double bailout_squared = bailout * bailout;
    const double log_log_bailout = log(log(bailout));
    const double log_2 = log(2.0);

    double final_magnitude = 0.0;

    memset(histogram, 0, (max_iterations + 1) * sizeof(int));

    for (int pixel_y = 0; pixel_y < image_height; ++pixel_y) {
        const double y0 = y_top - height * ((double) pixel_y / (double) image_height);

        for (int pixel_x = 0; pixel_x < image_width; ++pixel_x) {
            const double x0 = x_left + width * ((double) pixel_x / (double) image_width);

            double x = 0.0;
            double y = 0.0;

            // iteration, will be from 1 to max_iterations once the loop is done
            int iter = 0;

            while (iter < max_iterations) {
                const double x_squared = x*x;
                const double y_squared = y*y;

                if (x_squared + y_squared >= bailout_squared) {
                    final_magnitude = sqrt(x_squared + y_squared);
                    break;
                }

                const double xtemp = x_squared - y_squared + x0;
                y = 2.0*x*y + y0;
                x = xtemp;

                ++iter;
            }

            const int pixel = pixel_y * image_width + pixel_x;

            if (iter < max_iterations) {
                smoothed_distances_to_next_iteration_per_pixel[pixel] = 1.0 - fmin(1.0, (log(log(final_magnitude)) - log_log_bailout) / log_2);
                ++histogram[iter];  // no need to count histogram[max_iterations]
            }

            iterations_per_pixel[pixel] = iter;  // 1 .. max_iterations
        }
    }
}

void mandelbrot_colorize(int image_width, int image_height, int max_iterations, gradient_t *gradient,
                         unsigned char *image_data, int *histogram, int *iterations_per_pixel, float *smoothed_distances_to_next_iteration_per_pixel, float *normalized_colors)
{
    // Sum all iterations, not counting the last one at position histogram[max_iterations] (which
    // are points in the Mandelbrot Set).
    int total_iterations = 0;

    for (int i = 1; i < max_iterations; ++i)
        total_iterations += histogram[i];

    // Normalize the colors (0.0 .. 1.0) based on how often they are used in the image, not counting
    // histogram[max_iterations] (which are points in the Mandelbrot Set).
    int running_total = 0;

    for (int i = 1; i < max_iterations; ++i) {
        running_total += histogram[i];
        normalized_colors[i] = (double) running_total / (double) total_iterations;
    }

    for (int pixel_y = 0; pixel_y < image_height; ++pixel_y) {
        for (int pixel_x = 0; pixel_x < image_width; ++pixel_x) {
            int pixel = pixel_y * image_width + pixel_x;
            int iter = iterations_per_pixel[pixel];  // 1 .. max_iterations

            if (iter == max_iterations) {
                // pixels with max. iterations (aka. inside the Mandelbrot Set) are always black
                image_data[3 * pixel + 0] = 0;
                image_data[3 * pixel + 1] = 0;
                image_data[3 * pixel + 2] = 0;
            } else {
                // we use the color of the previous iteration in order to cover the full gradient range
                float r, g, b;
                float color_of_previous_iter = normalized_colors[iter - 1];
                float color_of_current_iter  = normalized_colors[iter];

                float smoothed_distance_to_next_iteration = smoothed_distances_to_next_iteration_per_pixel[pixel];  // 0 .. <1.0
                float pos_in_gradient = color_of_previous_iter + smoothed_distance_to_next_iteration * (color_of_current_iter - color_of_previous_iter);

                color_from_gradient(gradient, pos_in_gradient, &r, &g, &b);

                image_data[3 * pixel + 0] = (unsigned char) (255.0f * r);
                image_data[3 * pixel + 1] = (unsigned char) (255.0f * g);
                image_data[3 * pixel + 2] = (unsigned char) (255.0f * b);
            }
        }
    }
}

int save_image(const char *filename, const unsigned char *image_data, int width, int height)
{
    FILE *fp;

    if (!(fp = fopen(filename, "w")))
        return -1;

    fwrite(image_data, sizeof(unsigned char), 3 * width * height, fp);
    fclose(fp);

    return 0;
}

int cmp_doubles_func(const void *p1, const void *p2)
{
    double a = *((const double *) p1);
    double b = *((const double *) p2);

    if (equal_enough(a, b))
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

    if (!(sorted_values = malloc(num_values * sizeof(double))))
        die(ERROR_ALLOC_MEMORY);

    memcpy(sorted_values, values, num_values * sizeof(double));
    qsort(sorted_values, num_values, sizeof(double), cmp_doubles_func);

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

    if (!(sorted_values = malloc(num_values * sizeof(double))))
        die(ERROR_ALLOC_MEMORY);

    memcpy(sorted_values, values, num_values * sizeof(double));
    qsort(sorted_values, num_values, sizeof(double), cmp_doubles_func);

    if (!(buf = malloc(buf_len * sizeof(char))))
        die(ERROR_ALLOC_MEMORY);

    sprintf(buf, "[");
    pos = strlen(buf);

    for (int i = 0; i < num_values; ++i) {
        sprintf(tmp, "%s%f", (i > 0) ? ", " : "", sorted_values[i]);
        tmp_len = strlen(tmp);

        if ((pos + tmp_len + 1) > buf_len) {
            buf_len += 256;

            if (!(buf = realloc(buf, buf_len)))
                die(ERROR_ALLOC_MEMORY);
        }

        memcpy(buf + pos, tmp, tmp_len + 1);
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
    struct timeval tv;

    if (gettimeofday(&tv, NULL) != 0)
        die(ERROR_GETTIME);

    return (double) tv.tv_sec + (double) tv.tv_usec / 1000000.0;
}

int eval_int_arg(const char *s, int min, int max)
{
    int value;

    errno = 0;
    value = strtol(s, NULL, 0);

    if (errno == ERANGE || value < min || value > max)
        die(ERROR_EVAL_ARGS);

    return value;
}

double eval_double_arg(const char *s, double min, double max)
{
    double value;

    errno = 0;
    value = strtod(s, NULL);

    if (errno == ERANGE || value < min || value > max || isnan(value) || isinf(value))
        die(ERROR_EVAL_ARGS);

    return value;
}

int eval_args(int argc, char **argv, int *image_width, int *image_height, int *max_iterations, int *repetitions,
              double *center_x, double *center_y, double *height, char **colors, char **filename)
{
    if (argc < 10)
        return -1;

    *image_width    = eval_int_arg(argv[1], 1, 100000);
    *image_height   = eval_int_arg(argv[2], 1, 100000);
    *max_iterations = eval_int_arg(argv[3], 1, 1000000000);
    *repetitions    = eval_int_arg(argv[4], 1, 1000000);
    *center_x       = eval_double_arg(argv[5], -100.0, 100.0);
    *center_y       = eval_double_arg(argv[6], -100.0, 100.0);
    *height         = eval_double_arg(argv[7], -100.0, 100.0);
    *colors         = argv[8];
    *filename       = argv[9];

    return 0;
}

void go(int image_width, int image_height, int max_iterations, double center_x, double center_y, double height, gradient_t *gradient, unsigned char *image_data, double *durations, int repetitions)
{
    int *histogram;
    int *iterations_per_pixel;
    float *smoothed_distances_to_next_iteration_per_pixel;
    float *normalized_colors;

    // histogram & normalized_colors: for simplicity we only use indices [1] .. [max_iterations], [0] is unused
    if (!(histogram = malloc((max_iterations + 1) * sizeof(int))))
        die(ERROR_ALLOC_MEMORY);

    if (!(normalized_colors = malloc((max_iterations + 1) * sizeof(float))))
        die(ERROR_ALLOC_MEMORY);

    if (!(iterations_per_pixel = malloc(image_width * image_height * sizeof(int))))
        die(ERROR_ALLOC_MEMORY);

    if (!(smoothed_distances_to_next_iteration_per_pixel = malloc(image_width * image_height * sizeof(float))))
        die(ERROR_ALLOC_MEMORY);

    for (int i = 0; i < repetitions; ++i) {
        double t1, t2;

        t1 = gettime();
        mandelbrot_calc(image_width, image_height, max_iterations, center_x, center_y, height, histogram, iterations_per_pixel, smoothed_distances_to_next_iteration_per_pixel);
        mandelbrot_colorize(image_width, image_height, max_iterations, gradient, image_data, histogram, iterations_per_pixel, smoothed_distances_to_next_iteration_per_pixel, normalized_colors);
        t2 = gettime();

        durations[i] = t2 - t1;
    }

    free(iterations_per_pixel);
    free(smoothed_distances_to_next_iteration_per_pixel);
    free(histogram);
    free(normalized_colors);
}

int main(int argc, char **argv)
{
    int image_width, image_height;
    int max_iterations, repetitions;
    double center_x, center_y, height;
    double *durations;
    char *gradient_filename, *filename;
    unsigned char *image_data;
    gradient_t *gradient;

    if (eval_args(argc, argv, &image_width, &image_height, &max_iterations, &repetitions, &center_x, &center_y, &height, &gradient_filename, &filename) < 0)
        die(ERROR_EVAL_ARGS);

    if (!(gradient = load_gradient(gradient_filename)))
        die(ERROR_LOAD_GRADIENT);

    if (!(image_data = malloc(image_width * image_height * 3 * sizeof(unsigned char))))
        die(ERROR_ALLOC_MEMORY);

    if (!(durations = malloc(repetitions * sizeof(double))))
        die(ERROR_ALLOC_MEMORY);

    go(image_width, image_height, max_iterations, center_x, center_y, height, gradient, image_data, durations, repetitions);

    if (save_image(filename, image_data, image_width, image_height) < 0)
        die(ERROR_SAVE_IMAGE);

    show_summary(durations, repetitions);

    free(image_data);
    free(durations);
    free_gradient(gradient);

    return 0;
}
