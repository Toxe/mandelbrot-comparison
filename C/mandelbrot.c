/*
 * C Mandelbrot.
 *
 * Usage: ./mandelbrot <image_width> <image_height> <iterations> <repetitions (1+)> <center x> <center y> <section height> <gradient filename> <output filename>
 * Example: ./mandelbrot 320 200 20 1 -0.5 0.0 2.0 blue.gradient mandelbrot.raw
 *
 * Gradient file example:
 *   0.0: 0.0, 0.0, 0.0
 *   0.5: 0.0, 0.0, 1.0
 *   1.0: 1.0, 1.0, 1.0
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>

typedef enum error_t {
    ERROR_ALLOC_MEMORY = 1,
    ERROR_EVAL_ARGS,
    ERROR_LOAD_GRADIENT,
    ERROR_SAVE_IMAGE,
    ERROR_GETTIME
} error_t;

typedef struct
{
    double pos;
    double r, g, b;
} gradient_color_t;

typedef struct
{
    gradient_color_t *colors;

    int num_colors;
} gradient_t;


void die(error_t error)
{
    fprintf(stderr, "Error: %d\n", error);
    exit(error);
}

gradient_color_t *gradient_get_color_at_position(gradient_t *gradient, double pos)
{
    for (int i = 0; i < gradient->num_colors; ++i)
        if (gradient->colors[i].pos == pos)
            return &gradient->colors[i];

    return NULL;
}

int cmp_color_pos_func(const void *p1, const void *p2)
{
    const gradient_color_t *a = (const gradient_color_t *) p1;
    const gradient_color_t *b = (const gradient_color_t *) p2;

    if (a->pos == b->pos)
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
        double pos, r, g, b;
        gradient_color_t *col;

        if (sscanf(buf, "%lf: %lf, %lf, %lf", &pos, &r, &g, &b) != 4)
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

void color_from_gradient_range(gradient_color_t *left, gradient_color_t *right, double pos, double *r, double *g, double *b)
{
    if (left->pos == pos || left->pos == right->pos) {
        *r = left->r;
        *g = left->g;
        *b = left->b;
        return;
    }

    if (right->pos == pos) {
        *r = right->r;
        *g = right->g;
        *b = right->b;
        return;
    }

    double d = pos - left->pos;
    double pos2 = d / (right->pos - left->pos);

    *r = ((right->r - left->r) * pos2) + left->r;
    *g = ((right->g - left->g) * pos2) + left->g;
    *b = ((right->b - left->b) * pos2) + left->b;
}

int color_from_gradient(gradient_t *gradient, double pos, double *r, double *g, double *b)
{
    if (pos < 0.0)
        pos = 0.0;

    if (pos > 1.0)
        pos = 1.0;

    gradient_color_t *left, *right;

    left = &gradient->colors[0];

    for (int i = 0; i < gradient->num_colors; ++i) {
        right = &gradient->colors[i];

        if (pos >= left->pos && pos <= right->pos) {
            color_from_gradient_range(left, right, pos, r, g, b);
            return 0;
        }

        left = right;
    }

    return -1;
}

void mandelbrot(int image_width, int image_height, int max_iterations, double center_x, double center_y, double height, gradient_t *gradient, unsigned char *image_data)
{
    double width = height * ((double) image_width / (double) image_height);

    double x_left   = center_x - width / 2.0;
 // double x_right  = center_x + width / 2.0;
    double y_top    = center_y + height / 2.0;
 // double y_bottom = center_y - height / 2.0;

    int pixel_x, pixel_y;
    int iter;
    double x0, y0;
    double x, y;
    double xtemp;
    double x_squared, y_squared;

    int *histogram;
    int *image_iterations_per_pixel;
    double *normalized_colors;

    if (!(histogram = calloc(max_iterations, sizeof(int))))
        die(ERROR_ALLOC_MEMORY);

    if (!(image_iterations_per_pixel = malloc(image_width * image_height * sizeof(int))))
        die(ERROR_ALLOC_MEMORY);

    if (!(normalized_colors = calloc(max_iterations, sizeof(double))))
        die(ERROR_ALLOC_MEMORY);

    for (pixel_y = 0; pixel_y < image_height; ++pixel_y) {
        for (pixel_x = 0; pixel_x < image_width; ++pixel_x) {
            x0 = x_left + width * ((double) pixel_x / (double) image_width);
            y0 = y_top - height * ((double) pixel_y / (double) image_height);

            x = 0;
            y = 0;

            iter = 0;

            while (iter < max_iterations) {
                x_squared = x*x;
                y_squared = y*y;

                if (x_squared + y_squared >= 4.0)
                    break;

                xtemp = x_squared - y_squared + x0;
                y = 2*x*y + y0;
                x = xtemp;

                ++iter;
            }

            // count iterations (0 .. max - 1) in histogram and store iterations per pixel
            ++histogram[iter - 1];
            image_iterations_per_pixel[pixel_y * image_width + pixel_x] = iter - 1;
        }
    }

    // sum all iterations (not counting the last one (max. iteration - 1))
    int total_iterations = 0;

    for (int i = 0; i < max_iterations - 1; ++i)
        total_iterations += histogram[i];

    // normalize the colors (0.0 .. 1.0) based on how often they are used in the image (not counting max. iteration)
    int running_total = 0;

    for (int i = 0; i < max_iterations - 1; ++i) {
        running_total += histogram[i];
        normalized_colors[i] = (double) running_total / (double) total_iterations;
    }

    for (pixel_y = 0; pixel_y < image_height; ++pixel_y) {
        for (pixel_x = 0; pixel_x < image_width; ++pixel_x) {
            int iter = image_iterations_per_pixel[pixel_y * image_width + pixel_x];
            double r, g, b;

            if (iter == (max_iterations - 1)) {
                // pixels with max. iterations are always black
                r = 0.0;
                g = 0.0;
                b = 0.0;
            } else {
                color_from_gradient(gradient, normalized_colors[iter], &r, &g, &b);
            }

            image_data[3 * (pixel_y * image_width + pixel_x) + 0] = (unsigned char) (255.0 * r);
            image_data[3 * (pixel_y * image_width + pixel_x) + 1] = (unsigned char) (255.0 * g);
            image_data[3 * (pixel_y * image_width + pixel_x) + 2] = (unsigned char) (255.0 * b);
        }
    }

    free(image_iterations_per_pixel);
    free(histogram);
    free(normalized_colors);
}

int save_image(const unsigned char *image_data, int width, int height)
{
    FILE *fp;

    if (!(fp = fopen("mandelbrot.raw", "w")))
        return -1;

    fwrite(image_data, sizeof(unsigned char), 3 * width * height, fp);
    fclose(fp);

    return 0;
}

int cmp_doubles_func(const void *p1, const void *p2)
{
    double a = *((const double *) p1);
    double b = *((const double *) p2);

    if (a == b)
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

    if (num_values % 2 == 1)
        return sorted_values[(num_values - 1) / 2];

    double d = (sorted_values[num_values/2 - 1] + sorted_values[num_values/2]) / 2;

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

int eval_args(int argc, char **argv, int *image_width, int *image_height, int *iterations, int *repetitions,
               double *center_x, double *center_y, double *height, char **colors, char **filename)
{
    if (argc < 10)
        return -1;

    *image_width  = atoi(argv[1]);
    *image_height = atoi(argv[2]);
    *iterations   = atoi(argv[3]);
    *repetitions  = atoi(argv[4]);
    *center_x     = atof(argv[5]);
    *center_y     = atof(argv[6]);
    *height       = atof(argv[7]);
    *colors       = argv[8];
    *filename     = argv[9];

    return 0;
}

int main(int argc, char **argv)
{
    int image_width, image_height;
    int iterations, repetitions;
    double center_x, center_y, height;
    double *durations;
    char *gradient_filename, *filename;
    unsigned char *image_data;
    gradient_t *gradient;

    if (eval_args(argc, argv, &image_width, &image_height, &iterations, &repetitions, &center_x, &center_y, &height, &gradient_filename, &filename) < 0)
        die(ERROR_EVAL_ARGS);

    if (!(gradient = load_gradient(gradient_filename)))
        die(ERROR_LOAD_GRADIENT);

    if (!(image_data = malloc(image_width * image_height * 3 * sizeof(unsigned char))))
        die(ERROR_ALLOC_MEMORY);

    if (!(durations = malloc(repetitions * sizeof(double))))
        die(ERROR_ALLOC_MEMORY);

    for (int i = 0; i < repetitions; ++i) {
        double t1, t2;

        t1 = gettime();
        mandelbrot(image_width, image_height, iterations, center_x, center_y, height, gradient, image_data);
        t2 = gettime();

        durations[i] = t2 - t1;
    }

    if (save_image(image_data, image_width, image_height) < 0)
        die(ERROR_SAVE_IMAGE);

    show_summary(durations, repetitions);

    free(image_data);
    free(durations);
    free_gradient(gradient);

    return 0;
}
