/*
 * C Mandelbrot.
 *
 * Usage: ./mandelbrot <image_width> <image_height> <iterations> <repetitions (1+)> <center x> <center y> <section height> <output filename>
 * Example: ./mandelbrot 320 200 20 1 -0.5 0.0 2.0 mandelbrot.raw
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>

void mandelbrot(int image_width, int image_height, int max_iter, double center_x, double center_y, double height, unsigned char *image_data)
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

    for (pixel_y = 0; pixel_y < image_height; ++pixel_y) {
        for (pixel_x = 0; pixel_x < image_width; ++pixel_x) {
            x0 = x_left + width * ((double) pixel_x / (double) image_width);
            y0 = y_top - height * ((double) pixel_y / (double) image_height);

            x = 0;
            y = 0;

            iter = 0;

            while (iter < max_iter) {
                x_squared = x*x;
                y_squared = y*y;

                if (x_squared + y_squared >= 4.0)
                    break;

                xtemp = x_squared - y_squared + x0;
                y = 2*x*y + y0;
                x = xtemp;

                ++iter;
            }

            if (iter == max_iter) {
                image_data[3 * (pixel_y * image_width + pixel_x) + 0] = 0;
                image_data[3 * (pixel_y * image_width + pixel_x) + 1] = 0;
                image_data[3 * (pixel_y * image_width + pixel_x) + 2] = 0;
            } else {
                image_data[3 * (pixel_y * image_width + pixel_x) + 0] = 255;
                image_data[3 * (pixel_y * image_width + pixel_x) + 1] = 255;
                image_data[3 * (pixel_y * image_width + pixel_x) + 2] = 255;
            }
        }
    }
}

void die(int error)
{
    fprintf(stderr, "Error: %d\n", error);
    exit(error);
}

void save_image(const unsigned char *image_data, int width, int height)
{
    FILE *fp;

    if (!(fp = fopen("mandelbrot.raw", "w")))
        die(5);

    fwrite(image_data, sizeof(unsigned char), 3 * width * height, fp);
    fclose(fp);
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
        die(4);

    memcpy(sorted_values, values, num_values * sizeof(double));
    qsort(sorted_values, num_values, sizeof(double), cmp_doubles_func);

    if (num_values % 2 == 1)
        return sorted_values[(num_values - 1) / 2];

    double d = (sorted_values[num_values/2 - 1] + sorted_values[num_values/2]) / 2;

    free(sorted_values);
    return d;
}

void print_durations(const double *values, int num_values)
{
    double *sorted_values;

    if (!(sorted_values = malloc(num_values * sizeof(double))))
        die(3);

    memcpy(sorted_values, values, num_values * sizeof(double));
    qsort(sorted_values, num_values, sizeof(double), cmp_doubles_func);

    printf("[");

    for (int i = 0; i < num_values; ++i)
        printf("%s%f", (i > 0) ? ", " : "", sorted_values[i]);

    printf("]");

    free(sorted_values);
}

void show_summary(const double *durations, int repetitions)
{
    if (repetitions == 1) {
        printf("%f s\n", durations[0]);
    } else {
        printf("mean: %f s, median: %f s (repetitions=%d) ", mean(durations, repetitions), median(durations, repetitions), repetitions);
        print_durations(durations, repetitions);
        printf("\n");
    }
}

double gettime()
{
    struct timeval tv;

    if (gettimeofday(&tv, NULL) != 0)
        die(2);

    return (double) tv.tv_sec + (double) tv.tv_usec / 1000000.0;
}

void eval_args(int argc, char **argv, int *image_width, int *image_height, int *iterations, int *repetitions,
               double *center_x, double *center_y, double *height, char **filename)
{
    if (argc < 9)
        die(1);

    *image_width  = atoi(argv[1]);
    *image_height = atoi(argv[2]);
    *iterations   = atoi(argv[3]);
    *repetitions  = atoi(argv[4]);
    *center_x     = atof(argv[5]);
    *center_y     = atof(argv[6]);
    *height       = atof(argv[7]);
    *filename     = argv[8];
}

int main(int argc, char **argv)
{
    int image_width, image_height;
    int iterations, repetitions;
    double center_x, center_y, height;
    double *durations;
    char *filename;
    unsigned char *image_data;

    eval_args(argc, argv, &image_width, &image_height, &iterations, &repetitions, &center_x, &center_y, &height, &filename);

    image_data = malloc(image_width * image_height * 3 * sizeof(unsigned char));
    durations = malloc(repetitions * sizeof(double));

    for (int i = 0; i < repetitions; ++i) {
        double t1, t2;

        t1 = gettime();
        mandelbrot(image_width, image_height, iterations, center_x, center_y, height, image_data);
        t2 = gettime();

        durations[i] = t2 - t1;
    }

    save_image(image_data, image_width, image_height);
    show_summary(durations, repetitions);

    free(image_data);
    free(durations);

    return 0;
}
