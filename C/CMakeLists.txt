cmake_minimum_required(VERSION 3.14)
project("Mandelbrot Comparison (C)" LANGUAGES C)

add_executable(mandelbrot_c mandelbrot.c)
target_compile_features(mandelbrot_c PUBLIC c_std_11)
target_compile_options(mandelbrot_c PRIVATE -Wall -Wextra -Wmost -pedantic -Wfloat-equal)
set_target_properties(mandelbrot_c PROPERTIES C_EXTENSIONS OFF)
