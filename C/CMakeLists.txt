cmake_minimum_required(VERSION 3.14)
project("Mandelbrot Comparison (C)" LANGUAGES C)

if(CMAKE_C_COMPILER_ID MATCHES "MSVC")
    string(REGEX REPLACE " /W[0-4]" "" CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
endif()

add_executable(mandelbrot_c mandelbrot.c)
target_compile_features(mandelbrot_c PUBLIC c_std_11)
set_target_properties(mandelbrot_c PROPERTIES C_EXTENSIONS OFF)

if(CMAKE_C_COMPILER_ID MATCHES "Clang")
    target_compile_options(mandelbrot_c PRIVATE -Wall -Wextra -Wmost -pedantic -Wfloat-equal -fcolor-diagnostics)
elseif(CMAKE_C_COMPILER_ID MATCHES "MSVC")
    target_compile_options(mandelbrot_c PRIVATE /W4 /wd4996 /wd4706 /wd4267)
endif()

SET(NEEDS_LIBM Clang GNU)

if(CMAKE_C_COMPILER_ID IN_LIST NEEDS_LIBM)
    target_link_libraries(mandelbrot_c PRIVATE m)
endif()
