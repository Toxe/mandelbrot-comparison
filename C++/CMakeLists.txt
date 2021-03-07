cmake_minimum_required(VERSION 3.14)
project("Mandelbrot Comparison (C++)" LANGUAGES CXX)

if(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    string(REGEX REPLACE " /W[0-4]" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
endif()

add_executable(mandelbrot_cpp mandelbrot.cpp)
target_compile_features(mandelbrot_cpp PUBLIC cxx_std_20)
set_target_properties(mandelbrot_cpp PROPERTIES CXX_EXTENSIONS OFF)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(mandelbrot_cpp PRIVATE -Wall -Wextra -Wmost -pedantic -Wconversion -Wfloat-equal -Wold-style-cast -fcolor-diagnostics)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    target_compile_options(mandelbrot_cpp PRIVATE /W4)
endif()
