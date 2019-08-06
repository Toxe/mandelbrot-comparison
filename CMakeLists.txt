cmake_minimum_required(VERSION 3.14)
project("Mandelbrot Comparison" LANGUAGES C CXX)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

add_subdirectory(C)
add_subdirectory(C++)

if(CMAKE_HOST_APPLE)
    add_subdirectory(Swift)
endif()
