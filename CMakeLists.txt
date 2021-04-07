cmake_minimum_required(VERSION 3.14)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

project("Mandelbrot Comparison" LANGUAGES C CXX)

# default compiler options and warnings
include(cmake/DefaultCompilerOptionsAndWarnings.cmake)

# static analyzers
include(cmake/StaticAnalyzers.cmake)

# sanitizers
include(cmake/Sanitizers.cmake)

add_subdirectory(C)
add_subdirectory(C++)

if(CMAKE_HOST_APPLE)
    add_subdirectory(Swift)
endif()
