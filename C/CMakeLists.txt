cmake_minimum_required(VERSION 3.14)
project("Mandelbrot Comparison (C)" LANGUAGES C)

if(CMAKE_C_COMPILER_ID MATCHES "MSVC")
    string(REGEX REPLACE " /W[0-4]" "" CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
endif()

add_executable(mandelbrot_c mandelbrot.c)
set_target_properties(mandelbrot_c PROPERTIES C_EXTENSIONS OFF)
target_compile_features(mandelbrot_c PUBLIC c_std_11)
target_compile_options(mandelbrot_c PRIVATE ${SANITIZER_FLAGS} ${DEFAULT_COMPILER_OPTIONS} ${DEFAULT_COMPILER_WARNINGS})
target_link_libraries(mandelbrot_c PRIVATE ${SANITIZER_FLAGS})

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND NOT CMAKE_CXX_COMPILER_FRONTEND_VARIANT STREQUAL "MSVC")
    # Clang
    target_link_libraries(mandelbrot_c PRIVATE m)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    # GCC
    target_link_libraries(mandelbrot_c PRIVATE m)
endif ()
