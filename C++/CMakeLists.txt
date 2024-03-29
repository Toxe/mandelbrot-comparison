cmake_minimum_required(VERSION 3.14)
project("Mandelbrot Comparison (C++)" LANGUAGES CXX)

find_package(Threads REQUIRED)

set(ALL_TARGETS
    mandelbrot_cpp
    mandelbrot_cpp_async
)

add_executable(mandelbrot_cpp mandelbrot.cpp)
add_executable(mandelbrot_cpp_async mandelbrot_async.cpp)

foreach(target ${ALL_TARGETS})
    set_target_properties(${target} PROPERTIES CXX_EXTENSIONS OFF)
    target_compile_features(${target} PUBLIC cxx_std_20)
    target_compile_options(${target} PRIVATE ${SANITIZER_FLAGS} ${DEFAULT_COMPILER_OPTIONS} ${DEFAULT_COMPILER_WARNINGS})
    target_link_libraries(${target} PRIVATE ${SANITIZER_FLAGS})
endforeach()

target_link_libraries(mandelbrot_cpp_async PRIVATE Threads::Threads)
