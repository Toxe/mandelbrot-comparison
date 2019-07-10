add_executable(mandelbrot_cpp mandelbrot.cpp)
target_compile_features(mandelbrot_cpp PUBLIC cxx_std_17)
target_compile_options(mandelbrot_cpp PRIVATE -Wall -Wextra -Wmost -pedantic -Wconversion -Wfloat-equal -Wold-style-cast)
set_target_properties(mandelbrot_cpp PROPERTIES CXX_EXTENSIONS OFF)
