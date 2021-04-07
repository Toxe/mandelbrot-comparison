# Define default compiler options and warnings for Clang, GCC and MSVC.

# shared options between Clang and GCC
set(COMMON_CLANG_AND_GCC_WARNINGS
    -Wall
    -Wextra
    -pedantic
    -Wconversion
    -Wsign-conversion
    -Wfloat-equal
    -Wold-style-cast
    -Wshadow
    -Wunused
)

# Clang specific options
set(DEFAULT_CLANG_OPTIONS
    ${COMMON_CLANG_AND_GCC_WARNINGS}
    -Wmost
    -fcolor-diagnostics
)

# GCC specific options
set(DEFAULT_GCC_OPTIONS
    ${COMMON_CLANG_AND_GCC_WARNINGS}
    -fdiagnostics-color
)

# MSVC specific options
set(DEFAULT_MSVC_OPTIONS
    /W4
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(DEFAULT_COMPILER_OPTIONS_AND_WARNINGS ${DEFAULT_CLANG_OPTIONS})
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(DEFAULT_COMPILER_OPTIONS_AND_WARNINGS ${DEFAULT_GCC_OPTIONS})
elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    string(REGEX REPLACE " /W[0-4]" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    set(DEFAULT_COMPILER_OPTIONS_AND_WARNINGS ${DEFAULT_MSVC_OPTIONS})
endif()
