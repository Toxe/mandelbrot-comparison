# Define default compiler options for Clang, GCC and MSVC.

# Clang specific options
set(DEFAULT_CLANG_OPTIONS
        -fcolor-diagnostics
)

# GCC specific options
set(DEFAULT_GCC_OPTIONS
        -fdiagnostics-color
)

# MSVC specific options
set(DEFAULT_MSVC_OPTIONS
        /permissive-
        /utf-8
)

# Clang-cl specific options
set(DEFAULT_CLANG_CL_OPTIONS
        ${DEFAULT_MSVC_OPTIONS}
        /EHsc  # enable exception handling
)

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND NOT CMAKE_CXX_COMPILER_FRONTEND_VARIANT STREQUAL "MSVC")
    # Clang
    set(DEFAULT_COMPILER_OPTIONS ${DEFAULT_CLANG_OPTIONS})
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND CMAKE_CXX_COMPILER_FRONTEND_VARIANT STREQUAL "MSVC")
    # Clang-cl
    set(DEFAULT_COMPILER_OPTIONS ${DEFAULT_CLANG_CL_OPTIONS})
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    # GCC
    set(DEFAULT_COMPILER_OPTIONS ${DEFAULT_GCC_OPTIONS})
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    # MSVC
    set(DEFAULT_COMPILER_OPTIONS ${DEFAULT_MSVC_OPTIONS})
endif ()
