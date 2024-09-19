# Define default compiler warnings for Clang, GCC and MSVC.

# https://clang.llvm.org/docs/DiagnosticsReference.html
# https://gcc.gnu.org/onlinedocs/gcc/Warning-Options.html
# https://docs.microsoft.com/en-us/cpp/preprocessor/compiler-warnings-that-are-off-by-default?view=msvc-170

# shared options between Clang, Clang-cl and GCC
set(SHARED_CLANG_AND_CLANG_CL_AND_GCC_OPTIONS
        -Wcast-align
        -Wconversion
        -Wdouble-promotion
        -Wextra
        -Wfloat-equal
        -Wformat=2
        -Wimplicit-fallthrough
        -Wmisleading-indentation
        -Wnon-virtual-dtor
        -Wnull-dereference
        -Wold-style-cast
        -Woverloaded-virtual
        -Wshadow
        -Wsign-conversion
        -Wunused
)

# Clang specific options
set(DEFAULT_CLANG_OPTIONS
        ${SHARED_CLANG_AND_CLANG_CL_AND_GCC_OPTIONS}
        -Wall
        -Wmost
        -Wpedantic
        -Wshadow-all
)

# GCC specific options
set(DEFAULT_GCC_OPTIONS
        ${SHARED_CLANG_AND_CLANG_CL_AND_GCC_OPTIONS}
        -Wall
        -Wduplicated-branches
        -Wduplicated-cond
        -Wlogical-op
        -Wpedantic
        -Wuseless-cast
)

# MSVC specific options
set(DEFAULT_MSVC_OPTIONS
        /W4
        /w14062  # C4062: enumerator 'identifier' in switch of enum 'enumeration' is not handled
        /w14242  # C4242: 'identifier': conversion from 'type1' to 'type2', possible loss of data
        /w14254  # C4254: 'operator': conversion from 'type1' to 'type2', possible loss of data
        /w14263  # C4263: 'function': member function does not override any base class virtual member function
        /w14265  # C4265: 'class': class has virtual functions, but destructor is not virtual
        /w14287  # C4287: 'operator': unsigned/negative constant mismatch
        /w14289  # C4289: nonstandard extension used : 'var' : loop control variable declared in the for-loop is used outside the for-loop scope
        /w14296  # C4296: 'operator': expression is always false
        /w14311  # C4311: 'variable': pointer truncation from 'type' to 'type'
        /w14545  # C4545: expression before comma evaluates to a function which is missing an argument list
        /w14546  # C4546: function call before comma missing argument list
        /w14547  # C4547: 'operator': operator before comma has no effect; expected operator with side-effect
        /w14549  # C4549: 'operator1': operator before comma has no effect; did you intend 'operator2'?
        /w14555  # C4555: expression has no effect; expected expression with side-effect
        /w14619  # C4619: #pragma warning: there is no warning number 'number'
        /w14640  # C4640: 'instance': construction of local static object is not thread-safe
        /w14826  # C4826: Conversion from 'type1' to 'type2' is sign-extended. This may cause unexpected runtime behavior.
        /w14905  # C4905: wide string literal cast to 'LPSTR'
        /w14906  # C4906: string literal cast to 'LPWSTR'
        /w14928  # C4928: illegal copy-initialization; more than one user-defined conversion has been implicitly applied
)

# Clang-cl specific options
set(DEFAULT_CLANG_CL_OPTIONS
        ${SHARED_CLANG_AND_CLANG_CL_AND_GCC_OPTIONS}
        ${DEFAULT_MSVC_OPTIONS}
        -Wmost
        -Wpedantic
        -Wshadow-all
)

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND NOT CMAKE_CXX_COMPILER_FRONTEND_VARIANT STREQUAL "MSVC")
    # Clang
    set(DEFAULT_COMPILER_WARNINGS ${DEFAULT_CLANG_OPTIONS})
    set(DEFAULT_COMPILER_WARNINGS_DISABLED "-w")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND CMAKE_CXX_COMPILER_FRONTEND_VARIANT STREQUAL "MSVC")
    # Clang-cl
    set(DEFAULT_COMPILER_WARNINGS ${DEFAULT_CLANG_CL_OPTIONS})
    set(DEFAULT_COMPILER_WARNINGS_DISABLED "/w")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    # GCC
    set(DEFAULT_COMPILER_WARNINGS ${DEFAULT_GCC_OPTIONS})
    set(DEFAULT_COMPILER_WARNINGS_DISABLED "-w")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    # MSVC
    string(REGEX REPLACE " /W[0-4]" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    set(DEFAULT_COMPILER_WARNINGS ${DEFAULT_MSVC_OPTIONS})
    set(DEFAULT_COMPILER_WARNINGS_DISABLED "/w")
endif ()
