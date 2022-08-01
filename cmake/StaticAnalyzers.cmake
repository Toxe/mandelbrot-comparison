option(USE_CLANG_TIDY "Enable clang-tidy static analysis" OFF)

if(USE_CLANG_TIDY)
    find_program(CLANG_TIDY NAMES ${CMAKE_CXX_CLANG_TIDY} clang-tidy)

    if(CLANG_TIDY)
        set(CMAKE_CXX_CLANG_TIDY ${CLANG_TIDY})
        message(STATUS "Using clang-tidy: " ${CLANG_TIDY})
    else()
        message(SEND_ERROR "clang-tidy requested but not found!")
    endif()
endif()
