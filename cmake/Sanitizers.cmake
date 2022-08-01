# unsupported sanitizers for GCC
set(GCC_UNSUPPORTED_SANITIZERS
    memory
)

option(USE_SANITIZER "Enable sanitizer. Options are: OFF, address, leak, thread, undefined, memory" OFF)

string(TOLOWER ${USE_SANITIZER} USE_SANITIZER)

if(USE_SANITIZER STREQUAL "address")
    set(SANITIZER_FLAGS -fsanitize=address -fno-omit-frame-pointer)
elseif(USE_SANITIZER STREQUAL "leak")
    set(SANITIZER_FLAGS -fsanitize=leak)
elseif(USE_SANITIZER STREQUAL "thread")
    set(SANITIZER_FLAGS -fsanitize=thread)
elseif(USE_SANITIZER STREQUAL "undefined")
    set(SANITIZER_FLAGS -fsanitize=undefined -fno-omit-frame-pointer)
elseif(USE_SANITIZER STREQUAL "memory")
    set(SANITIZER_FLAGS -fsanitize=memory -fsanitize-memory-track-origins -fno-omit-frame-pointer -fPIE -pie)
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    if(USE_SANITIZER IN_LIST GCC_UNSUPPORTED_SANITIZERS)
        message(SEND_ERROR "GCC does not support sanitizer: ${USE_SANITIZER}")
    endif()
endif()

if(SANITIZER_FLAGS)
    message(STATUS "Using sanitizer: " ${USE_SANITIZER})
endif()
