cmake_minimum_required(VERSION 3.14)

SET(SWIFTC_OPTIONS -g)

if(CMAKE_BUILD_TYPE STREQUAL "Release")
    SET(SWIFTC_OPTIONS -O -whole-module-optimization)
endif()

if(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    SET(SWIFTC_OPTIONS -O -whole-module-optimization -g)
endif()

if(CMAKE_BUILD_TYPE STREQUAL "MinSizeRel")
    SET(SWIFTC_OPTIONS -Osize -whole-module-optimization)
endif()

add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/mandelbrot_swift
                   COMMAND swiftc ${SWIFTC_OPTIONS} main.swift mandelbrot.swift -o ${CMAKE_CURRENT_BINARY_DIR}/mandelbrot_swift
                   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_custom_target(mandelbrot_swift_binary ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/mandelbrot_swift)
