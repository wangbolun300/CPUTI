cmake_minimum_required(VERSION 3.14)
project(FAST_TI)

################################################################################
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(FAST_TI_TOPLEVEL_PROJECT ON)
else()
  set(FAST_TI_TOPLEVEL_PROJECT OFF)
endif()

OPTION(FAST_TI_USE_DOUBLE_PRECISION  "Use double precision"   ON)
OPTION(FAST_TI_CALCULATE_ERROR_BOUND "Calculate errors"       ON)


### Configuration
set(SIMPLE_IT_EXTERNAL "${CMAKE_CURRENT_SOURCE_DIR}/external")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

include(Warnings)
include(UseColors)

add_library(fast_ti
        src/cpu_ti/cpu_ti.cpp
        src/cpu_ti/timer.hpp
        src/cpu_ti/type.hpp
)
target_include_directories(fast_ti PUBLIC src)

if(FAST_TI_USE_DOUBLE_PRECISION)
        target_compile_definitions(fast_ti PUBLIC GPUTI_USE_DOUBLE_PRECISION)
endif()

if(FAST_TI_CALCULATE_ERROR_BOUND)
        target_compile_definitions(fast_ti PUBLIC CALCULATE_ERROR_BOUND)
endif()

include(tbb)
target_link_libraries(fast_ti PUBLIC TBB::tbb warnings::all)

# # Figure out AVX level support
# message(STATUS "Searching for AVX...")
# find_package(AVX)
# # Add SSE, AVX, and FMA flags to compiler flags
# string(REPLACE " " ";" SIMD_FLAGS "${AVX_FLAGS}")
# target_compile_options(test_bin PRIVATE ${SIMD_FLAGS})



if(FAST_TI_TOPLEVEL_PROJECT)
        add_executable(
                test_bin 
                app/test.cpp 
                app/io.cpp 
                app/read_rational_csv.cpp 
        )
        target_link_libraries(test_bin PUBLIC fast_ti)

        # gmp
        find_package(GMP)
        IF(NOT ${GMP_FOUND})
                MESSAGE(FATAL_ERROR "Cannot find GMP")
        ENDIF()
        target_include_directories(test_bin PUBLIC ${GMP_INCLUDE_DIRS})
        target_link_libraries(test_bin PUBLIC ${GMP_LIBRARIES})
endif()
