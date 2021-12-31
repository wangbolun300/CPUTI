cmake_minimum_required(VERSION 3.8)
project(SIMPLE_TI)

################################################################################
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

### Configuration
set(SIMPLE_IT_EXTERNAL "${CMAKE_CURRENT_SOURCE_DIR}/external")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

include(Warnings)
include(UseColors)
include(SIMPLE_TI_Utils)

include(${PROJECT_NAME}_Dependencies)


# inclusion CCD
# add_subdirectory(src)
# prepend_current_path(SIMPLETI_SOURCES)
# cputi_copy_headers(${SIMPLETI_SOURCES})
# cputi_set_source_group(${SIMPLETI_SOURCES})
# add_library(simple_ti ${SIMPLETI_SOURCES})
# target_compile_definitions(simple_ti PUBLIC NOMINMAX)
# target_include_directories(simple_ti PUBLIC ${PROJECT_BINARY_DIR}/include)
#target_link_libraries(simple_ti PUBLIC Eigen3::Eigen)
# target_link_libraries(simple_ti PRIVATE warnings::all)

#Optional
add_executable(
        test_bin 
        app/test.cpp 
        app/io.cpp 
        app/read_rational_csv.cpp 
        app/new_ccd.cpp)
# target_link_libraries(test_bin PRIVATE simple_ti)

# gmp
find_package(GMP)
IF(NOT ${GMP_FOUND})
        MESSAGE(FATAL_ERROR "Cannot find GMP")
ENDIF()
target_include_directories(test_bin PUBLIC ${GMP_INCLUDE_DIRS})
target_link_libraries(test_bin PUBLIC ${GMP_LIBRARIES})

include(tbb)
target_link_libraries(test_bin PUBLIC TBB::tbb)
