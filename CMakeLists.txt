cmake_minimum_required(VERSION 3.16)
project(barnes_hut_simulation VERSION 1.0.0 LANGUAGES CXX)

option(BARNES_HUT_ENABLE_OPENMP "Enable OpenMP parallelization when available" ON)
option(BARNES_HUT_BUILD_TESTS "Build the test executable" ON)

function(barnes_hut_enable_warnings target)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)
    endif()
endfunction()

add_library(barnes_hut
    src/simulation.cpp
)
target_include_directories(barnes_hut PUBLIC include)
target_compile_features(barnes_hut PUBLIC cxx_std_17)
barnes_hut_enable_warnings(barnes_hut)

if(BARNES_HUT_ENABLE_OPENMP)
    find_package(OpenMP QUIET)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(barnes_hut PUBLIC OpenMP::OpenMP_CXX)
        message(STATUS "OpenMP enabled")
    else()
        message(STATUS "OpenMP not found; building serial version")
    endif()
endif()

add_executable(barnes-hut src/main.cpp)
target_link_libraries(barnes-hut PRIVATE barnes_hut)
barnes_hut_enable_warnings(barnes-hut)

if(BARNES_HUT_BUILD_TESTS)
    enable_testing()
    add_executable(barnes-hut-tests tests/test_simulation.cpp)
    target_link_libraries(barnes-hut-tests PRIVATE barnes_hut)
    barnes_hut_enable_warnings(barnes-hut-tests)
    add_test(NAME simulation COMMAND barnes-hut-tests)
endif()
