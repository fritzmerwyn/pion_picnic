cmake_minimum_required(VERSION 3.15)
project(maintest2)


set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED YES)   # See below (1)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -O3")
add_executable(maintest2
                 main_test.cpp
                 Dyson_test_2_motherang.cpp
                 help_numerics.cpp
                 bse_root.cpp
                 lookuptable.cpp)


find_package(Eigen3 3.3 REQUIRED)    # See below (2)
find_package(GSL REQUIRED)
find_package(/usr/include/spectra/cmake/spectra)
target_link_libraries(maintest2 GSL::gsl GSL::gslcblas Eigen3::Eigen Spectra::Spectra)
