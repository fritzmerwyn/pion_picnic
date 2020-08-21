cmake_minimum_required(VERSION 3.15)
project(hyrdogen_matrixfree)

### set 'export' variables ###

set(PETSC_ARCH complex-mpi-nodebug)
set(SLEPC_DIR /home/felix/slepc-3.13.3)
set(PETSC_DIR /home/felix/petsc)

# find_package(Eigen3 3.3 REQUIRED)    # See below (2)
# find_package(GSL REQUIRED)
# find_package(Spectra REQUIRED)


### include slepc/petsc diretctories ###

include_directories(${SLEPC_DIR}/include)
include_directories(${SLEPC_DIR}/${PETSC_ARCH}/include)
include_directories(${PETSC_DIR}/include)
include_directories(${PETSC_DIR}/${PETSC_ARCH}/include)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED YES)   # See below (1)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -O3")


add_executable(hydrogen_matrixfree
                 hydrogen_matrixfree.cpp)

### link directories/libraries to target file ###

target_link_directories( hydrogen_matrixfree PUBLIC ${SLEPC_DIR}/${PETSC_ARCH}/lib PUBLIC ${PETSC_DIR}/${PETSC_ARCH}/lib )
target_link_libraries(hydrogen_matrixfree petsc slepc)
