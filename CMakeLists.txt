cmake_minimum_required (VERSION 3.1)

# compiler options
set(CMAKE_C_COMPILER "mpicc")
set(CMAKE_CXX_COMPILER "mpic++")
set(CMAKE_VERBOSE_MAKEFILE off)
set(CMAKE_BUILD_TYPE Debug)

project(Coral)

# standard flags
set(CMAKE_INSTALL_PREFIX /home/brbass/research/coral/bin)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -o3 -fopenmp")

# gsl flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -L/home/software/rhel6/gsl/1.15/lib -lgsl -lgslcblas -I/home/software/rhel6/gsl/1.15/include")

# trilinos flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lepetra -lamesos -laztecoo -lepetraext -ltriutils -lepetra -lteuchosremainder -lteuchosnumerics -lteuchoscomm -lteuchosparameterlist -lteuchoscore -I/usr/cac/rhel6/trilinos/12.0.1/include -L/usr/cac/rhel6/trilinos/12.0.1/lib")

file(GLOB coral_src "src/*.cc")

set(matrix_main test/matrix_solution.cc)

include_directories(src)

add_executable(matrix_solution ${coral_src} ${matrix_main})

install(TARGETS matrix_solution DESTINATION .)
