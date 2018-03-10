cmake_minimum_required(VERSION 3.9)
project(opbe_gh)

set(CMAKE_CXX_STANDARD 11)

include_directories(include)

add_executable(opbe_gh
        src/averagingproblem.cpp
        src/deltaproblem.cpp
        src/density.cpp
        src/fixedicproblem.cpp
        src/gsldriver.cpp
        src/hermitepolynomial.cpp
        src/main.cpp
        src/mkproblem.cpp
        src/modeindex.cpp
        src/opbeparameter.cpp
        src/problem.cpp
        src/runcontrol.cpp
        src/system.cpp)

