cmake_minimum_required(VERSION 3.4...3.18)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(scf_pb VERSION 0.0.1)

file(GLOB SRCFILES src/*.cpp)

add_subdirectory(pybind11)

pybind11_add_module(_scf_pb ${SRCFILES})