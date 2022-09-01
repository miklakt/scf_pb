cmake_minimum_required(VERSION 3.4...3.18)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(scf_pb VERSION 0.0.1)

add_subdirectory(pybind11)
pybind11_add_module(_scf_pb src/pybind.cpp)

#set_target_properties(_scf_pb
#    PROPERTIES
#    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/scf_pb"
#    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/scf_pb"
#    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/scf_pb"
#)

#set_target_properties(_scf_pb PROPERTIES OUTPUT_NAME scf_pb._scf_pb)