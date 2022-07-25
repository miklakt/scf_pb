//$ c++ -O3 -Wall -shared -std=c++14 -fPIC $(python3 -m pybind11 --includes) py.cpp -o py$(python3-config --extension-suffix)
#include <pybind11/pybind11.h>
#include "ss_scf_common.hpp"
#include "normalization_condition.hpp"

namespace py = pybind11;

int add(int i, int j) {
    return i + j;
}

PYBIND11_MODULE(scf_pb, m){
    m.doc() = "Analytical self-consistent filed for polymer brushes";

    m.def("add", &add, "A function that adds two numbers");

    //m.def("")
}
