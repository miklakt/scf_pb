//$ c++ -O3 -Wall -shared -std=c++14 -fPIC $(python3 -m pybind11 --includes) py.cpp -o py$(python3-config --extension-suffix)
//#pragma once
#include <pybind11/pybind11.h>
#include "ss_scf_common.hpp"
#include "normalization_condition.hpp"
#include "profiles.hpp"
#include "topology.hpp"
#include "particles.hpp"
#include "energy.hpp"

namespace py = pybind11;

double D_eff(const double N, const double sigma, const double chi, const double chi_PC, const double a0, const double a1, const double particle_width, const double particle_height){
    double kappa = topology::kappa(N);
    BrushProfilePlanar brush(chi, N, sigma, kappa);
    particle::Cylinder particle(particle_height, particle_width);
    auto gamma_phi = make_function::gamma_phi(a0, a1, chi, chi_PC);
    double d_eff = effective_diffusion_coefficient(&brush, &particle, gamma_phi);
    return d_eff;
}

PYBIND11_MODULE(scf_pb, m){
    m.doc() = "Analytical self-consistent filed for polymer brushes";

    m.def("D_eff", &D_eff, "Effective diffusion coefficient");
}
