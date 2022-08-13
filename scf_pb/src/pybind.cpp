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

double D_eff_mobility(const double N, const double sigma, const double chi, const double chi_PC, const double a0, const double a1, const double particle_width, const double particle_height){
    double kappa = topology::kappa(N);
    BrushProfilePlanar brush(chi, N, sigma, kappa);
    particle::Cylinder particle(particle_height, particle_width);
    auto gamma_phi = make_function::gamma_phi(a0, a1, chi, chi_PC);
    auto mobility = make_function::mobility_phi(particle_width);
    double d_eff = effective_diffusion_coefficient(&brush, &particle, gamma_phi, mobility);
    return d_eff;
}

double phi(const double N, const double sigma, const double chi, const double z){
    double kappa = topology::kappa(N);
    BrushProfilePlanar brush(chi, N, sigma, kappa);
    return brush.phi_z(z);
}

double D(const double N, const double sigma, const double chi){
    double kappa = topology::kappa(N);
    BrushProfilePlanar brush(chi, N, sigma, kappa);
    return brush.D();
}

PYBIND11_MODULE(scf_pb, m){
    m.doc() = "Analytical self-consistent filed for polymer brushes, calculate polymer density profile, insertion free energy penalty and diffusion coefficient";

    m.def("D_eff", &D_eff, "Effective diffusion coefficient");
    m.def("D_eff_corrected", &D_eff_mobility, "Corrected effective diffusion coefficient");
    m.def("phi", &phi, "Polymer density profile");
    m.def("D", &D, "Brush thickness");
}
