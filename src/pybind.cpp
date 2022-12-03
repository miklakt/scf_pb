//$ c++ -O3 -Wall -shared -std=c++14 -fPIC $(python3 -m pybind11 --includes) py.cpp -o py$(python3-config --extension-suffix)
//#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "ss_scf_common.hpp"
#include "normalization_condition.hpp"
#include "brush_profiles.hpp"
#include "topology.hpp"
#include "particle_geometry.hpp"
#include "particle_energy.hpp"
#include "vector"

namespace py = pybind11;

double conc_profile_fixed_source(const double N, const double sigma, const double chi, const double chi_PC, const double a0, const double a1, const double particle_width, const double particle_height, const double k_smooth, const double source_dist, const double z){
    double kappa = topology::kappa(N);
    BrushProfilePlanar brush(chi, N, sigma, kappa);
    particle::Cylinder particle(particle_height, particle_width);
    auto gamma_phi = surface_interaction_coefficient::gamma_2poly_model(a0, a1, chi, chi_PC);
    ParticleBrushInteractionEnergy particle_in_brush{&particle, &brush, gamma_phi};
    auto mobility_phi = particle_mobility::mobility_phi(particle_width, k_smooth);
    double conc = particle_in_brush.particle_concentration(mobility_phi, source_dist, z, 1.0);
    return conc;
}

double conc_profile(const double N, const double sigma, const double chi, const double chi_PC, const double a0, const double a1, const double particle_width, const double particle_height, const double k_smooth, const double l, const double z){
    double kappa = topology::kappa(N);
    BrushProfilePlanar brush(chi, N, sigma, kappa);
    particle::Cylinder particle(particle_height, particle_width);
    auto gamma_phi = surface_interaction_coefficient::gamma_2poly_model(a0, a1, chi, chi_PC);
    ParticleBrushInteractionEnergy particle_in_brush{&particle, &brush, gamma_phi};
    auto mobility_phi = particle_mobility::mobility_phi(particle_width, k_smooth);
    double source_dist = brush.D()+l;
    double conc = particle_in_brush.particle_concentration(mobility_phi, source_dist, z, 1.0);
    return conc;
}

double D_z(const double N, const double sigma, const double chi, const double chi_PC, const double a0, const double a1, const double particle_width, const double particle_height, const double k_smooth, const double z){
    double kappa = topology::kappa(N);
    BrushProfilePlanar brush(chi, N, sigma, kappa);
    particle::Cylinder particle(particle_height, particle_width);
    auto gamma_phi = surface_interaction_coefficient::gamma_2poly_model(a0, a1, chi, chi_PC);
    ParticleBrushInteractionEnergy particle_in_brush{&particle, &brush, gamma_phi};
    auto mobility_phi = particle_mobility::mobility_phi(particle_width, k_smooth);
    double mobility_factor = particle_in_brush.mobility_factor(mobility_phi, z);
    return mobility_factor;
}

double D_eff(const double a, const double b, const double N, const double sigma, const double chi, const double chi_PC, const double a0, const double a1, const double particle_width, const double particle_height, const double k_smooth){
    double kappa = topology::kappa(N);
    BrushProfilePlanar brush(chi, N, sigma, kappa);
    particle::Cylinder particle(particle_height, particle_width);
    auto gamma_phi = surface_interaction_coefficient::gamma_2poly_model(a0, a1, chi, chi_PC);
    ParticleBrushInteractionEnergy particle_in_brush{&particle, &brush, gamma_phi};
    auto mobility_phi = particle_mobility::mobility_phi(particle_width, k_smooth);
    double d_eff = particle_in_brush.diffusion_coefficient(mobility_phi, a, b);
    return d_eff;
}

double D_eff_uncorrected(const double a, const double b, const double N, const double sigma, const double chi, const double chi_PC, const double a0, const double a1, const double particle_width, const double particle_height){
    double kappa = topology::kappa(N);
    BrushProfilePlanar brush(chi, N, sigma, kappa);
    particle::Cylinder particle(particle_height, particle_width);
    auto gamma_phi = surface_interaction_coefficient::gamma_2poly_model(a0, a1, chi, chi_PC);
    ParticleBrushInteractionEnergy particle_in_brush{&particle, &brush, gamma_phi};
    double d_eff = particle_in_brush.diffusion_coefficient(a, b);
    return d_eff;
}

double D_eff_no_energy(const double N, const double sigma, const double chi, const double chi_PC, const double a0, const double a1, const double particle_width, const double particle_height, const double k_smooth){
    double kappa = topology::kappa(N);
    BrushProfilePlanar brush(chi, N, sigma, kappa);
    particle::Cylinder particle(particle_height, particle_width);
    auto gamma_phi = surface_interaction_coefficient::gamma_2poly_model(a0, a1, chi, chi_PC);
    ParticleBrushInteractionEnergy particle_in_brush{&particle, &brush, gamma_phi};
    auto mobility_phi = particle_mobility::mobility_phi(particle_width, k_smooth);
    double d_eff = particle_in_brush.diffusion_coefficient_no_energy(mobility_phi);
    return d_eff;
}


double PC(const double a, const double b, const double N, const double sigma, const double chi, const double chi_PC, const double a0, const double a1, const double particle_width, const double particle_height){
    double kappa = topology::kappa(N);
    BrushProfilePlanar brush(chi, N, sigma, kappa);
    particle::Cylinder particle(particle_height, particle_width);
    auto gamma_phi = surface_interaction_coefficient::gamma_2poly_model(a0, a1, chi, chi_PC);
    ParticleBrushInteractionEnergy particle_in_brush{&particle, &brush, gamma_phi};
    double pc = particle_in_brush.partition_coefficient(a, b);
    return pc;
}

double PC_perfect_sink_fixed_source(const double N, const double sigma, const double chi, const double chi_PC, const double a0, const double a1, const double particle_width, const double particle_height, const double k_smooth, const double source_dist){
    double kappa = topology::kappa(N);
    BrushProfilePlanar brush(chi, N, sigma, kappa);
    particle::Cylinder particle(particle_height, particle_width);
    auto gamma_phi = surface_interaction_coefficient::gamma_2poly_model(a0, a1, chi, chi_PC);
    ParticleBrushInteractionEnergy particle_in_brush{&particle, &brush, gamma_phi};
    auto mobility_phi = particle_mobility::mobility_phi(particle_width, k_smooth);
    double pc = particle_in_brush.partition_coefficient_perfect_sink(mobility_phi,source_dist, 1.0);
    return pc;
}

double PC_perfect_sink(const double N, const double sigma, const double chi, const double chi_PC, const double a0, const double a1, const double particle_width, const double particle_height, const double k_smooth, const double l){
    double kappa = topology::kappa(N);
    BrushProfilePlanar brush(chi, N, sigma, kappa);
    particle::Cylinder particle(particle_height, particle_width);
    auto gamma_phi = surface_interaction_coefficient::gamma_2poly_model(a0, a1, chi, chi_PC);
    ParticleBrushInteractionEnergy particle_in_brush{&particle, &brush, gamma_phi};
    auto mobility_phi = particle_mobility::mobility_phi(particle_width, k_smooth);
    double source_dist = brush.D()+l;
    double pc = particle_in_brush.partition_coefficient_perfect_sink(mobility_phi,source_dist, 1.0);
    return pc;
}



double D_eff_external(const std::vector<double> phi, const double chi, const double chi_PC, const double a0, const double a1, const double particle_width, const double particle_height, const double k_smooth, const double a, const double b){
    BrushProfileExternal brush(phi, chi);
    particle::Cylinder particle(particle_height, particle_width);
    auto gamma_phi = surface_interaction_coefficient::gamma_2poly_model(a0, a1, chi, chi_PC);
    ParticleBrushInteractionEnergy particle_in_brush{&particle, &brush, gamma_phi};
    auto mobility_phi = particle_mobility::mobility_phi(particle_width, k_smooth);
    double d_eff = particle_in_brush.diffusion_coefficient(mobility_phi, a, b);
    return d_eff;
}

double free_energy(const double N, const double sigma, const double chi, const double chi_PC, const double a0, const double a1, const double particle_width, const double particle_height, const double z){
    double kappa = topology::kappa(N);
    BrushProfilePlanar brush(chi, N, sigma, kappa);
    particle::Cylinder particle(particle_height, particle_width);
    auto gamma_phi = surface_interaction_coefficient::gamma_2poly_model(a0, a1, chi, chi_PC);
    ParticleBrushInteractionEnergy particle_in_brush{&particle, &brush, gamma_phi};
    double fe = particle_in_brush.total_free_energy(z);
    return fe;
}

auto free_energy_all(const double N, const double sigma, const double chi, const double chi_PC, const double a0, const double a1, const double particle_width, const double particle_height, const double z){
    double kappa = topology::kappa(N);
    BrushProfilePlanar brush(chi, N, sigma, kappa);
    particle::Cylinder particle(particle_height, particle_width);
    auto gamma_phi = surface_interaction_coefficient::gamma_2poly_model(a0, a1, chi, chi_PC);
    ParticleBrushInteractionEnergy particle_in_brush{&particle, &brush, gamma_phi};
    double fe_osm = particle_in_brush.osmotic_free_energy(z);
    double fe_sur = particle_in_brush.surface_free_energy(z);
    double fe = fe_osm+fe_sur;
    return std::make_tuple(fe, fe_osm, fe_sur);
}

double free_energy_external(const std::vector<double> phi, const double chi, const double chi_PC, const double a0, const double a1, const double particle_width, const double particle_height, const double z){
    BrushProfileExternal brush(phi, chi);
    particle::Cylinder particle(particle_height, particle_width);
    auto gamma_phi = surface_interaction_coefficient::gamma_2poly_model(a0, a1, chi, chi_PC);
    ParticleBrushInteractionEnergy particle_in_brush{&particle, &brush, gamma_phi};
    double fe = particle_in_brush.total_free_energy(z);
    return fe;
}

double phi(const double N, const double sigma, const double chi, const double z){
    double kappa = topology::kappa(N);
    BrushProfilePlanar brush(chi, N, sigma, kappa);
    return brush.phi_z(z);
}

double Pi_(const double N, const double sigma, const double chi, const double z){
    double kappa = topology::kappa(N);
    BrushProfilePlanar brush(chi, N, sigma, kappa);
    return brush.Pi_z(z);
}

double D(const double N, const double sigma, const double chi){
    double kappa = topology::kappa(N);
    BrushProfilePlanar brush(chi, N, sigma, kappa);
    return brush.D();
}

using namespace pybind11::literals;
PYBIND11_MODULE(_scf_pb, m){

    m.doc() = "Analytical self-consistent filed for polymer brushes, calculate polymer density profile, insertion free energy penalty and diffusion coefficient";

    m.def("phi", &phi, py::call_guard<py::gil_scoped_release>(), "Polymer density profile (cxx)",  py::kw_only{}, "N"_a, "sigma"_a, "chi"_a, "z"_a);
    m.def("Pi", &Pi_, py::call_guard<py::gil_scoped_release>(), "Polymer osmotic profile (cxx)",  py::kw_only{}, "N"_a, "sigma"_a, "chi"_a, "z"_a);
    m.def("D", &D, py::call_guard<py::gil_scoped_release>(), "Polymer brush thickness (cxx)", py::kw_only{}, "N"_a, "sigma"_a, "chi"_a);

    m.def("D_z", &D_z, py::call_guard<py::gil_scoped_release>(), "Local diffusion coefficient (cxx)", py::kw_only{}, "N"_a, "sigma"_a, "chi"_a, "chi_PC"_a, "a0"_a, "a1"_a, "particle_width"_a, "particle_height"_a, "k_smooth"_a, "z"_a);
    m.def("conc_profile_fixed_source", &conc_profile_fixed_source, py::call_guard<py::gil_scoped_release>(), "Concentration profile when perfect sink grafting surface (cxx)", py::kw_only{}, "N"_a, "sigma"_a, "chi"_a, "chi_PC"_a, "a0"_a, "a1"_a, "particle_width"_a, "particle_height"_a, "k_smooth"_a, "source_dist"_a, "z"_a);
    m.def("conc_profile", &conc_profile, py::call_guard<py::gil_scoped_release>(), "Concentration profile when perfect sink grafting surface (cxx)", py::kw_only{}, "N"_a, "sigma"_a, "chi"_a, "chi_PC"_a, "a0"_a, "a1"_a, "particle_width"_a, "particle_height"_a, "k_smooth"_a, "l"_a, "z"_a);

    m.def("free_energy", &free_energy, py::call_guard<py::gil_scoped_release>(), "Insertion free energy profile (cxx)", py::kw_only{}, "N"_a, "sigma"_a, "chi"_a, "chi_PC"_a, "a0"_a, "a1"_a, "particle_width"_a, "particle_height"_a, "z"_a);
    m.def("free_energy_all", &free_energy_all, py::call_guard<py::gil_scoped_release>(), "Insertion free energy profile (cxx)", py::kw_only{}, "N"_a, "sigma"_a, "chi"_a, "chi_PC"_a, "a0"_a, "a1"_a, "particle_width"_a, "particle_height"_a, "z"_a);
    m.def("free_energy_external", &free_energy_external, py::call_guard<py::gil_scoped_release>(), "Insertion free energy profile (cxx)", py::kw_only{}, "phi"_a, "chi"_a, "chi_PC"_a, "a0"_a, "a1"_a, "particle_width"_a, "particle_height"_a, "z"_a);

    //m.def("mobility", &particle_mobility::mobility_factor, py::call_guard<py::gil_scoped_release>(), "Mobility factor", py::kw_only{}, "phi"_a, "d"_a, "k_smooth"_a);

    //m.def("D_eff_uncorrected", &D_eff_uncorrected, py::call_guard<py::gil_scoped_release>(), "Uncorrected effective diffusion coefficient (cxx)",  "N"_a, "sigma"_a, "chi"_a, py::kw_only{}, "chi_PC"_a, "a0"_a, "a1"_a, "particle_width"_a, "particle_height"_a);
    m.def("D_eff", &D_eff, py::call_guard<py::gil_scoped_release>(), "Effective diffusion coefficient through polymer brush membrane (cxx)", "a"_a, "b"_a, "N"_a, "sigma"_a, "chi"_a, py::kw_only{}, "chi_PC"_a, "a0"_a, "a1"_a, "particle_width"_a, "particle_height"_a, "k_smooth"_a);
    //m.def("D_eff_open", &D_eff_open, py::call_guard<py::gil_scoped_release>(), "Effective diffusion coefficient through polymer brush membrane (cxx)",  "N"_a, "sigma"_a, "chi"_a, py::kw_only{}, "chi_PC"_a, "a0"_a, "a1"_a, "particle_width"_a, "particle_height"_a, "k_smooth"_a);

    m.def("D_eff_uncorrected", &D_eff_uncorrected, py::call_guard<py::gil_scoped_release>(), "Effective diffusion coefficient through polymer brush membrane (cxx)", "a"_a, "b"_a, "N"_a, "sigma"_a, "chi"_a, py::kw_only{}, "chi_PC"_a, "a0"_a, "a1"_a, "particle_width"_a, "particle_height"_a);
    //m.def("D_eff_uncorrected_open", &D_eff_uncorrected_open, py::call_guard<py::gil_scoped_release>(), "Effective diffusion coefficient through polymer brush membrane (cxx)",  "N"_a, "sigma"_a, "chi"_a, py::kw_only{}, "chi_PC"_a, "a0"_a, "a1"_a, "particle_width"_a, "particle_height"_a);

    m.def("D_eff_no_energy", &D_eff_no_energy, py::call_guard<py::gil_scoped_release>(), "Effective diffusion coefficient through polymer brush membrane (cxx)",  "N"_a, "sigma"_a, "chi"_a, py::kw_only{}, "chi_PC"_a, "a0"_a, "a1"_a, "particle_width"_a, "particle_height"_a, "k_smooth"_a);
    m.def("D_eff_external", &D_eff_external, py::call_guard<py::gil_scoped_release>(), "Effective diffusion coefficient through polymer brush membrane (cxx)", py::kw_only{},  "phi"_a, "chi"_a, "chi_PC"_a, "a0"_a, "a1"_a, "particle_width"_a, "particle_height"_a, "k_smooth"_a, "a"_a, "b"_a);

    //m.def("PC_open", &PC_open, py::call_guard<py::gil_scoped_release>(), "Partition coefficient of particles in polymer brush and a semi-infinite solution (cxx)",  "N"_a, "sigma"_a, "chi"_a, py::kw_only{}, "chi_PC"_a, "a0"_a, "a1"_a, "particle_width"_a, "particle_height"_a);
    m.def("PC", &PC, py::call_guard<py::gil_scoped_release>(), "Partition coefficient of particles in polymer brush and a semi-infinite solution (cxx)", "a"_a, "b"_a, "N"_a, "sigma"_a, "chi"_a, py::kw_only{}, "chi_PC"_a, "a0"_a, "a1"_a, "particle_width"_a, "particle_height"_a);
    m.def("PC_perfect_sink_fixed_source", &PC_perfect_sink_fixed_source, py::call_guard<py::gil_scoped_release>(), "Partition coefficient of particles in polymer brush and a semi-infinite solution when perfect sink grafting surface (cxx)", py::kw_only{}, "N"_a, "sigma"_a, "chi"_a, "chi_PC"_a, "a0"_a, "a1"_a, "particle_width"_a, "particle_height"_a, "k_smooth"_a, "source_dist"_a);
    m.def("PC_perfect_sink", &PC_perfect_sink, py::call_guard<py::gil_scoped_release>(), "Partition coefficient of particles in polymer brush and a semi-infinite solution when perfect sink grafting surface (cxx)", py::kw_only{}, "N"_a, "sigma"_a, "chi"_a, "chi_PC"_a, "a0"_a, "a1"_a, "particle_width"_a, "particle_height"_a, "k_smooth"_a, "l"_a);

    m.def("mobility_factor", &particle_mobility::mobility_factor, py::call_guard<py::gil_scoped_release>(), "Polymer network mobility factor", py::kw_only{}, "phi"_a, "d"_a, "k_smooth"_a);
}
