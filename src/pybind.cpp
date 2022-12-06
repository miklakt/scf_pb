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

double conc_profile(const double N, const double sigma, const double chi, const double chi_PC, const double a0, const double a1, const double particle_width, const double particle_height, const double k_smooth, const double l, const double z){
    double kappa = topology::kappa(N);
    BrushProfilePlanar brush(chi, N, sigma, kappa);
    particle::Sphere particle(particle_height/2);
    auto gamma_phi = surface_interaction_coefficient::gamma_2poly_model(a0, a1, chi, chi_PC);
    ParticleBrushInteractionEnergy particle_in_brush{&particle, &brush, gamma_phi};
    auto mobility_phi = particle_mobility::mobility_phi(particle_width, k_smooth);
    double source_dist = brush.D()+l;
    double conc = particle_in_brush.particle_concentration(mobility_phi, source_dist, z);
    return conc;
}

double D_z(const double N, const double sigma, const double chi, const double chi_PC, const double a0, const double a1, const double particle_width, const double particle_height, const double k_smooth, const double z){
    double kappa = topology::kappa(N);
    BrushProfilePlanar brush(chi, N, sigma, kappa);
    particle::Sphere particle(particle_height/2);
    auto gamma_phi = surface_interaction_coefficient::gamma_2poly_model(a0, a1, chi, chi_PC);
    ParticleBrushInteractionEnergy particle_in_brush{&particle, &brush, gamma_phi};
    auto mobility_phi = particle_mobility::mobility_phi(particle_width, k_smooth);
    double mobility_factor = particle_in_brush.mobility_factor(mobility_phi, z);
    return mobility_factor;
}

double D_eff(const double a, const double b, const double N, const double sigma, const double chi, const double chi_PC, const double a0, const double a1, const double particle_width, const double particle_height, const double k_smooth){
    double kappa = topology::kappa(N);
    BrushProfilePlanar brush(chi, N, sigma, kappa);
    particle::Sphere particle(particle_height/2);
    auto gamma_phi = surface_interaction_coefficient::gamma_2poly_model(a0, a1, chi, chi_PC);
    ParticleBrushInteractionEnergy particle_in_brush{&particle, &brush, gamma_phi};
    auto mobility_phi = particle_mobility::mobility_phi(particle_width, k_smooth);
    double d_eff = particle_in_brush.diffusion_coefficient(mobility_phi, a, b);
    return d_eff;
}

double D_eff_uncorrected(const double a, const double b, const double N, const double sigma, const double chi, const double chi_PC, const double a0, const double a1, const double particle_width, const double particle_height){
    double kappa = topology::kappa(N);
    BrushProfilePlanar brush(chi, N, sigma, kappa);
    particle::Sphere particle(particle_height/2);
    auto gamma_phi = surface_interaction_coefficient::gamma_2poly_model(a0, a1, chi, chi_PC);
    ParticleBrushInteractionEnergy particle_in_brush{&particle, &brush, gamma_phi};
    double d_eff = particle_in_brush.diffusion_coefficient(a, b);
    return d_eff;
}

double D_eff_no_energy(const double a, const double b, const double N, const double sigma, const double chi, const double chi_PC, const double a0, const double a1, const double particle_width, const double particle_height, const double k_smooth){
    double kappa = topology::kappa(N);
    BrushProfilePlanar brush(chi, N, sigma, kappa);
    particle::Sphere particle(particle_height/2);
    auto gamma_phi = surface_interaction_coefficient::gamma_2poly_model(a0, a1, chi, chi_PC);
    ParticleBrushInteractionEnergy particle_in_brush{&particle, &brush, gamma_phi};
    auto mobility_phi = particle_mobility::mobility_phi(particle_width, k_smooth);
    double d_eff = particle_in_brush.diffusion_coefficient_no_energy(mobility_phi, a, b);
    return d_eff;
}


double PC(const double a, const double b, const double N, const double sigma, const double chi, const double chi_PC, const double a0, const double a1, const double particle_width, const double particle_height){
    double kappa = topology::kappa(N);
    BrushProfilePlanar brush(chi, N, sigma, kappa);
    particle::Sphere particle(particle_height/2);
    auto gamma_phi = surface_interaction_coefficient::gamma_2poly_model(a0, a1, chi, chi_PC);
    ParticleBrushInteractionEnergy particle_in_brush{&particle, &brush, gamma_phi};
    double pc = particle_in_brush.partition_coefficient(a, b);
    return pc;
}

double free_energy(const double N, const double sigma, const double chi, const double chi_PC, const double a0, const double a1, const double particle_width, const double particle_height, const double z){
    double kappa = topology::kappa(N);
    BrushProfilePlanar brush(chi, N, sigma, kappa);
    particle::Sphere particle(particle_height/2);
    auto gamma_phi = surface_interaction_coefficient::gamma_2poly_model(a0, a1, chi, chi_PC);
    ParticleBrushInteractionEnergy particle_in_brush{&particle, &brush, gamma_phi};
    double fe = particle_in_brush.total_free_energy(z);
    return fe;
}

auto free_energy_all(const double N, const double sigma, const double chi, const double chi_PC, const double a0, const double a1, const double particle_width, const double particle_height, const double z){
    double kappa = topology::kappa(N);
    BrushProfilePlanar brush(chi, N, sigma, kappa);
    particle::Sphere particle(particle_height/2);
    auto gamma_phi = surface_interaction_coefficient::gamma_2poly_model(a0, a1, chi, chi_PC);
    ParticleBrushInteractionEnergy particle_in_brush{&particle, &brush, gamma_phi};
    double fe_osm = particle_in_brush.osmotic_free_energy(z);
    double fe_sur = particle_in_brush.surface_free_energy(z);
    double fe = fe_osm+fe_sur;
    return std::make_tuple(fe, fe_osm, fe_sur);
}

double free_energy_external(const std::vector<double> phi, const double chi, const double chi_PC, const double a0, const double a1, const double particle_width, const double particle_height, const double z){
    BrushProfileExternal brush(phi, chi);
    particle::Sphere particle(particle_height/2);
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

//Temporary here
double d_critical(const double a, const double b, const double N, const double sigma, const double chi, const double chi_PC, const double a0, const double a1, const double k_smooth, const double min, const double max){
    double kappa = topology::kappa(N);
    BrushProfilePlanar brush(chi, N, sigma, kappa);
    auto gamma_phi = surface_interaction_coefficient::gamma_2poly_model(a0, a1, chi, chi_PC);
    auto fsolve = [=](const double particle_height){
        particle::Sphere particle(particle_height/2);
        ParticleBrushInteractionEnergy particle_in_brush{&particle, &brush, gamma_phi};
        auto mobility_phi = particle_mobility::mobility_phi(particle_height, k_smooth);
        double d_eff = particle_in_brush.diffusion_coefficient(mobility_phi, a, b);
        return std::log(d_eff)-1.0;
    };

    const double maxit = 50;
    std::uintmax_t it = maxit;
    int digits = std::numeric_limits<double>::digits;
    int get_digits = digits - 3;
    eps_tolerance<double> tol(get_digits);


    std::pair<double, double> r = toms748_solve(fsolve, min, max, tol, it);

    return r.first + (r.second - r.first)/2;
}

double chi_PC_critical(const double a, const double b, const double N, const double sigma, const double chi, const double a0, const double a1, const double d, const double k_smooth, const double min, const double max){
    double kappa = topology::kappa(N);
    BrushProfilePlanar brush(chi, N, sigma, kappa);
    particle::Sphere particle(d/2);
    auto fsolve = [=](const double chi_PC){
        auto gamma_phi = surface_interaction_coefficient::gamma_2poly_model(a0, a1, chi, chi_PC);
        ParticleBrushInteractionEnergy particle_in_brush{&particle, &brush, gamma_phi};
        auto mobility_phi = particle_mobility::mobility_phi(d, k_smooth);
        double d_eff = particle_in_brush.diffusion_coefficient(mobility_phi, a, b);
        return std::log(d_eff)-1.0;
    };

    const double maxit = 50;
    std::uintmax_t it = maxit;
    int digits = std::numeric_limits<double>::digits;
    int get_digits = digits - 3;
    eps_tolerance<double> tol(get_digits);

    const double chi_crit = 6.0 * std::log(5.0 / 6.0);
    const double D = brush->D();
    const double phi_av = N*sigma/D;
    const double chi_PC_max = chi_crit+chi_PS*(1-phi_av);


    std::pair<double, double> r = toms748_solve(fsolve, min, max, tol, it);

    return r.first + (r.second - r.first)/2;
}

using namespace pybind11::literals;
PYBIND11_MODULE(_scf_pb, m){

    m.doc() = "Analytical self-consistent filed for polymer brushes, calculate polymer density profile, insertion free energy penalty and diffusion coefficient";

    m.def("phi", &phi, py::call_guard<py::gil_scoped_release>(), "Polymer density profile (cxx)",  py::kw_only{}, "N"_a, "sigma"_a, "chi"_a, "z"_a);
    m.def("Pi", &Pi_, py::call_guard<py::gil_scoped_release>(), "Polymer osmotic profile (cxx)",  py::kw_only{}, "N"_a, "sigma"_a, "chi"_a, "z"_a);
    m.def("D", &D, py::call_guard<py::gil_scoped_release>(), "Polymer brush thickness (cxx)", py::kw_only{}, "N"_a, "sigma"_a, "chi"_a);

    m.def("D_z", &D_z, py::call_guard<py::gil_scoped_release>(), "Local diffusion coefficient (cxx)", py::kw_only{}, "N"_a, "sigma"_a, "chi"_a, "chi_PC"_a, "a0"_a, "a1"_a, "particle_width"_a, "particle_height"_a, "k_smooth"_a, "z"_a);
    m.def("conc_profile", &conc_profile, py::call_guard<py::gil_scoped_release>(), "Concentration profile when perfect sink grafting surface (cxx)", py::kw_only{}, "N"_a, "sigma"_a, "chi"_a, "chi_PC"_a, "a0"_a, "a1"_a, "particle_width"_a, "particle_height"_a, "k_smooth"_a, "l"_a, "z"_a);

    m.def("free_energy", &free_energy, py::call_guard<py::gil_scoped_release>(), "Insertion free energy profile (cxx)", py::kw_only{}, "N"_a, "sigma"_a, "chi"_a, "chi_PC"_a, "a0"_a, "a1"_a, "particle_width"_a, "particle_height"_a, "z"_a);
    m.def("free_energy_all", &free_energy_all, py::call_guard<py::gil_scoped_release>(), "Insertion free energy profile (cxx)", py::kw_only{}, "N"_a, "sigma"_a, "chi"_a, "chi_PC"_a, "a0"_a, "a1"_a, "particle_width"_a, "particle_height"_a, "z"_a);
    m.def("free_energy_external", &free_energy_external, py::call_guard<py::gil_scoped_release>(), "Insertion free energy profile (cxx)", py::kw_only{}, "phi"_a, "chi"_a, "chi_PC"_a, "a0"_a, "a1"_a, "particle_width"_a, "particle_height"_a, "z"_a);


    m.def("D_eff", &D_eff, py::call_guard<py::gil_scoped_release>(), "Effective diffusion coefficient through polymer brush membrane (cxx)", "a"_a, "b"_a, "N"_a, "sigma"_a, "chi"_a, py::kw_only{}, "chi_PC"_a, "a0"_a, "a1"_a, "particle_width"_a, "particle_height"_a, "k_smooth"_a);

    m.def("D_eff_uncorrected", &D_eff_uncorrected, py::call_guard<py::gil_scoped_release>(), "Effective diffusion coefficient through polymer brush membrane (cxx)", "a"_a, "b"_a, "N"_a, "sigma"_a, "chi"_a, py::kw_only{}, "chi_PC"_a, "a0"_a, "a1"_a, "particle_width"_a, "particle_height"_a);

    m.def("D_eff_no_energy", &D_eff_no_energy, py::call_guard<py::gil_scoped_release>(), "Effective diffusion coefficient through polymer brush membrane (cxx)", "a"_a, "b"_a, "N"_a, "sigma"_a, "chi"_a, py::kw_only{}, "chi_PC"_a, "a0"_a, "a1"_a, "particle_width"_a, "particle_height"_a, "k_smooth"_a);

    m.def("PC", &PC, py::call_guard<py::gil_scoped_release>(), "Partition coefficient of particles in polymer brush and a semi-infinite solution (cxx)", "a"_a, "b"_a, "N"_a, "sigma"_a, "chi"_a, py::kw_only{}, "chi_PC"_a, "a0"_a, "a1"_a, "particle_width"_a, "particle_height"_a);

    m.def("mobility_factor", &particle_mobility::mobility_factor, py::call_guard<py::gil_scoped_release>(), "Polymer network mobility factor", py::kw_only{}, "phi"_a, "d"_a, "k_smooth"_a);


    m.def("d_critical", &d_critical, py::call_guard<py::gil_scoped_release>(), "Effective diffusion coefficient through polymer brush membrane (cxx)", "a"_a, "b"_a, "N"_a, "sigma"_a, "chi"_a, py::kw_only{}, "chi_PC"_a, "a0"_a, "a1"_a, "k_smooth"_a, "dmin"_a, "dmax"_a);
    m.def("chi_PC_critical", &chi_PC_critical, py::call_guard<py::gil_scoped_release>(), "Effective diffusion coefficient through polymer brush membrane (cxx)", "a"_a, "b"_a, "N"_a, "sigma"_a, "chi"_a, py::kw_only{}, "a0"_a, "a1"_a, "d"_a, "k_smooth"_a, "chi_PC_min"_a, "chi_PC_max"_a);

}
