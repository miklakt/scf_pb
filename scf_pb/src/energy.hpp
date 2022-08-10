#pragma once
#include "ss_scf_common.hpp"
#include "normalization_condition.hpp"
#include <boost/math/quadrature/trapezoidal.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include "profiles.hpp"
#include "particles.hpp"



const double chi_crit = 6.0*std::log(5.0/6.0);
double gamma_2poly_model(const double a1, const double a2, const double chi, const double chi_PC, const double phi){
    double chi_ads = chi_PC - chi*(1-phi);
    double psi = a1*phi + a2*phi*phi;
    double gamma = (chi_ads-chi_crit)*psi;
    return gamma;
}


namespace make_function{
auto gamma_phi = [](const double a1, const double a2, const double chi, const double chi_PC){return [=](const double phi){return gamma_2poly_model(a1, a2, chi, chi_PC, phi);};};

template<typename BrushType, typename ParticleType>
auto osmotic_free_energy_func(const BrushType* brush, const ParticleType* particle){
    
    auto integrand = [brush, particle](const double particle_center){
        return [brush, particle, particle_center](const double z){
            return brush->Pi_z(particle_center - particle->height/2+z)*particle->volume_integrand(z);
            };
        };
    auto integral = [particle, integrand](const double particle_center){
        return boost::math::quadrature::gauss_kronrod<double, 31>::integrate(integrand(particle_center), 0.0, particle->height);
        };

    return integral;
}

template<typename BrushType, typename ParticleType, typename SurfaceCoefficient>
auto surface_free_energy_func(const BrushType* brush, const ParticleType* particle, SurfaceCoefficient gamma)
{
    const double A0 = particle->surface_edges()[0];
    const double A1 = particle->surface_edges()[1];
    auto phi_z = [brush](const double arg){return brush->phi_z(arg);};
    auto gamma_z = [phi_z, gamma](const double arg){return gamma(phi_z(arg));};
    auto surface_integrand = [particle](const double arg){return particle->surface_integrand(arg);};
    
    auto free_energy = [=](const double particle_position){
        double z0 = particle_position - particle->height/2;
        double z1 = particle_position + particle->height/2;
        double phi0 = phi_z(z0);
        double phi1 = phi_z(z1);

        auto integrand = [phi_z, gamma_z, surface_integrand, z0](const double z){return phi_z(z+z0)*gamma_z(z+z0)*surface_integrand(z);};
        double I = boost::math::quadrature::gauss_kronrod<double, 31>::integrate(integrand, 0.0, particle->height); 
        return I + A0*phi0*gamma(phi0) + A1*phi1*gamma(phi1); 
    };
    return free_energy;
}

template<typename BrushType, typename ParticleType, typename SurfaceCoefficient>
auto total_free_energy_func(const BrushType* brush, const ParticleType* particle, SurfaceCoefficient gamma){
    auto free_energy = [=](const double  particle_position){return osmotic_free_energy_func(brush, particle)(particle_position)+surface_free_energy_func(brush, particle, gamma)(particle_position);};
    return free_energy;
}

template<typename BrushType, typename ParticleType, typename SurfaceCoefficient>
auto effective_diffusion_coefficient_ab(const BrushType* brush, const ParticleType* particle, SurfaceCoefficient gamma){
    auto free_energy_z = make_function::total_free_energy_func(brush, particle, gamma);
    auto integrand = [free_energy_z](const double z){return std::exp(free_energy_z(z));};
    auto d_eff = [=](const double a, const double b){return (b-a)/boost::math::quadrature::gauss_kronrod<double, 31>::integrate(integrand, a, b);};
    return d_eff;
}
}

template<typename BrushType, typename ParticleType, typename SurfaceCoefficient>
auto effective_diffusion_coefficient(const BrushType* brush, const ParticleType* particle, SurfaceCoefficient gamma){
    return make_function::effective_diffusion_coefficient_ab(brush, particle, gamma)(particle->height/2, brush->D());
}

