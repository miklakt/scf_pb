#pragma once
#include "ss_scf_common.hpp"
#include "normalization_condition.hpp"
#include "particle_geometry.hpp"
#include "brush_profiles.hpp"

#include "iostream"

#define INTEGRATE_FUNC(F, a, b) boost::math::quadrature::gauss_kronrod<double, 31>::integrate(F, a, b);

namespace surface_interaction_coefficient{
    const double chi_crit = 6.0*std::log(5.0/6.0);
    
    double gamma_2poly_model(const double a0, const double a1, const double chi, const double chi_PC, const double phi){
        double chi_ads = chi_PC - chi*(1-phi);
        double psi = a0*phi + a1*phi*phi;
        double gamma = (chi_ads-chi_crit)*psi;
    return gamma;
    }

    auto gamma_phi = [](const double a0, const double a1, const double chi, const double chi_PC){return [=](const double phi){return gamma_2poly_model(a0, a1, chi, chi_PC, phi);};};
}//namespace surface_interaction_coefficient

namespace particle_mobility{
    double mobility_factor(const double phi, const double d, const double k){
        if (phi == 0){
            return 1;
        }
        double eps = 1/phi;
        double x = eps*eps/(d*d);
        return x/std::pow((1+std::pow(x,k)), 1/k);
    }
    
    auto mobility_phi = [](const double d, const double k){return [=](const double phi){return mobility_factor(phi, d, k);};};
}//namespace particle_mobility


template<typename ParticleType, typename BrushType, typename SurfaceInteractionModel>
class ParticleBrushInteractionEnergy{
    private:
        const ParticleType& particle;
        const BrushType& brush;
        const SurfaceInteractionModel gamma_phi;
    public:
        ParticleBrushInteractionEnergy(const ParticleType& particle_, const BrushType& brush_, SurfaceInteractionModel gamma_phi_)
        : particle(particle_), brush(brush_), gamma_phi(gamma_phi_){}

        double osmotic_free_energy(const double particle_center){
            auto integrand = [this, &particle_center](const double z){
                    const double z0 = particle_center-particle->height/2;
                    return brush->Pi_z(z0+z)*particle->volume_integrand(z);
            };
            return INTEGRATE_FUNC(integrand, 0, particle->height);
        }

        double surface_free_energy(const double particle_center)
        {
            const double A0 = particle->surface_edges()[0];
            const double A1 = particle->surface_edges()[1];

            const double z0 = particle_center - particle->height/2;
            const double z1 = particle_center + particle->height/2;

            const double phi0 = brush->phi_z(z0);
            const double phi1 = brush->phi_z(z1);

            const double gamma0 = gamma_phi(phi0);
            const double gamma1 = gamma_phi(phi1);

            const double free_energy_at_boundaries = A0*phi0*gamma0+A1*phi1*gamma1;

            auto integrand = [this, &z0](const double arg){
                const double phi = brush->phi_z(arg+z0);
                const double gamma = gamma_phi(phi);
                return phi*gamma*particle->surface_integrand(arg);
                };

            return INTEGRATE_FUNC(integrand, 0.0, particle->height) + free_energy_at_boundaries;
        }

        double total_free_energy(const double particle_center)
        {
            return surface_free_energy(particle_center) + osmotic_free_energy(particle_center);
        }

        double diffusion_coefficient(const double a, const double b){
            auto integrand = [this](const double z){return std::exp(total_free_energy(z));};
            return (b-a)/INTEGRATE_FUNC(integrand, a, b);
        }

        double diffusion_coefficient(){
            const double a = particle->height/2;
            const double b = brush->D()+particle->height/2;
            return diffusion_coefficient(a, b);
        }

        template<typename MobilityFunc>
        double diffusion_coefficient(const MobilityFunc mobility_phi, const double a, const double b){
            auto integrand = [this, mobility_phi](const double z){
                const double phi = brush->phi_z(z);
                const double mobility_factor = mobility_phi(phi);
                return std::exp(total_free_energy(z))/mobility_factor;
                };
            return (b-a)/INTEGRATE_FUNC(integrand, a, b);
        }

        template<typename MobilityFunc>
        double diffusion_coefficient(const MobilityFunc mobility_phi){
            const double a = particle->height/2;
            const double b = brush->D()+particle->height/2;
            return diffusion_coefficient(mobility_phi, a, b);
        }
};
