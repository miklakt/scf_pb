#pragma once
#include "ss_scf_common.hpp"
#include "normalization_condition.hpp"
#include "particles.hpp"
#include "profiles.hpp"

#define INTEGRATE_PROFILE(F, a, b) boost::math::quadrature::gauss_kronrod<double, 31>::integrate(F, a, b);

const double chi_crit = 6.0*std::log(5.0/6.0);
double gamma_2poly_model(const double a1, const double a2, const double chi, const double chi_PC, const double phi){
    double chi_ads = chi_PC - chi*(1-phi);
    double psi = a1*phi + a2*phi*phi;
    double gamma = (chi_ads-chi_crit)*psi;
    return gamma;
}

const double mobility_factor(const double phi, const double d, const double k){
    if (phi == 0){
        return 1;
    }
    double eps = 1/phi;
    double x = eps*eps/(d*d);
    return x/std::pow((1+std::pow(x,k)), 1/k);
}

class InsertedParticleProfiles{
    public:
        const Particle* particle;
        const BrushProfile* brush;
        InsertedParticleProfiles(const Particle* particle_, const BrushProfile* brush_): particle(particle_), brush(brush_){}
        double osmotic_free_energy(const double particle_center){
            auto integrand = [this, particle_center](const double z){
                    return brush->Pi_z(particle_center - particle->height/2+z)*particle->volume_integrand(z);
            };
            return INTEGRATE_PROFILE(integrand, 0, particle->height);
        }

};