#pragma once
#include <algorithm>
#include <stdexcept>

#include "ss_scf_common.hpp"
#include "normalization_condition.hpp"
#include "particle_geometry.hpp"
#include "brush_profiles.hpp"

#define INTEGRATE_FUNC(F, a, b) boost::math::quadrature::gauss_kronrod<double, 31>::integrate(F, a, b)

namespace surface_interaction_coefficient
{
    const double chi_crit = 6.0 * std::log(5.0 / 6.0);

    double gamma_2poly_model(const double a0, const double a1, const double chi, const double chi_PC, const double phi)
    {
        double chi_ads = chi_PC - chi * (1 - phi);
        double psi = a0 * phi + a1 * phi * phi;
        double gamma = (chi_ads - chi_crit) * psi;
        return gamma;
    }

    auto gamma_phi = [](const double a0, const double a1, const double chi, const double chi_PC)
    { return [=](const double phi)
      { return gamma_2poly_model(a0, a1, chi, chi_PC, phi); }; };
} // namespace surface_interaction_coefficient

namespace particle_mobility
{
    double mobility_factor(const double phi, const double d, const double k)
    {
        if (phi == 0)
        {
            return 1;
        }
        double eps = 1 / phi;
        double x = eps * eps / (d * d);
        return x / std::pow((1 + std::pow(x, k)), 1 / k);
    }

    auto mobility_phi = [](const double d, const double k)
    { return [=](const double phi)
      { return mobility_factor(phi, d, k); }; };
} // namespace particle_mobility


template <typename ParticleType, typename BrushType, typename SurfaceInteractionModel>
class ParticleBrushInteractionEnergy
{
private:
    const ParticleType *particle;
    const BrushType *brush;
    const SurfaceInteractionModel gamma_phi;

public:
    ParticleBrushInteractionEnergy(const ParticleType *particle_, const BrushType *brush_, SurfaceInteractionModel gamma_phi_)
        : particle(particle_), brush(brush_), gamma_phi(gamma_phi_) {}

    double osmotic_free_energy(const double particle_center);

    double surface_free_energy(const double particle_center);

    double total_free_energy(const double particle_center)
    {
        return surface_free_energy(particle_center) + osmotic_free_energy(particle_center);
    }

    double diffusion_coefficient(const double a, const double b);

    double diffusion_coefficient();

    template <typename MobilityFunc>
    double diffusion_coefficient(const MobilityFunc mobility_phi, const double a, const double b);

    template <typename MobilityFunc>
    double diffusion_coefficient(const MobilityFunc mobility_phi);

    double partition_coefficient(const double a, const double b);
    double partition_coefficient(const double b);
    double partition_coefficient();

    double partition_coefficient_open();
    double partition_coefficient_slit(const double R);
};

// template implementation file
#include "particle_energy.tpp"
// TODO:
// ParticleBrushInteractionEnergy<Cylinder, BrushProfilePlanar> specialization
// Considering the fact of integrating over constant cross-section area
