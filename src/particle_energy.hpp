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

    auto gamma_2poly_model(const double a0, const double a1, const double chi, const double chi_PC)
    {
        auto gama_phi = [=](const double phi)
        {
            double chi_ads = chi_PC - chi * (1 - phi);
            double psi = a0 * phi + a1 * phi * phi;
            double gamma = (chi_ads - chi_crit) * psi;
            return gamma;
        };
        return gama_phi;
    }
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

    double osmotic_free_energy(const double z0) const;

    double surface_free_energy(const double z0) const;

    double total_free_energy(const double z0) const
    {
        return surface_free_energy(z0) + osmotic_free_energy(z0);
    }

    template <typename MobilityFunc>
    double mobility_factor(const MobilityFunc mobility_phi, const double z0) const;

    // template <typename MobilityFunc>
    // double particle_concentration(const double z0);

    double partition_coefficient(const double a, const double b) const;

    double diffusion_coefficient(const double a, const double b) const;

    template <typename MobilityFunc>
    double diffusion_coefficient(const MobilityFunc mobility_phi, const double a, const double b) const;

    template <typename MobilityFunc>
    double particle_concentration(const MobilityFunc mobility_phi, const double l, const double z0, const double c_bulk) const;

    template <typename MobilityFunc>
    double diffusion_coefficient_no_energy(const MobilityFunc mobility_phi) const;
};

// template implementation file
#include "particle_energy.tpp"
