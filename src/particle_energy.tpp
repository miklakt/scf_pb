template <typename ParticleType, typename BrushType, typename SurfaceInteractionModel>
double ParticleBrushInteractionEnergy<ParticleType, BrushType, SurfaceInteractionModel>::
    osmotic_free_energy(const double z0) const
{
    // particle is outside the brush
    if (z0 > brush->D())
    {
        return 0.0;
    }

    // const double right_limit = std::min(particle->height, brush->D()-z0);
    const double right_limit = particle->height;

    auto integrand = [this, &z0](const double arg)
    {
        return brush->Pi_z(arg + z0);
    };

    return particle::integrators::integrate_over_volume(particle, integrand, right_limit);
}

template <typename ParticleType, typename BrushType, typename SurfaceInteractionModel>
double ParticleBrushInteractionEnergy<ParticleType, BrushType, SurfaceInteractionModel>::
    surface_free_energy(const double z0) const
{
    // particle is outside the brush
    if (z0 > brush->D())
    {
        return 0.0;
    }

    // const double right_limit = std::min(particle->height, brush->D() - z0);
    const double right_limit = particle->height;

    auto integrand = [this, &z0](const double arg)
    {
        const double phi = brush->phi_z(arg + z0);
        const double gamma = this->gamma_phi(phi);
        return gamma; // bugfix return phi*gamma;!!!
    };

    return particle::integrators::integrate_over_surface(particle, integrand, right_limit);
}

template <typename ParticleType, typename BrushType, typename SurfaceInteractionModel>
template <typename MobilityFunc>
double ParticleBrushInteractionEnergy<ParticleType, BrushType, SurfaceInteractionModel>::
    mobility_factor(const MobilityFunc mobility_phi, const double z0) const
{
    const double phi = brush->phi_z(z0);
    return mobility_phi(phi);
}

template <typename ParticleType, typename BrushType, typename SurfaceInteractionModel>
double ParticleBrushInteractionEnergy<ParticleType, BrushType, SurfaceInteractionModel>::
    diffusion_coefficient(const double a, const double b) const
{
    if (a >= b)
    {
        throw std::invalid_argument("a>=b");
    }
    auto integrand = [this](const double z)
    { return std::exp(total_free_energy(z)); };
    double I = INTEGRATE_FUNC(integrand, a, b);
    return (b - a) / I;
}

template <typename ParticleType, typename BrushType, typename SurfaceInteractionModel>
template <typename MobilityFunc>
double ParticleBrushInteractionEnergy<ParticleType, BrushType, SurfaceInteractionModel>::
    diffusion_coefficient(const MobilityFunc mobility_phi, const double a, const double b) const
{
    auto integrand = [this, mobility_phi](const double z)
    {
        return std::exp(total_free_energy(z)) / mobility_factor(mobility_phi, z);
    };
    double I = INTEGRATE_FUNC(integrand, a, b);
    return (b - a) / I;
}

template <typename ParticleType, typename BrushType, typename SurfaceInteractionModel>
double ParticleBrushInteractionEnergy<ParticleType, BrushType, SurfaceInteractionModel>::
    partition_coefficient(const double a, const double b) const
{
    auto integrand = [this](const double z)
    { return std::exp(-total_free_energy(z)); };
    double I = INTEGRATE_FUNC(integrand, a, b);
    return I / b - a;
}

template <typename ParticleType, typename BrushType, typename SurfaceInteractionModel>
template <typename MobilityFunc>
double ParticleBrushInteractionEnergy<ParticleType, BrushType, SurfaceInteractionModel>::
    particle_concentration(const MobilityFunc mobility_phi, const double source_dist, const double z0) const
{
    if (z0 >= source_dist)
    {
        return 1.0;
    }
    auto integrand = [this, mobility_phi](const double z)
    {
        return std::exp(total_free_energy(z)) / mobility_factor(mobility_phi, z);
    };
    const double psi = INTEGRATE_FUNC(integrand, 0.0, z0)/INTEGRATE_FUNC(integrand, 0.0, source_dist);
    return psi / std::exp(total_free_energy(z0));
}

template <typename ParticleType, typename BrushType, typename SurfaceInteractionModel>
template <typename MobilityFunc>
double ParticleBrushInteractionEnergy<ParticleType, BrushType, SurfaceInteractionModel>::
    diffusion_coefficient_no_energy(MobilityFunc mobility_phi) const
{
    const double a = 0;
    const double b = brush->D() - particle->height;
    auto integrand = [this, mobility_phi](const double z)
    {
        return 1.0 / mobility_factor(mobility_phi, z);
    };
    return (b - a) / INTEGRATE_FUNC(integrand, a, b);
}