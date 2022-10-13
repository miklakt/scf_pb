
template <typename ParticleType, typename BrushType, typename SurfaceInteractionModel>
double ParticleBrushInteractionEnergy<ParticleType, BrushType, SurfaceInteractionModel>::osmotic_free_energy(const double particle_center)
{
    auto integrand = [this, &particle_center](const double z)
    {
        const double z0 = particle_center - particle->height / 2;
        return brush->Pi_z(z0 + z) * particle->volume_integrand(z);
    };
    return INTEGRATE_FUNC(integrand, 0, particle->height);
}

template <typename ParticleType, typename BrushType, typename SurfaceInteractionModel>
double ParticleBrushInteractionEnergy<ParticleType, BrushType, SurfaceInteractionModel>::surface_free_energy(const double particle_center)
{
    const double A0 = particle->surface_edges()[0];
    const double A1 = particle->surface_edges()[1];

    const double z0 = particle_center - particle->height / 2;
    const double z1 = particle_center + particle->height / 2;

    const double phi0 = brush->phi_z(z0);
    const double phi1 = brush->phi_z(z1);

    const double gamma0 = this->gamma_phi(phi0);
    const double gamma1 = this->gamma_phi(phi1);

    const double free_energy_at_boundaries = A0 * phi0 * gamma0 + A1 * phi1 * gamma1;

    auto integrand = [this, &z0](const double arg)
    {
        const double phi = brush->phi_z(arg + z0);
        const double gamma = this->gamma_phi(phi);
        return phi * gamma * particle->surface_integrand(arg);
    };

    return INTEGRATE_FUNC(integrand, 0.0, particle->height) + free_energy_at_boundaries;
}

template <typename ParticleType, typename BrushType, typename SurfaceInteractionModel>
double ParticleBrushInteractionEnergy<ParticleType, BrushType, SurfaceInteractionModel>::diffusion_coefficient(const double a, const double b)
{
    auto integrand = [this](const double z)
    { return std::exp(total_free_energy(z)); };
    return (b - a) / INTEGRATE_FUNC(integrand, a, b);
}

template <typename ParticleType, typename BrushType, typename SurfaceInteractionModel>
double ParticleBrushInteractionEnergy<ParticleType, BrushType, SurfaceInteractionModel>::diffusion_coefficient()
{
    const double a = particle->height / 2;
    const double b = brush->D() + particle->height / 2;
    return diffusion_coefficient(a, b);
}

template <typename ParticleType, typename BrushType, typename SurfaceInteractionModel>
template <typename MobilityFunc>
double ParticleBrushInteractionEnergy<ParticleType, BrushType, SurfaceInteractionModel>::diffusion_coefficient(const MobilityFunc mobility_phi, const double a, const double b)
{
    auto integrand = [this, mobility_phi](const double z)
    {
        const double phi = brush->phi_z(z);
        const double mobility_factor = mobility_phi(phi);
        return std::exp(total_free_energy(z)) / mobility_factor;
    };
    return (b - a) / INTEGRATE_FUNC(integrand, a, b);
}

template <typename ParticleType, typename BrushType, typename SurfaceInteractionModel>
template <typename MobilityFunc>
double ParticleBrushInteractionEnergy<ParticleType, BrushType, SurfaceInteractionModel>::diffusion_coefficient(const MobilityFunc mobility_phi)
{
    const double a = particle->height / 2;
    const double b = brush->D() + particle->height / 2;
    return diffusion_coefficient(mobility_phi, a, b);
}

template <typename ParticleType, typename BrushType, typename SurfaceInteractionModel>
double ParticleBrushInteractionEnergy<ParticleType, BrushType, SurfaceInteractionModel>::partition_coefficient(const double a, const double b)
{
    auto integrand = [this](const double z)
    { return std::exp(-total_free_energy(z)); };
    const double left = a - particle->height / 2;
    const double right = b - particle->height / 2;
    return INTEGRATE_FUNC(integrand, a, b) / (b - a);
}

template <typename ParticleType, typename BrushType, typename SurfaceInteractionModel>
double ParticleBrushInteractionEnergy<ParticleType, BrushType, SurfaceInteractionModel>::partition_coefficient(const double b)
{
    const double a = 0;
    return partition_coefficient(a, b);
}

template <typename ParticleType, typename BrushType, typename SurfaceInteractionModel>
double ParticleBrushInteractionEnergy<ParticleType, BrushType, SurfaceInteractionModel>::partition_coefficient()
{
    const double a = 0;
    const double b = brush->D() - particle->height/2;
    return partition_coefficient(a, b);
}

template <typename ParticleType, typename BrushType, typename SurfaceInteractionModel>
double ParticleBrushInteractionEnergy<ParticleType, BrushType, SurfaceInteractionModel>::partition_coefficient_open()
{
    const double a = 0;
    const double b = brush->D() + particle->height/2;
    return partition_coefficient(a, b);
}

template <typename ParticleType, typename BrushType, typename SurfaceInteractionModel>
double ParticleBrushInteractionEnergy<ParticleType, BrushType, SurfaceInteractionModel>::partition_coefficient_slit(const double R)
{
    const double a = 0;
    const double b = std::min(R, brush->D() + particle->height/2);
    return partition_coefficient(a, b);
}