template <typename ParticleType, typename BrushType, typename SurfaceInteractionModel>
double ParticleBrushInteractionEnergy<ParticleType, BrushType, SurfaceInteractionModel>::osmotic_free_energy(const double z0) const
{
    //particle is outside the brush
    if (z0 > brush->D()){
        return 0.0;
    }

    const double right_limit = std::min(particle->height, brush->D()-z0);

    return particle::integrators::integrate_over_volume(particle, [brush=brush, z0](const double z){return brush->Pi_z(z0 + z);}, right_limit);
}

template <typename ParticleType, typename BrushType, typename SurfaceInteractionModel>
double ParticleBrushInteractionEnergy<ParticleType, BrushType, SurfaceInteractionModel>::surface_free_energy(const double z0) const
{
    //particle is outside the brush
    if (z0 > brush->D()){
        return 0.0;
    }

    const double right_limit = std::min(particle->height, brush->D()-z0);

    auto integrand = [this, &z0](const double arg)
    {
        const double phi = brush->phi_z(arg + z0);
        const double gamma = this->gamma_phi(phi);
        return phi * gamma;
    };

    return particle::integrators::integrate_over_surface(particle, integrand, right_limit);
}

template <typename ParticleType, typename BrushType, typename SurfaceInteractionModel>
template <typename MobilityFunc>
double ParticleBrushInteractionEnergy<ParticleType, BrushType, SurfaceInteractionModel>::mobility_factor(const MobilityFunc mobility_phi, const double z0) const
{
    const double phi = brush->phi_z(z0);
    return mobility_phi(phi);
}

template <typename ParticleType, typename BrushType, typename SurfaceInteractionModel>
double ParticleBrushInteractionEnergy<ParticleType, BrushType, SurfaceInteractionModel>::diffusion_coefficient(const double a, const double b) const
{
    if (a>=b){throw std::invalid_argument( "a>=b" );}
    auto integrand = [this](const double z)
    { return std::exp(total_free_energy(z));};
    return (b - a) / brush::integrators::integrate_over_z(brush, particle, integrand, a, b);
}

template <typename ParticleType, typename BrushType, typename SurfaceInteractionModel>
double ParticleBrushInteractionEnergy<ParticleType, BrushType, SurfaceInteractionModel>::diffusion_coefficient() const
{
    const double a = 0.0;
    const double b = brush->D() - particle->height;
    return diffusion_coefficient(a, b);
}

template <typename ParticleType, typename BrushType, typename SurfaceInteractionModel>
double ParticleBrushInteractionEnergy<ParticleType, BrushType, SurfaceInteractionModel>::diffusion_coefficient_open() const
{
    const double a = 0.0;
    const double b = brush->D();
    return diffusion_coefficient(a, b);
}

template <typename ParticleType, typename BrushType, typename SurfaceInteractionModel>
double ParticleBrushInteractionEnergy<ParticleType, BrushType, SurfaceInteractionModel>::diffusion_coefficient_slit(const double R) const
{
    const double a = 0.0;
    const double b = std::min(R, brush->D());
    return diffusion_coefficient(a, b);
}

template <typename ParticleType, typename BrushType, typename SurfaceInteractionModel>
template <typename MobilityFunc>
double ParticleBrushInteractionEnergy<ParticleType, BrushType, SurfaceInteractionModel>::diffusion_coefficient(const MobilityFunc mobility_phi, const double a, const double b) const
{
    auto integrand = [this, mobility_phi](const double z)
    {
        //const double phi = brush->phi_z(z);
        //const double mobility_factor = mobility_phi(phi);
        return std::exp(total_free_energy(z)) / mobility_factor(mobility_phi, z);
    };
    //return (b - a) / INTEGRATE_FUNC(integrand, a, b);
    return (b - a) / brush::integrators::integrate_over_z(brush, particle, integrand, a, b);
}



template <typename ParticleType, typename BrushType, typename SurfaceInteractionModel>
template <typename MobilityFunc>
double ParticleBrushInteractionEnergy<ParticleType, BrushType, SurfaceInteractionModel>::diffusion_coefficient(const MobilityFunc mobility_phi) const
{
    const double a = 0;
    const double b = brush->D() - particle->height;
    return diffusion_coefficient(mobility_phi, a, b);
}

template <typename ParticleType, typename BrushType, typename SurfaceInteractionModel>
template <typename MobilityFunc>
double ParticleBrushInteractionEnergy<ParticleType, BrushType, SurfaceInteractionModel>::diffusion_coefficient_open(const MobilityFunc mobility_phi) const
{
    const double a = 0;
    const double b = brush->D();
    return diffusion_coefficient(mobility_phi, a, b);
}

template <typename ParticleType, typename BrushType, typename SurfaceInteractionModel>
template <typename MobilityFunc>
double ParticleBrushInteractionEnergy<ParticleType, BrushType, SurfaceInteractionModel>::diffusion_coefficient_slit(const MobilityFunc mobility_phi, const double R) const
{
    const double a = 0;
    const double b = std::min(R, brush->D());
    return diffusion_coefficient(mobility_phi, a, b);
}

template <typename ParticleType, typename BrushType, typename SurfaceInteractionModel>
double ParticleBrushInteractionEnergy<ParticleType, BrushType, SurfaceInteractionModel>::partition_coefficient(const double a, const double b) const
{
    auto integrand = [this](const double z)
    { return std::exp(-total_free_energy(z)); };
    return brush::integrators::integrate_over_z(brush, particle, integrand, a, b)/(b - a);
}

template <typename ParticleType, typename BrushType, typename SurfaceInteractionModel>
double ParticleBrushInteractionEnergy<ParticleType, BrushType, SurfaceInteractionModel>::partition_coefficient() const
{
    const double a = 0;
    const double b = brush->D() - particle->height;
    return partition_coefficient(a, b);
}

template <typename ParticleType, typename BrushType, typename SurfaceInteractionModel>
double ParticleBrushInteractionEnergy<ParticleType, BrushType, SurfaceInteractionModel>::partition_coefficient_open() const
{
    const double a = 0;
    const double b = brush->D();
    return partition_coefficient(a, b);
}

template <typename ParticleType, typename BrushType, typename SurfaceInteractionModel>
double ParticleBrushInteractionEnergy<ParticleType, BrushType, SurfaceInteractionModel>::partition_coefficient_slit(const double R) const
{
    const double a = 0;
    const double b = std::min(R, brush->D());
    return partition_coefficient(a, b);
}


template <typename ParticleType, typename BrushType, typename SurfaceInteractionModel>
template <typename MobilityFunc>
double ParticleBrushInteractionEnergy<ParticleType, BrushType, SurfaceInteractionModel>::particle_concentration(const MobilityFunc mobility_phi, const double source_dist, const double z0, const double c_bulk) const
{   
    if (z0>=source_dist){return c_bulk;}
    auto integrand = [this, mobility_phi](const double z)
    {
        return std::exp(total_free_energy(z)) / mobility_factor(mobility_phi, z);
    };
    const double psi = c_bulk * brush::integrators::integrate_over_z(brush, particle, integrand, 0.0, z0)/brush::integrators::integrate_over_z(brush, particle, integrand, 0.0, source_dist);
    return psi/std::exp(total_free_energy(z0));
}


template <typename ParticleType, typename BrushType, typename SurfaceInteractionModel>
template <typename MobilityFunc>
double ParticleBrushInteractionEnergy<ParticleType, BrushType, SurfaceInteractionModel>::partition_coefficient_perfect_sink(const MobilityFunc mobility_phi, const double b, const double source_dist, const double c_bulk) const
{
    auto integrand = [=](const double z)
    {
        return particle_concentration(mobility_phi, source_dist, z, c_bulk);
    };
    return INTEGRATE_FUNC(integrand, 0.0, b)/b;
}

template <typename ParticleType, typename BrushType, typename SurfaceInteractionModel>
template <typename MobilityFunc>
double ParticleBrushInteractionEnergy<ParticleType, BrushType, SurfaceInteractionModel>::partition_coefficient_perfect_sink(const MobilityFunc mobility_phi, const double source_dist, const double c_bulk) const
{
    const double b = brush->D() - particle->height;
    return partition_coefficient_perfect_sink(mobility_phi, b, source_dist, c_bulk);
}

template <typename ParticleType, typename BrushType, typename SurfaceInteractionModel>
template <typename MobilityFunc>
double ParticleBrushInteractionEnergy<ParticleType, BrushType, SurfaceInteractionModel>::partition_coefficient_perfect_sink_open(const MobilityFunc mobility_phi, const double source_dist, const double c_bulk) const
{
    const double b = brush->D();
    return partition_coefficient_perfect_sink(mobility_phi, b, source_dist, c_bulk);
}


template <typename ParticleType, typename BrushType, typename SurfaceInteractionModel>
template <typename MobilityFunc>
double ParticleBrushInteractionEnergy<ParticleType, BrushType, SurfaceInteractionModel>::diffusion_coefficient_no_energy(MobilityFunc mobility_phi) const
{   
    const double a = 0;
    const double b = brush->D() - particle->height;
    auto integrand = [this, mobility_phi](const double z)
    {
        return 1.0 / mobility_factor(mobility_phi, z);
    };
    return (b - a) / brush::integrators::integrate_over_z(brush, particle, integrand, a, b);
}