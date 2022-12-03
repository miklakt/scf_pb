#include "../particle_geometry.hpp"

#define INTEGRATE_FUNC(F, a, b) boost::math::quadrature::gauss_kronrod<double, 31>::integrate(F, a, b)

bool test_volume_integral(const Particle* particle){
    auto volume_integrated = INTEGRATE_FUNC([particle](const double z){return particle->volume_integrand(z);}, 0, particle->height);
    auto surf_integrated = INTEGRATE_FUNC([particle](const double z){return particle->surface_integrand(z);}, 0, particle->height);
    surf_integrated += particle->surface_edges()[0] + particle->surface_edges()[1];

    std::cout << "Height:" << particle->height <<"\n";
    std::cout << "Width:" << particle->width <<"\n";
    std::cout << "Integrated volume:" << volume_integrated <<"\n";
    std::cout << "Volume" << particle->volume() <<"\n";

    std::cout << "Integrated surface:" << surf_integrated <<"\n";
    std::cout << "Surface" << particle->surface() <<"\n";

    return ((volume_integrated-particle->volume())<1e-6);
}

int main(int argc, char* argv[]){
    const double r = 20;
    particle::Cylinder particle{r, r};
    test_volume_integral(&particle);

    /* std::cout << "_______________________" <<"\n";
    particle::Sphere particle2{r};
    test_volume_integral(&particle2); */

    
    auto volume = particle::integrators::integrate_over_volume(&particle, [f = 1.0f](const double z){return f;});
    std::cout << "Integrated volume (integrator):" << volume <<"\n";
    auto surface = particle::integrators::integrate_over_surface(&particle, [f = 1.0f](const double z){return f;});
    std::cout << "Integrated surface (integrator):" << surface <<"\n";
}