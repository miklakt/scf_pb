#include "../ss_scf_common.hpp"
#include "../normalization_condition.hpp"
#include "../brush_profiles.hpp"
#include "../topology.hpp"
#include "../particle_geometry.hpp"
#include "../particle_energy.hpp"

#include <iostream>

int main(int argc, char* argv[]){
    double ph=10;
    double pw = 10;
    double chi = 0;
    double chi_PC = -2;
    double sigma = 0.02;
    double N = 1000;
    double a0 = 0.18;
    double a1 = -0.09;

    auto get_time = [](){return std::chrono::steady_clock::now();};
    auto const start_time = get_time();
    
    auto current_time = get_time();
    auto time_step = [&current_time, get_time](){
        auto delta = std::chrono::duration_cast<std::chrono::microseconds>(get_time() - current_time);
        current_time = get_time();
        return delta.count();
    };


    particle::Cylinder particle{ph,pw};
    BrushProfilePlanar brush{chi, N, sigma, topology::kappa(N)};

    auto gamma_phi = surface_interaction_coefficient::gamma_phi(a0, a1, chi, chi_PC);

    std::cout << "Gamma(phi=0.1): " << gamma_phi(0.1) << std::endl;

    ParticleBrushInteractionEnergy particle_in_brush{&particle, &brush, gamma_phi};

    double z = 10;
    std::cout << "Osmotic free energy: " << particle_in_brush.osmotic_free_energy(z) << std::endl;
    std::cout << "Surface free energy: " << particle_in_brush.surface_free_energy(z) << std::endl;
    std::cout << "Total free energy: " << particle_in_brush.total_free_energy(z) << std::endl;

    std::cout << "Apparent diffusion coefficient through the brush: " << particle_in_brush.diffusion_coefficient() << std::endl;

    std::cout << "Time passed:" << time_step() << std::endl;

}