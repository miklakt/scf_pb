#include <iostream>
#include <chrono>
#include <memory>
#include <cmath>
#include "../ss_scf_common.hpp"
#include "../normalization_condition.hpp"
#include "../profiles.hpp"
#include "../topology.hpp"
#include "../particles.hpp"
#include "../energy.hpp"



int main(int argc, char* argv[]){

    double chi = 0.0;
    double N=1000.0;
    double sigma = 0.02;
    double R = 300;
    double kappa = topology::kappa(N);

    auto start = std::chrono::high_resolution_clock::now();

    //BrushProfilePore brush(chi, N, sigma, kappa, R);

    BrushProfile *b_old = new BrushProfilePore(chi, N, sigma, kappa, R);

    //std::cout << "phi profile:";
    std::vector<double>phi_external;
    for (int i=0; i < b_old->D()+1; ++i){
        phi_external.push_back(b_old->phi_z(i));
        //std::cout << " " << b->phi_z(i) << std::endl;
    }
    std::cout << std::endl;

    auto first_smaller_it = std::lower_bound(phi_external.rbegin(), phi_external.rend(), 1e-4);
    double idx = std::distance(first_smaller_it, phi_external.rend());
    std::cout << "interpolation test D: " << idx << std::endl;

    BrushProfile *b = new BrushProfileExternal(phi_external, chi);

/*     for (int i=0; i < b->D()+1; ++i ){
        std::cout << " " << b->phi_z(i) << std::endl;
    } */

    std::cout << "Program tests brush density profile calculation in a pore " << std::endl;

    std::cout << "------------------------------------------" << std::endl;

    std::cout << "chi: " << chi << std::endl;
    std::cout << "N: " << N << std::endl;
    std::cout << "sigma: " << sigma << std::endl;
    std::cout << "kappa: " << kappa << std::endl;
    std::cout << "R: " << R << std::endl;

    std::cout << "------------------------------------------" << std::endl;

    std::cout << "H: " << b->D() << std::endl;
    std::cout << "phi_D: " << b->phi_D() << std::endl;
    std::cout << "phi_0: " << b->phi_z(0) << std::endl;

    std::cout << "------------------------------------------" << std::endl;

    std::cout << "Pi_D: " << b->Pi_z(b->D()) << std::endl;
    std::cout << "Pi_0: " << b->Pi_z(0) << std::endl;

    std::cout << "------------------------------------------" << std::endl;


    //double phi_integral = integrator([b](double z){return b->phi_z(z);}, 0.0, b->D());

    std::cout << "theta: " << b->phi_z_cumulative(b->D()) << std::endl;

    std::cout << "------------------------------------------" << std::endl;

    //init particle
    Particle *particle = new particle::Cylinder(4.0, 4.0);

    std::cout << "particle_surface_edges: " << particle->surface_edges()[0] <<"; "<<particle->surface_edges()[1]<< std::endl;
    std::cout << "particle_height: " << particle->width << std::endl;
    std::cout << "particle_width: " << particle->height << std::endl;
    std::cout << "particle_surface: " << particle->surface() << std::endl;
    std::cout << "particle_volume: " << particle->volume() << std::endl;

    std::cout << "------------------------------------------" << std::endl;

    double position = 10;
    double osmotic_free_energy = make_function::osmotic_free_energy_func(b, particle)(position);
    double a0 = 0.18; double a1 = -0.09;
    double chi_PC = -1.5;
    auto gamma_phi = make_function::gamma_phi(a0, a1, chi, chi_PC);
    auto surface_free_energy = make_function::surface_free_energy_func(b, particle, gamma_phi);
    auto total_free_energy = make_function::total_free_energy_func(b, particle, gamma_phi);

    std::cout << "osmotic energy(z=10): " << osmotic_free_energy << std::endl;
    std::cout << "surface energy(z=10): " << surface_free_energy(position) << std::endl;
    std::cout << "total energy(z=10): " << total_free_energy(position) << std::endl;

    std::cout << "------------------------------------------" << std::endl;

    //double d_eff = make_function::effective_diffusion_coefficient_ab(b, particle, gamma_phi)(particle->height/2, b->D());
    double d_eff = effective_diffusion_coefficient(b, particle, gamma_phi);

    std::cout << "D_eff: " << d_eff << std::endl;

    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start);
    std::cout << "wall time: " << duration.count() << std::endl;

}