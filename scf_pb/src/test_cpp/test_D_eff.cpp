#include "../ss_scf_common.hpp"
#include "../normalization_condition.hpp"
#include "../profiles.hpp"
#include "../topology.hpp"
#include "../particles.hpp"
#include "../energy.hpp"

double D_eff(const double N, const double sigma, const double chi, const double chi_PC, const double a0, const double a1, const double particle_width, const double particle_height){
    double kappa = topology::kappa(N);
    BrushProfilePlanar brush(chi, N, sigma, kappa);
    particle::Cylinder particle(particle_height, particle_width);
    auto gamma_phi = make_function::gamma_phi(a0, a1, chi, chi_PC);
    double d_eff = effective_diffusion_coefficient(&brush, &particle, gamma_phi);
    return d_eff;
}

double D_eff_mobility(const double N, const double sigma, const double chi, const double chi_PC, const double a0, const double a1, const double particle_width, const double particle_height){
    double kappa = topology::kappa(N);
    BrushProfilePlanar brush(chi, N, sigma, kappa);
    particle::Cylinder particle(particle_height, particle_width);
    auto gamma_phi = make_function::gamma_phi(a0, a1, chi, chi_PC);
    auto mobility = make_function::mobility_phi(particle_width);
    double d_eff = effective_diffusion_coefficient(&brush, &particle, gamma_phi, mobility);
    return d_eff;
}

int main(int argc, char* argv[]){
    double N = 1000;
    double chi = 0.5;
    double sigma = 0.02;
    double chi_PC = -1.0;
    double a0 = 0.18;
    double a1 = -0.09;
    double width = 4.0;
    double height = 4.0;



    double d_eff = D_eff(N, sigma, chi, chi_PC, a0, a1, width, height);
    double d_eff_corrected = D_eff_mobility(N, sigma, chi, chi_PC, a0, a1, width, height);

    std::cout << "Effective diffusion coefficient: " << d_eff << std::endl;
    std::cout << "Corrected effective diffusion coefficient: " << d_eff_corrected << std::endl;
}
