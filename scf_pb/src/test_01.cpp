#include <iostream>
#include <vector>
#include <numeric>
#include <chrono>
#include "ss_scf_common.hpp"
#include "normalization_condition.hpp"
#include "phi.hpp"

int main(int argc, char* argv[]){
    double chi = 0.0;
    double kappa = 0.0015707963267948967;
    double N=1000.0;
    double sigma = 0.02;
    //auto norm = open_planar(chi, N*sigma, kappa);
    auto norm = pore::open(chi, N*sigma*2*M_PI*300, kappa, 300.0);
    double d = solve_normalization(norm, N*sigma, N);
    //double R = d-10.0;
    int zsize = 100;
    std::vector<double> z(zsize);
    std::iota(z.begin(), z.end(), 0);

    for (unsigned i =0; i<zsize; i++){
        z[i] = d/(zsize-1)*i;
    }

    std::vector<double> phi(zsize);
    double phi_D = phi_at_zero_Pi(chi);

    auto start = std::chrono::high_resolution_clock::now();

    auto phi_z_func = [=](double z){return phi_z(chi, phi_D, d, z, kappa);};
    for (unsigned i =0; i<zsize; i++){
        phi[i] = phi_z_func(z[i]);
        std::cout << phi[i] << std::endl;
    }
    auto norm_r = pore::R_opening(N*sigma*2*M_PI*272, kappa, 0.0);
    double R_open = solve_normalization(norm_r, 50.0, 300.0);
    //std::cout << phi_D_r;
    //auto stop = std::chrono::high_resolution_clock::now();
    //auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    //start = stop;
    //std::cout << duration.count() << std::endl;
    //for (unsigned i = 0; i < 100; i++){
    //    D_unrestricted(chi, 20.0, kappa, 100.0);
    //}
    //stop = std::chrono::high_resolution_clock::now();
    //duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    //std::cout << duration.count() << std::endl;
    std::cout << R_open << std::endl;
    std::cout << d << std::endl;
}