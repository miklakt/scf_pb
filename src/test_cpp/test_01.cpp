#include <iostream>
#include <vector>
#include <numeric>
#include <chrono>
#include "../ss_scf_common.hpp"
#include "../normalization_condition.hpp"
#include "../brush_profiles.hpp"
#include "../topology.hpp"

int main(int argc, char* argv[]){

    std::cout << "Free brush-----------------------------------------------------------------" << std::endl;

    double chi = 0.6;
    double N=1000.0;
    double sigma = 0.02;
    double kappa = topology::kappa(N);

    auto get_time = [](){return std::chrono::steady_clock::now();};
    auto const start_time = get_time();

    auto current_time = get_time();
    auto time_step = [&current_time, get_time](){
        auto delta = std::chrono::duration_cast<std::chrono::microseconds>(get_time() - current_time);
        current_time = get_time();
        return delta.count();
    };

    BrushProfilePlanar brush(chi, N, sigma, kappa);
    std::cout << "Brush profile instantiated" << std::endl;
    std::cout << time_step() << std::endl;

    double d = brush.D();
    std::cout << "Brush thickness D:" << d << std::endl;
    std::cout << "Volume concentration at the edge phi_D:" << brush.phi_D() << std::endl;
    time_step();

    int zsize = 11;
    std::vector<double> z(zsize);
    std::iota(z.begin(), z.end(), 0);
    for (unsigned i =0; i<zsize; i++){
        z[i] = d/(zsize-1)*i;
    }
    std::vector<double> phi(zsize);
    std::transform(z.begin(), z.end(), phi.begin() ,[brush](const double z_){return brush.phi_z(z_);});
    auto t = time_step();
    std::cout << "Profile is calculated:" << std::endl;
    std::cout << "z/D" << "\t" << "phi" << std::endl;
    for (unsigned i =0; i<zsize; i++){
        std::cout << z[i]/d << "\t" << phi[i] << std::endl;
    }
    std::cout << "Time passed:" << t << std::endl;

    std::cout << "Restricted brush------------------------------------------------------------" << std::endl;
    time_step();
    std::cout << "R opening: " << solve_normalization(planar::restricted::R_opening(chi, N*sigma, kappa), N*sigma, N) << std::endl;
    std::cout << "Time passed:" << time_step() << std::endl;
    double R = 50;

    std::cout << "Brush restricted at R=" << R << std::endl;
    std::cout << "chi opening: " << solve_normalization(planar::restricted::chi_opening(N*sigma, kappa, R), 0.0, 1.0) << std::endl;
    std::cout << "Time passed:" << time_step() << std::endl;

    BrushProfilePlanar brush_r(chi, N, sigma, kappa, R);
    std::cout << "Brush profile instantiated" << std::endl;
    std::cout << time_step() << std::endl;

    double d_r = brush_r.D();
    std::cout << "Brush thickness D:" << d_r << std::endl;
    std::cout << "Volume concentration at the edge phi_D:" << brush_r.phi_D() << std::endl;
    time_step();

    int zsize_r = 11;
    std::vector<double> z_r(zsize);
    std::iota(z_r.begin(), z_r.end(), 0);
    for (unsigned i =0; i<zsize_r; i++){
        z_r[i] = d_r/(zsize_r-1)*i;
    }
    std::vector<double> phi_r(zsize_r);
    std::transform(z_r.begin(), z_r.end(), phi_r.begin() ,[brush_r](const double z_){return brush_r.phi_z(z_);});
    t = time_step();
    std::cout << "Profile is calculated:" << std::endl;
    std::cout << "z/D" << "\t" << "phi" << std::endl;
    for (unsigned i =0; i<zsize_r; i++){
        std::cout << z_r[i]/d_r << "\t" << phi_r[i] << std::endl;
    }
    std::cout << "Time passed:" << t << std::endl;




}