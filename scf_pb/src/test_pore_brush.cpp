#include <iostream>
#include "ss_scf_common.hpp"
#include "normalization_condition.hpp"
#include "profiles.hpp"

int main(int argc, char* argv[]){
    
    double chi = 0.6;
    double kappa = 0.0015707963267948967;
    double N=1000.0;
    double sigma = 0.02;
    double R = 300;

    BrushPhiProfilePore brush(chi, N, sigma, kappa, R);

    std::cout << "Program tests brush density profile calculation in a pore " << std::endl;

    std::cout << "------------------------------------------" << std::endl;

    std::cout << "chi: " << chi << std::endl;
    std::cout << "N: " << N << std::endl;
    std::cout << "sigma: " << sigma << std::endl;
    std::cout << "kappa: " << kappa << std::endl;
    std::cout << "R: " << R << std::endl;

    std::cout << "------------------------------------------" << std::endl;

    std::cout << "H: " << brush.D() << std::endl;
    std::cout << "phi_D: " << brush.phi_D() << std::endl;

    std::cout << "------------------------------------------" << std::endl;
    

}