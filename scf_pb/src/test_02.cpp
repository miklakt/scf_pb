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
    //auto norm = pore::free::D(chi, N*sigma*2*M_PI*300, kappa, 300.0);
    //double d = solve_normalization(norm, N*sigma, N);
    //double R = d-10.0;
    BrushEdgePlanar brush(chi, N*sigma, kappa);
}