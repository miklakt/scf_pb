#pragma once
#include <math.h>
#include <stdexcept>

namespace topology
{

template <class T>
T kappa(T const &N, T const &eta = T(1.0)){
    return M_PI*eta/(2*N);
}

template <class T>
T eta(int const &g, int const &q){
    int qq = 1;
    for (int i=1; i<g+1; ++i){
        qq+=q*i;
    }
    switch (g)
    {
    case 1:
        return std::atan(1/std::sqrt(q))*2*qq/M_PI;
    case 2:
        return std::atan(1/std::sqrt(q*(q+2)))*2*qq/M_PI;
    default:
        throw std::invalid_argument("g other than 1 and 2 are not implemented");
    }
}

} // namespace topology
