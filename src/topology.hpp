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
    /* qq = 1+np.sum([q**i for i in range(1,g+1)])
    if g == 1:
        eta = math.atan(1/math.sqrt(q))*2/np.pi*qq
    if g == 2:
        eta =  math.atan(1/math.sqrt(q*(q+2)))*2/np.pi*qq
    return eta */
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


