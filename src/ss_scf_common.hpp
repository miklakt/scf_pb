#pragma once
#include <cmath>
#include <boost/math/tools/roots.hpp>
#include <boost/math/special_functions/lambert_w.hpp>
#include <boost/math/special_functions/pow.hpp>

#define DEFAULT_PHI_MAX 0.9999f
#define GUESS_MINUS_MIN 0.1f
//#define USE_PREVIOUS_GUESS false
#define DEFAULT_MU_MAX(CHI) -std::log(1-DEFAULT_PHI_MAX) - 2*CHI*DEFAULT_PHI_MAX

namespace ss_scf_common{

//Osmotic pressure
template <class T>
T Pi(T const &phi, T const &chi){
    T pi = -std::log(1-phi)-chi*phi*phi-phi;
    return pi;
}

//phi_D from the vanishing pressure condition
template <class T>
T phi_at_zero_Pi(T const &chi, T const &min, T const &max, T const &guess)
{
    if (chi <= 0.5f){
        T result = static_cast<T>(0);
        return result;
    }

    auto f = [=](T phi){return Pi(phi, chi);};
    auto fprime  = [=](T phi){return - 2*chi*phi - 1 + 1/(1-phi);};
    auto fprime2  = [=](T phi){return - 2*chi + 1/((1-phi)*(1-phi));};
    auto solve_expr = [=](T phi){return std::make_tuple(f(phi), fprime(phi), fprime2(phi));};

    using namespace boost::math::tools;

    const int digits = std::numeric_limits<T>::digits;  // Maximum possible binary digits accuracy for type T.
    // digits used to control how accurate to try to make the result.
    int get_digits = static_cast<int>(digits);          // Accuracy triples with each step, so stop when just
                                                        // over one third of the digits are correct.
    std::uintmax_t maxit = 50;

    T result = halley_iterate(solve_expr, guess, min, max, get_digits, maxit);

    return result;
}

template <class T>
T phi_at_zero_Pi(T const &chi, T const &guess, T const &max = DEFAULT_PHI_MAX)
{
    if (chi <= 0.5f){
        T result = static_cast<T>(0);
        return result;
    }
    const T min = (chi-0.5)/chi;
    return phi_at_zero_Pi(chi, min, max, guess);
}

template <class T>
T phi_at_zero_Pi(T const &chi)
{
    if (chi <= 0.5f){
        T result = static_cast<T>(0);
        return result;
    }
    const T min = (chi-0.5)/chi;
    const T max = DEFAULT_PHI_MAX;
    //if (USE_PREVIOUS_GUESS){
    //    static T guess = (min + max)/2;
    //    auto res =  phi_at_zero_Pi(chi, min, max, guess);
    //    guess = res;
    //    return res;
    //}
    //else{
        T guess = (min + max)/2;
        return phi_at_zero_Pi(chi, min, max, guess);
    //}
}


//Segment potential at given z
template <class T>
T mu_z(T const &mu_D, T const &d, T const &z, T const &kappa){
    return (3/2)*kappa*kappa*(d*d-z*z) + mu_D;
}

//Segment potential for given phi
template <class T>
T mu_phi(T const &phi, T const &chi){
    return -std::log(1-phi) - 2*chi*phi;
}

//inverse of mu with respect to phi
template <class T>
T mu_inv(T const &mu, T const &chi){
    using boost::math::lambert_w0;
    if (chi==0){
        return 1-std::exp(-mu);
    }
    double lambert_arg = -2*chi*std::exp(-2*chi-mu);
    if (lambert_arg <= -0.367879){return 0.0;} //near singularity point
    return 1+lambert_w0(lambert_arg)/(2*chi);
}

//Polymer density at given z
template <class T>
T phi_z(T const &chi, T const &phi_D, T const &d, T const &z, T const &kappa){
    if (z>d){return T(0);}
    T mu_D = mu_phi(phi_D, chi);
    T mu = mu_z(mu_D, d, z, kappa);
    T phi = mu_inv(mu, chi);
    return phi;
}

//const double chi_crit = 6.0*std::log(5.0/6.0);
//double gamma_2poly_model(const double a1, const double a2, const double chi, const double chi_PC, const double phi){
//    double chi_ads = chi_PC - chi*(1-phi);
//    double psi = a1*phi + a2*phi*phi;
//    double gamma = (chi_ads-chi_crit)*psi;
//    return gamma;
//}
//
//const double mobility_factor(const double phi, const double d, const double k){
//    if (phi == 0){
//        return 1;
//    }
//    double eps = 1/phi;
//    double x = eps*eps/(d*d);
//    return x/std::pow((1+std::pow(x,k)), 1/k);
//}

}