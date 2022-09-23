#pragma once
#include "ss_scf_common.hpp"
#include "normalization_condition.hpp"
#include <boost/math/quadrature/trapezoidal.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>

#define ALMOST_ONE 0.99999f

using namespace ss_scf_common;
using namespace boost::math::tools;
using namespace boost::math::quadrature;


template <class F, class T>
T solve_normalization(F integral, T const &min, T const &max, const size_t maxit = 50)
{
  std::uintmax_t it = maxit;
  int digits = std::numeric_limits<T>::digits;

  int get_digits = digits - 3;
  eps_tolerance<T> tol(get_digits);

  std::pair<T, T> r = toms748_solve(integral, min, max, tol, it);

  return r.first + (r.second - r.first)/2;
}

namespace planar
{
template <class T>
T normalization_integral(T const &chi, T const &phi_D, T const& d,  T const &kappa, T const &theta){
    auto integrand = [=](T z){return phi_z(chi, phi_D, d, z, kappa);};
    return gauss_kronrod<T, 31>::integrate(integrand, T(0), d) - theta;
}

namespace free{
template <class T>
auto D(T const &chi, T const &theta, T const &kappa){
    T phi_D = phi_at_zero_Pi(chi);
    auto integral = [=](T const& d){return normalization_integral(chi, phi_D, d, kappa, theta);};
    return integral;
}
} // namespace free

namespace restricted{
template <class T>
auto phi_D(T const &chi, T const &theta, T const &kappa, T const &R){
    auto integral = [=](T const& phi_D){return normalization_integral(chi, phi_D, R, kappa, theta);};
    return integral;
}

template <class T>
auto chi_opening(T const &theta, T const &kappa, T const &R){
    auto integral = [=](T const& chi){
        T phi_D = phi_at_zero_Pi(chi);
        return normalization_integral(chi, phi_D, R, kappa, theta);
    };
    return integral;
}

template <class T>
auto R_opening(T const &chi, T const &theta, T const &kappa){
    T phi_D = phi_at_zero_Pi(chi);
    auto integral = [=](T const& R){return normalization_integral(chi, phi_D, R, kappa, theta);};
    return integral;
}
} // namespace restricted
} // namespace planar

namespace pore{

template <class T>
double normalization_integral(T const &chi, T const &phi_D, T const& d,  T const &kappa, T const &theta, T const &pore_R){
    auto integrand = [=](T z){return phi_z(chi, phi_D, d, z, kappa)*abs(pore_R-z)*2*M_PI;};
    return gauss_kronrod<double, 31>::integrate(integrand, T(0), d) - theta;
}

namespace free{
template <class T>
auto D(T const &chi, T const &theta, T const &kappa, T const &pore_R){
    T phi_D = phi_at_zero_Pi(chi);
    auto integral = [=](T const& d){return normalization_integral(chi, phi_D, d, kappa, theta, pore_R);};
    return integral;
}
} // namespace free

namespace restricted{
template <class T>
auto phi_D(T const &chi, T const &theta, T const &kappa, T const &R){
    const T pore_R = R;
    auto integral = [=](T const& phi_D){return normalization_integral(chi, phi_D, R, kappa, theta, pore_R);};
    return integral;
}

template <class T>
auto chi_opening(T const &theta, T const &kappa, T const &R){
    const T pore_R = R;
    auto integral = [=](T const& chi){
        T phi_D = phi_at_zero_Pi(chi);
        return normalization_integral(chi, phi_D, R, kappa, theta, pore_R);
        };
    return integral;
}

template <class T>
auto R_opening(T const &theta, T const &kappa, T const &chi){
    T phi_D = phi_at_zero_Pi(chi);
    auto integral = [=](T const& R){
        const T pore_R = R;
        return normalization_integral(chi, phi_D, R, kappa, theta, pore_R);
        };
    return integral;
}
} // namespace restricted
} // namespace pore