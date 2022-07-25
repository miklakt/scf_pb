#pragma once
#include "ss_scf_common.hpp"
#include "normalization_condition.hpp"
#include <boost/math/quadrature/trapezoidal.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>

#define ALMOST_ONE 0.99999f

using namespace general_identities;

namespace planar
{
namespace free{
template <class T>
auto D(T const &chi, T const &theta, T const &kappa){
    using namespace std;
    using namespace boost::math::tools;
    using namespace boost::math::quadrature;

    T phi_D = phi_at_zero_Pi(chi);
    auto integrand = [chi, phi_D, kappa](T d){return [=](T z){return phi_z(chi, phi_D, d, z, kappa);};};
    auto integral = [integrand, theta](T d){return gauss_kronrod<double, 31>::integrate(integrand(d), T(0), d)- theta;};
    return integral;
}
} // namespace free

namespace restricted{
template <class T>
auto phi_D(T const &chi, T const &theta, T const &kappa, T const &R){
    using namespace std;
    using namespace boost::math::tools;
    using namespace boost::math::quadrature;

    auto integrand = [chi, R, kappa](T phi_D){return [=](T z){return phi_z(chi, phi_D, R, z, kappa);};};
    auto integral = [integrand, theta, R](T phi_D){return gauss_kronrod<double, 31>::integrate(integrand(phi_D), T(0), R)- theta;};
    return integral;
}
} // namespace restricted
} // namespace planar

namespace pore{
namespace free{
template <class T>
auto D(T const &chi, T const &theta, T const &kappa, T const &pore_R){
    using namespace std;
    using namespace boost::math::tools;
    using namespace boost::math::quadrature;

    T phi_D = phi_at_zero_Pi(chi);
    auto integrand = [chi, phi_D, kappa, pore_R](T d){return [=](T z){return phi_z(chi, phi_D, d, z, kappa)*abs(pore_R-z)*2*M_PI;};};
    auto integral = [integrand, theta](T d){return gauss_kronrod<double, 31>::integrate(integrand(d), T(0), d)- theta;};
    return integral;
}
} // namespace free

namespace restricted{
template <class T>
auto phi_D(T const &chi, T const &theta, T const &kappa, T const &R){
    using namespace std;
    using namespace boost::math::tools;
    using namespace boost::math::quadrature;

    const T pore_R = R;
    auto integrand = [chi, R, kappa, pore_R](T phi_D){return [=](T z){return phi_z(chi, phi_D, R, z, kappa)*abs(pore_R-z);};};
    auto integral = [integrand, theta, R](T phi_D){return 2*M_PI*gauss_kronrod<double, 31>::integrate(integrand(phi_D), T(0), R)- theta;};
    return integral;
}

template <class T>
auto chi_opening(T const &theta, T const &kappa, T const &R){
    using namespace std;
    using namespace boost::math::tools;
    using namespace boost::math::quadrature;

    const T pore_R = R;
    auto integrand = [R, kappa, pore_R](T chi){
        return [=](T z){
            return phi_z(chi, phi_at_zero_Pi(chi), R, z, kappa)*abs(pore_R-z);
        };
    };

    auto integral = [integrand, theta, R](T chi){
        return 2*M_PI*gauss_kronrod<double, 31>::integrate(integrand(chi), T(0), R) - theta;
        };
    return integral;
}

template <class T>
auto R_opening(T const &theta, T const &kappa, T const &chi){
    using namespace std;
    using namespace boost::math::tools;
    using namespace boost::math::quadrature;

    auto integrand = [kappa, chi](T R){
        return [=](T z){
            return phi_z(chi, phi_at_zero_Pi(chi), R, z, kappa)*abs(R-z);
        };
    };

    auto integral = [integrand, theta, chi](T R){
        return 2*M_PI*gauss_kronrod<double, 31>::integrate(integrand(R), T(0), R) - theta;
        };
    return integral;
}
} // namespace restricted
} // namespace pore

template <class F, class T>
T solve_normalization(F integral, T const &min, T const &max, const size_t maxit = 50)
{
  using namespace std;
  using namespace boost::math::tools;
  using namespace boost::math::quadrature;

  std::uintmax_t it = maxit;
  int digits = std::numeric_limits<T>::digits;

  int get_digits = digits - 3;
  eps_tolerance<T> tol(get_digits);

  std::pair<T, T> r = toms748_solve(integral, min, max, tol, it);

  return r.first + (r.second - r.first)/2;
}