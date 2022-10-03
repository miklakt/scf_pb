#pragma once
#include "ss_scf_common.hpp"
#include <vector>
#include <algorithm>
#include <iterator>
#include "normalization_condition.hpp"
#include <boost/math/quadrature/trapezoidal.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/interpolators/cubic_b_spline.hpp>

#define UNBOUND_MEMBER_FUNCTION(x) [this](double arg) { return this->x(arg); }
#define INTEGRATE_FUNC(F, a, b) boost::math::quadrature::gauss_kronrod<double, 31>::integrate(F, a, b)
#define MAKE_DEFAULT_CUMULATIVE_MEMBER_FUNCTION(F) \
    double F##_cumulative(const double z) const { return INTEGRATE_FUNC(UNBOUND_MEMBER_FUNCTION(F), 0.0, z); }

class BrushProfile
{
public:
    virtual double D() const = 0;
    virtual double phi_D() const = 0;

    virtual double phi_z(const double z) const = 0;
    virtual MAKE_DEFAULT_CUMULATIVE_MEMBER_FUNCTION(phi_z)

        virtual double Pi_z(const double z) const = 0;
    virtual MAKE_DEFAULT_CUMULATIVE_MEMBER_FUNCTION(Pi_z)
};

class BrushProfilePlanar : public BrushProfile
{
private:
    const double m_chi, m_kappa, m_R, m_N, m_sigma, m_theta;
    mutable double m_D, m_phi_D;

public:
    BrushProfilePlanar(double const &chi, double const &N, double const &sigma, double const &kappa, double const &R);
    BrushProfilePlanar(double const &chi, double const &N, double const &sigma, double const &kappa) : BrushProfilePlanar(chi, N, sigma, kappa, std::numeric_limits<double>::max()){};
    double phi_z(const double z) const override;
    double Pi_z(const double z) const override;
    double D() const override { return m_D; }
    double phi_D() const override { return m_phi_D; }
};

class BrushProfilePore : public BrushProfile
{
private:
    const double m_chi, m_kappa, m_R, m_N, m_sigma, m_theta;
    mutable double m_D, m_phi_D;

public:
    BrushProfilePore(double const &chi, double const &N, double const &sigma, double const &kappa, double const &R);
    double phi_z(const double z) const override;
    double Pi_z(const double z) const override;
    double D() const override { return m_D; }
    double phi_D() const override { return m_phi_D; }
};

class BrushProfileExternal : public BrushProfile
{
private:
    const std::vector<double> m_phi_profile;
    mutable double m_D, m_phi_D;
    const double m_dz, m_phi_cut, m_chi;
    const boost::math::cubic_b_spline<double> m_phi_spline;
public:
    BrushProfileExternal(const std::vector<double> &phi, const double chi, const double dz = 1, const double phi_cut = 1e-4);
    double phi_z(const double z) const override;
    double Pi_z(const double z) const override;
    double D() const override { return m_D; }
    double phi_D() const override { return m_phi_D; }
};