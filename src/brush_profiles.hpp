#pragma once
#include "ss_scf_common.hpp"
#include <vector>
#include <algorithm>
#include <iterator>
#include "normalization_condition.hpp"
#include <boost/math/quadrature/trapezoidal.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/interpolators/cubic_b_spline.hpp>

#define UNBOUND_MEMBER_FUNCTION(x) [this](double arg){return this->x(arg);}
#define INTEGRATE_FUNC(F, a, b) boost::math::quadrature::gauss_kronrod<double, 31>::integrate(F, a, b);
#define MAKE_DEFAULT_CUMULATIVE_MEMBER_FUNCTION(F) double F##_cumulative(const double z) const {return INTEGRATE_FUNC(UNBOUND_MEMBER_FUNCTION(F), 0.0, z);}

class BrushProfile
{
    public:
        virtual double D() const = 0;
        virtual double phi_D() const = 0;

        virtual double phi_z(const double z) const =0;
        virtual MAKE_DEFAULT_CUMULATIVE_MEMBER_FUNCTION(phi_z)

        virtual double Pi_z(const double z) const =0;
        virtual MAKE_DEFAULT_CUMULATIVE_MEMBER_FUNCTION(Pi_z)
};


class BrushProfilePlanar : public BrushProfile
{
    private:
    const double m_chi, m_kappa, m_R, m_N, m_sigma, m_theta;
    mutable double m_D, m_phi_D;

    public:
    BrushProfilePlanar(double const &chi, double const &N, double const &sigma, double const &kappa, double const &R)
    : m_chi(chi), m_N(N), m_sigma(sigma), m_kappa(kappa), m_R(R), m_theta(N*sigma)
    {
        double D_free = solve_normalization(planar::free::D(m_chi, m_theta, m_kappa), m_theta, m_N);
        //free brush
        if (D_free<m_R){
            m_D = D_free;
            m_phi_D = phi_at_zero_Pi(m_chi);
        }
        else{
            m_phi_D = solve_normalization(planar::restricted::phi_D(m_chi, m_theta, m_kappa, m_R), 0.0, double(ALMOST_ONE));
            m_D = m_R;
        }
    }
    BrushProfilePlanar(double const &chi, double const &N, double const &sigma, double const &kappa) : BrushProfilePlanar(chi, N, sigma, kappa, std::numeric_limits<double>::max()) {};

    double D() const override{
        return m_D;
    }

    double phi_D() const override{
        return m_phi_D;
    }

    double phi_z(const double z) const override{
        return ss_scf_common::phi_z(m_chi, m_phi_D, m_D, z, m_kappa);
    }

    double Pi_z(const double z) const override{
        return ss_scf_common::Pi(phi_z(z), m_chi);
    }
};

class BrushProfilePore : public BrushProfile
{
    private:
    const double m_chi, m_kappa, m_R, m_N, m_sigma, m_theta;
    mutable double m_D, m_phi_D;

    public:
    BrushProfilePore(double const &chi, double const &N, double const &sigma, double const &kappa, double const &R)
    : m_chi(chi), m_N(N), m_sigma(sigma), m_kappa(kappa), m_R(R), m_theta(2*M_PI*R*N*sigma)
    {
        double D_free = solve_normalization(pore::free::D(m_chi, m_theta, m_kappa, m_R), m_N*m_sigma, m_N);
        //free brush
        if (D_free<m_R){
            m_D = D_free;
            m_phi_D = phi_at_zero_Pi(m_chi);
        }
        else{
            m_phi_D = solve_normalization(pore::restricted::phi_D(m_chi, m_theta, m_kappa, m_R), 0.0, double(ALMOST_ONE));
            m_D = m_R;
        }
    }

    double D() const override{
        return m_D;
    }

    double phi_D() const override{
        return m_phi_D;
    }

    double phi_z(const double z) const override{
        return ss_scf_common::phi_z(m_chi, m_phi_D, m_D, z, m_kappa);
    }

    double Pi_z(const double z) const override{
        return ss_scf_common::Pi(phi_z(z), m_chi);
    }
};


class BrushProfileExternal : public BrushProfile{
    private:
        const std::vector<double> m_phi_profile;
        mutable double m_D, m_phi_D;
        const double m_dz, m_phi_cut, m_chi;
        const boost::math::cubic_b_spline<double> m_phi_spline;

    public:
    BrushProfileExternal(const std::vector<double> &phi, const double chi, const double dz=1, const double phi_cut = 1e-4)
    : m_phi_profile(phi), m_dz(dz), m_phi_cut(phi_cut), m_chi(chi), m_phi_spline(phi.begin(), phi.end(), 0, dz)
    {
        auto first_smaller_it = std::lower_bound(m_phi_profile.rbegin(), m_phi_profile.rend(), m_phi_cut);
        double idx = std::distance(first_smaller_it, m_phi_profile.rend())*m_dz;
        m_D = idx*dz;
        m_phi_D = m_phi_profile[idx];
        //m_D = boost::math::tools::bisect(m_phi_spline, 0, m_phi_profile.size()*dx, phi_cut);
    }

    double phi_z(const double z) const override{
        if (z>m_D)
            return 0.0;
        if (z<0)
            return phi_z(0);
        return m_phi_spline(z);
    }

    double Pi_z(const double z) const override{
        return ss_scf_common::Pi(phi_z(z), m_chi);
    }

    double D() const override{
        return m_D;
    }

    double phi_D() const override{
        return m_phi_D;
    }

};