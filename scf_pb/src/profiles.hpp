#pragma once
#include <algorithm>
#include "ss_scf_common.hpp"
#include <string>
#include "normalization_condition.hpp"
#include <boost/math/quadrature/trapezoidal.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>



class BrushPhiProfileBase
{
    public:
        virtual double D() const = 0;
        virtual double phi_D() const = 0;
        virtual double phi_z(double z) const =0;
};


class BrushPhiProfilePlanar : public BrushPhiProfileBase
{
    private:
    const double m_chi, m_kappa, m_R, m_N, m_sigma, m_theta;
    mutable double m_D, m_phi_D;

    public:
    BrushPhiProfilePlanar(double const &chi, double const &N, double const &sigma, double const &kappa, double const &R) : m_chi(chi), m_N(N), m_sigma(sigma), m_kappa(kappa), m_R(R), m_theta(N*sigma)
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
    BrushPhiProfilePlanar(double const &chi, double const &N, double const &sigma, double const &kappa) : BrushPhiProfilePlanar(chi, N, sigma, kappa, std::numeric_limits<double>::max()) {};
    
    double D() const override{
        return m_D;
    }

    double phi_D() const override{
        return m_phi_D;
    }

    double phi_z(const double z) const override{
        return ss_scf_common::phi_z(m_chi, m_phi_D, m_D, z, m_kappa);
    }
};

class BrushPhiProfilePore : public BrushPhiProfileBase
{
    private:
    const double m_chi, m_kappa, m_R, m_N, m_sigma, m_theta;
    mutable double m_D, m_phi_D;

    public:
    BrushPhiProfilePore(double const &chi, double const &N, double const &sigma, double const &kappa, double const &R) : m_chi(chi), m_N(N), m_sigma(sigma), m_kappa(kappa), m_R(R), m_theta(2*M_PI*R*N*sigma)
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
};

