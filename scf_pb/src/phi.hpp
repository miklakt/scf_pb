#pragma once
#include "ss_scf_common.hpp"
#include <string>
#include "normalization_condition.hpp"
#include <boost/math/quadrature/trapezoidal.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>

enum BRUSH_GEOMETRY: std::int8_t{
    PLANAR=0,
    PORE=1
};

template <typename T>
class BrushEdge{
    public:
        virtual T D() const;
        virtual T phi_D() const;
};

template <typename T>
class BrushEdgePlanar : public BrushEdge<T>{
    BrushEdgePlanar(T const &chi, T const &theta, T const &kappa, T const &R = std::numeric_limits<T>::max()){
        m_chi=chi;
        m_theta=theta;
        m_kappa = kappa;
        m_R = R;
        solve();
    };
    void solve() const{
        T D_free = planar::open(m_chi, m_theta, m_kappa);
        if (D_free<m_R){
            m_D = D_free;
            m_phi_D = phi_at_zero_Pi(m_chi);
        }
        else{
            m_phi_D = solve_normalization(planar::restricted(m_chi, m_theta, m_kappa, m_R));
            m_D = m_R;
        }
    };

    T phi_D(){
        return m_phi_D;
    }

    T D(){
        return m_D;
    }

    private:
    T m_chi, m_theta, m_kappa, m_R, m_D, m_phi_D;
};