 #include "brush_profiles.hpp"

BrushProfilePlanar::BrushProfilePlanar(double const &chi, double const &N, double const &sigma, double const &kappa, double const &R)
    : m_chi(chi), m_N(N), m_sigma(sigma), m_kappa(kappa), m_R(R), m_theta(N * sigma)
{
        double D_free = solve_normalization(planar::free::D(m_chi, m_theta, m_kappa), m_theta, m_N);
        // free brush
        if (D_free < m_R)
        {
                m_D = D_free;
                m_phi_D = phi_at_zero_Pi(m_chi);
        }
        else
        {
                m_phi_D = solve_normalization(planar::restricted::phi_D(m_chi, m_theta, m_kappa, m_R), 0.0, double(ALMOST_ONE));
                m_D = m_R;
        }
}

double BrushProfilePlanar::phi_z(const double z) const
{
        return ss_scf_common::phi_z(m_chi, m_phi_D, m_D, z, m_kappa);
}

double BrushProfilePlanar::Pi_z(const double z) const
{
        return ss_scf_common::Pi(phi_z(z), m_chi);
}

BrushProfilePore::BrushProfilePore(double const &chi, double const &N, double const &sigma, double const &kappa, double const &R)
    : m_chi(chi), m_N(N), m_sigma(sigma), m_kappa(kappa), m_R(R), m_theta(2 * M_PI * R * N * sigma)
{
        double D_free = solve_normalization(pore::free::D(m_chi, m_theta, m_kappa, m_R), m_N * m_sigma, m_N);
        // free brush
        if (D_free < m_R)
        {
                m_D = D_free;
                m_phi_D = phi_at_zero_Pi(m_chi);
        }
        else
        {
                m_phi_D = solve_normalization(pore::restricted::phi_D(m_chi, m_theta, m_kappa, m_R), 0.0, double(ALMOST_ONE));
                m_D = m_R;
        }
}

double BrushProfilePore::phi_z(const double z) const
{
        return ss_scf_common::phi_z(m_chi, m_phi_D, m_D, z, m_kappa);
}

double BrushProfilePore::Pi_z(const double z) const
{
        return ss_scf_common::Pi(phi_z(z), m_chi);
}

BrushProfileExternal::BrushProfileExternal(const std::vector<double> &phi, const double chi, const double dz, const double phi_cut)
    : m_phi_profile(phi), m_dz(dz), m_phi_cut(phi_cut), m_chi(chi), m_phi_spline(phi.begin(), phi.end(), 0, dz)
{
        auto first_smaller_it = std::lower_bound(m_phi_profile.rbegin(), m_phi_profile.rend(), m_phi_cut);
        double idx = std::distance(first_smaller_it, m_phi_profile.rend()) * m_dz;
        m_D = idx * dz;
        m_phi_D = m_phi_profile[idx];
        // m_D = boost::math::tools::bisect(m_phi_spline, 0, m_phi_profile.size()*dx, phi_cut);
}

double BrushProfileExternal::phi_z(const double z) const
{
        if (z > m_D)
                return 0.0;
        if (z < 0)
                return phi_z(0);
        return m_phi_spline(z);
}

double BrushProfileExternal::Pi_z(const double z) const
{
        return ss_scf_common::Pi(phi_z(z), m_chi);
}