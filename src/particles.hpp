#pragma once
#include "ss_scf_common.hpp"
#include "normalization_condition.hpp"
#include <boost/math/quadrature/trapezoidal.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>


class Particle{
    public:
        const double height, width;
        Particle(const double w, const double h) : width(w), height(h){}

        virtual double volume_integrand(const double z) const = 0;
        virtual double surface_integrand(const double z) const = 0;
        virtual const double* const surface_edges() const = 0;
        virtual double volume() const = 0;
        virtual double surface() const = 0;

        virtual double d_equivalent(){return std::cbrt(6*volume()/M_PI);}
};

namespace particle{
class Cylinder: public Particle
{
    private:
        const double m_volume, m_surface;
        const double A[2];
    public:
        Cylinder(const double h, const double w)
        : Particle(w, h), m_volume(M_PI*w*w/4*h), m_surface(M_PI*(w*h+w*w/2)), A{M_PI*w*w/4, M_PI*w*w/4}{}

        double volume_integrand(const double z) const override{
            if ((z<0)||(z>height))
            {
                return 0;
            }
            return M_PI*width*width/4;
        }

        double surface_integrand(const double z) const override{
            if ((z<0)||(z>height)){
                return 0;
            }
            return M_PI*width*2;
        }

        double volume() const override{
            return m_volume;
        }

        double surface() const override{
            return m_surface;
        }

        const double* const surface_edges() const override{
            return A;
        }
};
}