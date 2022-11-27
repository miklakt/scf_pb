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

        virtual double z0(const double particle_center) const {return particle_center - height/2;}
        virtual double z1(const double particle_center) const {return particle_center + height/2;}

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


class Sphere: public Particle{
    private:
    const double m_volume, m_surface;
    const double A[2];
    const double radius;
    double chord_r2(const double z) const{
        //double z_origin = 
        return (radius*radius-(z-radius)*(z-radius));
    }
    public:
    Sphere(const double r): Particle(2*r, 2*r), radius(r), m_volume(M_PI*r*r*r*4/3), m_surface(4*M_PI*r*r), A{0, 0}{}
    
    double volume_integrand(const double z) const override{
        if ((z<0)||(z>height))
        {
            return 0;
        }
        return M_PI*chord_r2(z);
    }

    double surface_integrand(const double z) const override{
        if ((z<0)||(z>height))
        {
            return 0;
        }
        //return 2*M_PI*std::sqrt(chord_r2(z));
        return 2*M_PI*radius;
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

}//namespace particle


namespace particle::integrators{
    template <typename ParticleType, typename FunctionToIntegrate>
    double integrate_over_volume(const ParticleType *particle, const FunctionToIntegrate func, const double a, const double b)
    {   
        if (a==b){return 0;}
        auto integrand = [=](const double z){
            return particle->volume_integrand(z)*func(z);
        };
        return boost::math::quadrature::gauss_kronrod<double, 31>::integrate(integrand, a, b);
    }
    
    template <typename ParticleType, typename FunctionToIntegrate>
    double integrate_over_volume(const ParticleType *particle, const FunctionToIntegrate func, const double b)
    {
        return integrate_over_volume(particle, func, 0.0, b);
    }

    template <typename ParticleType, typename FunctionToIntegrate>
    double integrate_over_volume(const ParticleType *particle, const FunctionToIntegrate func)
    {
        return integrate_over_volume(particle, func, 0.0, particle->height);
    }

    template <typename ParticleType, typename FunctionToIntegrate>
    double integrate_over_surface(const ParticleType *particle, const FunctionToIntegrate func, const double a, const double b)
    {   
        if (a==b){return 0;}
        auto integrand = [=](const double z){
            return particle->surface_integrand(z)*func(z);
        };
        auto A = particle->surface_edges();
        double caps = (a==0.0)*A[0]*func(a)+(b==particle->height)*A[1]*func(b);
        return boost::math::quadrature::gauss_kronrod<double, 31>::integrate(integrand, a, b) + caps;
    }

    template <typename ParticleType, typename FunctionToIntegrate>
    double integrate_over_surface(const ParticleType *particle, const FunctionToIntegrate func, const double b){
        return integrate_over_surface(particle, func, 0.0, b);
    }

    template <typename ParticleType, typename FunctionToIntegrate>
    double integrate_over_surface(const ParticleType *particle, const FunctionToIntegrate func){
        return integrate_over_surface(particle, func, 0.0, particle->height);
    }

    
    //simple integration over cylinder
    template <typename FunctionToIntegrate>
    double integrate_over_volume(const particle::Cylinder *particle, const FunctionToIntegrate func, const double a, const double b)
    {   
        if (a==b){return 0;}
        auto integrand = [=](const double z){
            return func(z);
        };
        return boost::math::quadrature::gauss_kronrod<double, 31>::integrate(integrand, a, b)*particle->surface_edges()[0];
    }

    template <typename FunctionToIntegrate>
    double integrate_over_surface(const particle::Cylinder *particle, const FunctionToIntegrate func, const double a, const double b)
    {   
        if (a==b){return 0;}
        auto integrand = [=](const double z){
            return func(z);
        };
        auto A = particle->surface_edges();
        double caps = (a==0.0)*A[0]*func(a)+(b==particle->height)*A[1]*func(b);
        return boost::math::quadrature::gauss_kronrod<double, 31>::integrate(integrand, a, b)*M_PI*particle->width*2 + caps;
    }
    
}//namespace integrators