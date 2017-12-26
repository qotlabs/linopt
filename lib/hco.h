#ifndef HCO_H
#define HCO_H

#include "matrix.h"
#include "bfgs.h"

using namespace linopt;

class hypercube
{
private:
    point c;
    point r;
public:
    explicit hypercube(const point &center = point(),
                       const point &radius = point()):
        c(center), r(radius) {}

    point lower_bound() const
    {
        return c - r;
    }

    hypercube &set_lower_bound(const point &lb)
    {
        point ub = c + r;
        c = (ub + lb)/2.;
        r = (ub - lb)/2.;
        return *this;
    }

    point upper_bound() const
    {
        return c + r;
    }

    hypercube &set_upper_bound(const point &ub)
    {
        point lb = c - r;
        c = (ub + lb)/2.;
        r = (ub - lb)/2.;
        return *this;
    }

    point center() const
    {
        return c;
    }

    hypercube &set_center(const point &center)
    {
        c = center;
        return *this;
    }

    point radius() const
    {
        return r;
    }

    hypercube &set_radius(const point &radius)
    {
        r = radius;
        return *this;
    }
    point random_point()
    {
        point x = point::Random(c.size());
        return x.cwiseProduct(r) + c;
    }
};

#include <iostream>

template<typename functor>
real_type hco(functor &f, point &x, const stop_criterion &crit = stop_criterion())
{
    const int dim = x.size();
    const int N = 1;
    stop_criterion bfgs_stop(1000, 0, 0, 1e-3, 0);
    point xbest = point::Ones(dim)/2.;
    real_type fbest = f(xbest);
    real_type fx;
    real_type s;
    hypercube hc(xbest, xbest);
    int i = 0;
    do
    {
        for(int n = 0; n < N; n++)
        {
            x = hc.random_point();
            fx = f(x);
            std::cout << fx << std::endl;
            if(fx > fbest)
            {
                xbest = x;
                fbest = fx;
            }
        }
        fbest = bfgs(f, xbest, bfgs_stop);
        std::cout << "fbest = " << fbest << "\tr = " << hc.radius().norm() << std::endl;
        s = ((xbest - hc.center()).cwiseQuotient(hc.radius())).norm();
        s /= 4.*std::sqrt(dim);
        s = 1. - 0.2*std::exp(-3.*s);
        hc.set_center((xbest + hc.center())/2.);
        hc.set_radius(hc.radius()*s);
        i++;
    }
    while(i < crit.max_iter);
    x = xbest;
    return fbest;
}

#endif // HCO_H

