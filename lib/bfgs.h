#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H

#include "states.h"
#include "matrix.h"

using namespace linopt;

template<typename functor>
point gradient(functor &f, point x, const real_type dx = 1e-8)
{
    point grad(x.size());
    real_type fx = f(x);
    real_type fxnew, tmp;
    for(int i = 0; i < x.size(); i++)
    {
        tmp = x[i];
        x[i] += dx;
        fxnew = f(x);
        x[i] = tmp;
        grad[i] = (fxnew - fx)/dx;
    }
    return grad;
}

template<typename functor>
real_type linesearch(functor &f, const real_type dt = 1e-8)
{
    const int imax = 20;
    const real_type c1 = 0.1;
    const real_type c2 = 0.9;
    real_type ft, dft;
    real_type f0 = f(0);
    real_type df0 = (f(dt) - f0)/dt;
    real_type t = 100.;
    int i = 0;
    do
    {
        t /= 2.;
        ft = f(t);
        dft = (f(t + dt) - ft)/dt;
        i++;
    }
    while(!(ft >= f0 + c1*t*df0 && dft <= c2*df0) && i < imax);
    return t;
}

class stop_criterion
{
public:
    enum status {
        CONTINUE = 0,
        ITER,
        GRAD_NORM,
        DELTA_F_ABS,
        DELTA_F_REL,
        DELTA_X_NORM};
    stop_criterion(int max_iter = 1000,
                   real_type max_grad_norm = 1e-4,
                   real_type max_delta_f_abs = 1e-10,
                   real_type max_delta_f_rel = 1e-10,
                   real_type max_delta_x_norm = 1e-10):
        max_iter(max_iter),
        max_grad_norm(max_grad_norm),
        max_delta_f_abs(max_delta_f_abs),
        max_delta_f_rel(max_delta_f_rel),
        max_delta_x_norm(max_delta_x_norm) {};
    int max_iter;
    real_type max_grad_norm;
    real_type max_delta_f_abs;
    real_type max_delta_f_rel;
    real_type max_delta_x_norm;
    status check_convergence(
            int iter,
            real_type grad_norm,
            real_type f1,
            real_type f2,
            real_type delta_x) const
    {
        if(iter >= max_iter)
            return ITER;
        if(grad_norm <= max_grad_norm)
            return GRAD_NORM;
        if(std::fabs(f2 - f1) <= max_delta_f_abs)
            return DELTA_F_ABS;
        if(std::fabs(f2/f1 - 1.) <= max_delta_f_rel)
            return DELTA_F_REL;
        if(delta_x <= max_delta_x_norm)
            return DELTA_X_NORM;
        return CONTINUE;
    }
};

template<typename functor>
real_type bfgs(functor &f, point &x, const stop_criterion &crit = stop_criterion())
{
    const int N = x.size();
    real_type t;
    point x1(N), x2(N);
    real_type f1, f2;
    point df1(N), df2(N);
    Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic> C(N, N);
    C.setIdentity();
    Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic> id = C;
    point p(N), y(N), s(N);
    real_type rho;
    int i = 0;
    x2 = x;
    f2 = f(x2);
    df2 = gradient(f, x2);
    do
    {
        x1 = x2;
        f1 = f2;
        df1 = df2;
        i++;
        p = C*df1;
        auto line_func = [&](real_type tt)
        {
            return f(x1 + tt*p);
        };
        t = linesearch(line_func);
        x2 = x1 + t*p;
        f2 = f(x2);
        df2 = gradient(f, x2);
        y = df2 - df1;
        s = x2 - x1;
        rho = 1./y.dot(s);
        C = (id - rho*s*y.transpose())*C*(id - rho*y*s.transpose()) - rho*s*s.transpose();
    }
    while(crit.check_convergence(i, df2.norm(), f1, f2, s.norm()) == stop_criterion::CONTINUE);
    x = x2;
    return f2;
}

#endif // OPTIMIZATION_H

