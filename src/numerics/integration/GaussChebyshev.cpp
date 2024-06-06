//
// Created by past12am on 3/2/23.
//

#include "../../../include/numerics/integration/GaussChebyshev.hpp"

#include <cmath>
#include <numbers>

#include <functional>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

gsl_complex GaussChebyshev::integrate_f_times_sqrt(std::function<gsl_complex(double)>& f)
{
    gsl_complex val = GSL_COMPLEX_ZERO;
    double a = M_PI/(n + 1.0);
    for (int i = 0; i < n; i++)
    {
        val = gsl_complex_add(val, gsl_complex_mul_real(f(a * i), (a * pow(sin(a * i), 2))));
    }
    return val;
}

GaussChebyshev::GaussChebyshev(int n) : n(n)
{

}
