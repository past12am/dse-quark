//
// Created by past12am on 3/2/23.
//

#include "../../../include/numerics/integration/GaussChebyshev.hpp"

#include <cmath>
#include <numbers>

#include <functional>

double GaussChebyshev::integrate_f_times_sqrt(std::function<double(double)>& f)
{
    double val = 0;
    double a = M_PI/(n + 1.0);
    for (int i = 0; i < n; i++)
    {
        val += (a * pow(sin(a * i), 2)) * f(cos(a * i));
    }
    return val;
}

GaussChebyshev::GaussChebyshev(int n) : n(n)
{

}
