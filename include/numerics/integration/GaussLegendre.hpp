//
// Created by past12am on 3/2/23.
//

#ifndef QUARKDSE_GAUSSLEGENDRE_HPP
#define QUARKDSE_GAUSSLEGENDRE_HPP

#include <cstddef>
#include <tuple>

#include <functional>

#include <gsl/gsl_complex.h>

class GaussLegendre
{
    private:
        int n;

        double* w_arr;
        double* x_arr;

        static std::tuple<double*, double*> generageWeights(int n);

    public:
        GaussLegendre(int n);
        virtual ~GaussLegendre();

        gsl_complex integrate(std::function<gsl_complex(gsl_complex)>& f, double a, double b);
};

#endif //QUARKDSE_GAUSSLEGENDRE_HPP
