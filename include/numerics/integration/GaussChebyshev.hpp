//
// Created by past12am on 3/2/23.
//

#ifndef QUARKDSE_GAUSSCHEBYSHEV_HPP
#define QUARKDSE_GAUSSCHEBYSHEV_HPP

#include <stddef.h>
#include <tuple>

#include <functional>

#include <gsl/gsl_complex.h>

class GaussChebyshev
{
    private:
        int n;

    public:
        GaussChebyshev(int n);

        gsl_complex integrate_f_times_sqrt(std::function<gsl_complex(double)>& f);
};


#endif //QUARKDSE_GAUSSCHEBYSHEV_HPP
