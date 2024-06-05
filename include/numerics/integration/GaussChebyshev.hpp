//
// Created by past12am on 3/2/23.
//

#ifndef QUARKDSE_GAUSSCHEBYSHEV_HPP
#define QUARKDSE_GAUSSCHEBYSHEV_HPP

#include <stddef.h>
#include <tuple>

#include <functional>

class GaussChebyshev
{
    private:
        int n;

    public:
        GaussChebyshev(int n);

        double integrate_f_times_sqrt(std::function<double(double)>& f);
};


#endif //QUARKDSE_GAUSSCHEBYSHEV_HPP
