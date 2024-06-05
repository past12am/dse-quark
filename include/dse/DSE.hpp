//
// Created by past12am on 04/06/24.
//

#ifndef QUARKDSE_V3_DSE_HPP
#define QUARKDSE_V3_DSE_HPP



#include "../numerics/integration/GaussChebyshev.hpp"
#include "../numerics/integration/GaussLegendre.hpp"
#include "../qcd/QuarkSelfEnergy.hpp"
#include "../qcd/QuarkPropagator.hpp"

class DSE
{
    private:
        double Z_2 = 1.0;
        double L2;

        MomentumGrid* p2Grid;
        QuarkPropagator* quarkPropagator;
        QuarkSelfEnergy* selfEnergy;

        GaussLegendre* gaussLegendreIntegrator;
        GaussChebyshev* gaussChebyshevIntegrator;

        double F(double p2, double q2, double k2, double z);
        double g(double k2);
        double alpha(double k2);
        double calc_k2(double p2, double q2, double z);

        double selfEnergyIntegralKernelSigma_A(double p2, double q2, double z);
        double selfEnergyIntegralKernelSigma_M(double p2, double q2, double z);

        double performIntegration_Sigma_M(double p2);
        double performIntegration_Sigma_A(double p2);

        double zIntegral(double p2, double q2, const std::function<double(double, double, double)>& f);
        double q2Integral(double p2, const std::function<double(double, double, double)> &f, double upper_cutoff);

    public:
        DSE(double L2, MomentumGrid* p2Grid);

        void solveDSE();

        QuarkPropagator *getQuarkPropagator() const;
};



#endif //QUARKDSE_V3_DSE_HPP
