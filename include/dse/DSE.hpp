//
// Created by past12am on 04/06/24.
//

#ifndef QUARKDSE_V3_DSE_HPP
#define QUARKDSE_V3_DSE_HPP



#include <gsl/gsl_complex_math.h>
#include "../numerics/integration/GaussChebyshev.hpp"
#include "../numerics/integration/GaussLegendre.hpp"
#include "../qcd/QuarkSelfEnergy.hpp"
#include "../qcd/QuarkPropagator.hpp"
#include "../grid/ComplexMomentumGrid.hpp"

class DSE
{
    private:
        gsl_complex Z_2 = GSL_COMPLEX_ONE;
        double L2;

        ComplexMomentumGrid* p2Grid;
        QuarkPropagator* quarkPropagator;
        QuarkSelfEnergy* selfEnergy;

        GaussLegendre* gaussLegendreIntegrator;
        GaussChebyshev* gaussChebyshevIntegrator;

        gsl_complex F(gsl_complex p2, gsl_complex q2, gsl_complex k2, double z);
        gsl_complex g(gsl_complex k2);
        gsl_complex alpha(gsl_complex k2);
        gsl_complex calc_k2(gsl_complex p2, gsl_complex q2, double z);

        gsl_complex selfEnergyIntegralKernelSigma_A(gsl_complex p2, gsl_complex q2, double z, InterpAccels* interpAccels);
        gsl_complex selfEnergyIntegralKernelSigma_M(gsl_complex p2, gsl_complex q2, double z, InterpAccels* interpAccels);

        gsl_complex performIntegration_Sigma_M(gsl_complex p2, InterpAccels* interpAccels);
        gsl_complex performIntegration_Sigma_A(gsl_complex p2, InterpAccels* interpAccels);

        gsl_complex zIntegral(gsl_complex p2, gsl_complex q2, const std::function<gsl_complex(gsl_complex , gsl_complex , double)>& f);
        gsl_complex q2Integral(gsl_complex p2, const std::function<gsl_complex(gsl_complex, gsl_complex, double)> &f, double upper_cutoff);

        InterpAccels allocInterpAccels();
        void freeInterpAccels(InterpAccels interpAccels);

    public:
        DSE(double L2, ComplexMomentumGrid* p2Grid);

        void solveDSE();

        QuarkPropagator *getQuarkPropagator() const;
};



#endif //QUARKDSE_V3_DSE_HPP
