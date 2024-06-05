//
// Created by past12am on 04/06/24.
//

#ifndef QUARKDSE_V3_QUARKPROPAGATOR_HPP
#define QUARKDSE_V3_QUARKPROPAGATOR_HPP



#include "QuarkSelfEnergy.hpp"



class QuarkPropagator
{
    private:
        QuarkSelfEnergy* selfEnergy;

    public:
        QuarkPropagator(QuarkSelfEnergy* selfEnergy);

        double M(double p2);
        double A(double p2);

        double M(double p2, gsl_interp_accel* interpAccel_M, gsl_interp_accel* interpAccel_A);
        double A(double p2, gsl_interp_accel* interpAccel_A);

        double sigma_v(double p2);
        double sigma_s(double p2);

        double sigma_v(double p2, gsl_interp_accel* interpAccel_M, gsl_interp_accel* interpAccel_A);
        double sigma_s(double p2, gsl_interp_accel* interpAccel_M, gsl_interp_accel* interpAccel_A);
};



#endif //QUARKDSE_V3_QUARKPROPAGATOR_HPP
