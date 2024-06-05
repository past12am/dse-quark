//
// Created by past12am on 04/06/24.
//

#include "../../include/qcd/QuarkPropagator.hpp"
#include "../../include/Definitions.hpp"

#include <cmath>

QuarkPropagator::QuarkPropagator(QuarkSelfEnergy* selfEnergy) : selfEnergy(selfEnergy)
{

}

double QuarkPropagator::sigma_v(double p2)
{
    return 1.0 / (A(p2) * (p2 + pow(M(p2), 2)));
}

double QuarkPropagator::sigma_v(double p2, gsl_interp_accel* interpAccel_M, gsl_interp_accel* interpAccel_A)
{
    return 1.0 / (A(p2, interpAccel_A) * (p2 + pow(M(p2, interpAccel_M, interpAccel_A), 2)));
}

double QuarkPropagator::sigma_s(double p2)
{
    return sigma_v(p2) * M(p2);
}

double QuarkPropagator::sigma_s(double p2, gsl_interp_accel* interpAccel_M, gsl_interp_accel* interpAccel_A)
{
    return sigma_v(p2, interpAccel_M, interpAccel_A) * M(p2, interpAccel_M, interpAccel_A);
}

double QuarkPropagator::M(double p2)
{
    return 1.0/A(p2) * (QUARK_MASS + selfEnergy->Sigma_M(p2) - selfEnergy->Sigma_M(RENORM_POINT2));
}

double QuarkPropagator::M(double p2, gsl_interp_accel* interpAccel_M, gsl_interp_accel* interpAccel_A)
{
    return 1.0/A(p2, interpAccel_A) * (QUARK_MASS + selfEnergy->Sigma_M(p2, interpAccel_M) - selfEnergy->Sigma_M(RENORM_POINT2, interpAccel_M));
}

double QuarkPropagator::A(double p2)
{
    return 1.0 + selfEnergy->Sigma_A(p2) - selfEnergy->Sigma_A(RENORM_POINT2);
}

double QuarkPropagator::A(double p2, gsl_interp_accel* interpAccel_A)
{
    return 1.0 + selfEnergy->Sigma_A(p2, interpAccel_A) - selfEnergy->Sigma_A(RENORM_POINT2, interpAccel_A);
}