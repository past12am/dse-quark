//
// Created by past12am on 04/06/24.
//

#include "../../include/qcd/QuarkSelfEnergy.hpp"
#include "../../include/Definitions.hpp"

#include <cmath>
#include <gsl/gsl_spline.h>

QuarkSelfEnergy::QuarkSelfEnergy(MomentumGrid* p2Grid) : p2Grid(p2Grid)
{
    Sigma_A_val_grid = new double[p2Grid->getNumPoints()];
    Sigma_M_val_grid = new double[p2Grid->getNumPoints()];

    splineSigma_A = gsl_spline_alloc(gsl_interp_cspline, p2Grid->getNumPoints());
    splineSigma_M = gsl_spline_alloc(gsl_interp_cspline, p2Grid->getNumPoints());
}

double *QuarkSelfEnergy::getSigmaAValGrid() const
{
    return Sigma_A_val_grid;
}

double *QuarkSelfEnergy::getSigmaMValGrid() const
{
    return Sigma_M_val_grid;
}

void QuarkSelfEnergy::reclacSigma_A_Spline()
{
    gsl_spline_init(splineSigma_A, p2Grid->getMomentumGrid(), Sigma_A_val_grid, p2Grid->getNumPoints());
}

void QuarkSelfEnergy::reclacSigma_M_Spline()
{
    gsl_spline_init(splineSigma_M, p2Grid->getMomentumGrid(), Sigma_M_val_grid, p2Grid->getNumPoints());
}

double QuarkSelfEnergy::Sigma_M(double p2)
{
    interpAccelSigma_M = gsl_interp_accel_alloc();
    double Sigma_M_at_p2_val = Sigma_M(p2, interpAccelSigma_M);
    gsl_interp_accel_free(interpAccelSigma_M);

    return Sigma_M_at_p2_val;
}

double QuarkSelfEnergy::Sigma_A(double p2)
{
    interpAccelSigma_A = gsl_interp_accel_alloc();
    double Sigma_A_at_p2_val = Sigma_A(p2, interpAccelSigma_A);
    gsl_interp_accel_free(interpAccelSigma_A);

    return Sigma_A_at_p2_val;
}

double QuarkSelfEnergy::Sigma_M(double p2, gsl_interp_accel* splineInterpAccel)
{
    return gsl_spline_eval(splineSigma_M, p2, splineInterpAccel);
}

double QuarkSelfEnergy::Sigma_A(double p2, gsl_interp_accel* splineInterpAccel)
{
    return gsl_spline_eval(splineSigma_A, p2, splineInterpAccel);
}

void QuarkSelfEnergy::allocateSplineAccel_M()
{
    interpAccelSigma_M = gsl_interp_accel_alloc();
}

void QuarkSelfEnergy::allocateSplineAccel_A()
{
    interpAccelSigma_A = gsl_interp_accel_alloc();
}

void QuarkSelfEnergy::freeSplineAccel_M()
{
    gsl_interp_accel_free(interpAccelSigma_M);
}

void QuarkSelfEnergy::freeSplineAccel_A()
{
    gsl_interp_accel_free(interpAccelSigma_A);
}

void QuarkSelfEnergy::setSigmaMAt(int idx, double val)
{
    Sigma_M_val_grid[idx] = val;
}

void QuarkSelfEnergy::setSigmaAAt(int idx, double val)
{
    Sigma_A_val_grid[idx] = val;
}
