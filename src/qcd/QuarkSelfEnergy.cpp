//
// Created by past12am on 04/06/24.
//

#include "../../include/qcd/QuarkSelfEnergy.hpp"
#include "../../include/Definitions.hpp"

#include <cmath>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_complex_math.h>

QuarkSelfEnergy::QuarkSelfEnergy(ComplexMomentumGrid* p2Grid) : p2Grid(p2Grid)
{
    Sigma_A_imag_val_grid = new double[p2Grid->getNumPoints()];
    Sigma_A_real_val_grid = new double[p2Grid->getNumPoints()];
    Sigma_M_imag_val_grid = new double[p2Grid->getNumPoints()];
    Sigma_M_real_val_grid = new double[p2Grid->getNumPoints()];


    splineSigma_A_real = gsl_spline2d_alloc(gsl_interp2d_bilinear, p2Grid->getNumRealPoints(), p2Grid->getNumImagPoints());
    splineSigma_M_real = gsl_spline2d_alloc(gsl_interp2d_bilinear, p2Grid->getNumRealPoints(), p2Grid->getNumImagPoints());
    splineSigma_A_imag = gsl_spline2d_alloc(gsl_interp2d_bilinear, p2Grid->getNumRealPoints(), p2Grid->getNumImagPoints());
    splineSigma_M_imag = gsl_spline2d_alloc(gsl_interp2d_bilinear, p2Grid->getNumRealPoints(), p2Grid->getNumImagPoints());
}

gsl_complex QuarkSelfEnergy::getSigmaAValGrid(int imag_idx, int real_idx) const
{
    return gsl_complex_rect(gsl_spline2d_get(splineSigma_A_real, Sigma_A_real_val_grid, real_idx, imag_idx),
                            gsl_spline2d_get(splineSigma_A_imag, Sigma_A_imag_val_grid, real_idx, imag_idx));
}

gsl_complex QuarkSelfEnergy::getSigmaMValGrid(int imag_idx, int real_idx) const
{
    return gsl_complex_rect(gsl_spline2d_get(splineSigma_M_real, Sigma_M_real_val_grid, real_idx, imag_idx),
                            gsl_spline2d_get(splineSigma_M_imag, Sigma_M_imag_val_grid, real_idx, imag_idx));
}

void QuarkSelfEnergy::reclacSigma_A_Spline()
{
    // for real and imag part of Sigma, initialize a spline
    gsl_spline2d_init(splineSigma_A_imag, p2Grid->getRealGrid(), p2Grid->getImagGrid(), Sigma_A_imag_val_grid, p2Grid->getNumRealPoints(), p2Grid->getNumImagPoints());
    gsl_spline2d_init(splineSigma_A_real, p2Grid->getRealGrid(), p2Grid->getImagGrid(), Sigma_A_real_val_grid, p2Grid->getNumRealPoints(), p2Grid->getNumImagPoints());
}

void QuarkSelfEnergy::reclacSigma_M_Spline()
{
    gsl_spline2d_init(splineSigma_M_imag, p2Grid->getRealGrid(), p2Grid->getImagGrid(), Sigma_M_imag_val_grid, p2Grid->getNumRealPoints(), p2Grid->getNumImagPoints());
    gsl_spline2d_init(splineSigma_M_real, p2Grid->getRealGrid(), p2Grid->getImagGrid(), Sigma_M_real_val_grid, p2Grid->getNumRealPoints(), p2Grid->getNumImagPoints());
}

gsl_complex QuarkSelfEnergy::Sigma_M(gsl_complex p2)
{
    InterpAccels accels = {
        .interpAccel_M_X_real = gsl_interp_accel_alloc(),
        .interpAccel_M_Y_real = gsl_interp_accel_alloc(),
        .interpAccel_A_X_real = gsl_interp_accel_alloc(),
        .interpAccel_A_Y_real = gsl_interp_accel_alloc(),
        .interpAccel_M_X_imag = gsl_interp_accel_alloc(),
        .interpAccel_M_Y_imag = gsl_interp_accel_alloc(),
        .interpAccel_A_X_imag = gsl_interp_accel_alloc(),
        .interpAccel_A_Y_imag = gsl_interp_accel_alloc(),
    };

    gsl_complex Sigma_M_at_p2_val = Sigma_M(p2, &accels);

    gsl_interp_accel_free(accels.interpAccel_M_X_real);
    gsl_interp_accel_free(accels.interpAccel_M_Y_real);
    gsl_interp_accel_free(accels.interpAccel_A_X_real);
    gsl_interp_accel_free(accels.interpAccel_A_Y_real);
    gsl_interp_accel_free(accels.interpAccel_M_X_imag);
    gsl_interp_accel_free(accels.interpAccel_M_Y_imag);
    gsl_interp_accel_free(accels.interpAccel_A_X_imag);
    gsl_interp_accel_free(accels.interpAccel_A_Y_imag);

    return Sigma_M_at_p2_val;
}

gsl_complex QuarkSelfEnergy::Sigma_A(gsl_complex p2)
{
    InterpAccels accels = {
    .interpAccel_M_X_real = gsl_interp_accel_alloc(),
    .interpAccel_M_Y_real = gsl_interp_accel_alloc(),
    .interpAccel_A_X_real = gsl_interp_accel_alloc(),
    .interpAccel_A_Y_real = gsl_interp_accel_alloc(),
    .interpAccel_M_X_imag = gsl_interp_accel_alloc(),
    .interpAccel_M_Y_imag = gsl_interp_accel_alloc(),
    .interpAccel_A_X_imag = gsl_interp_accel_alloc(),
    .interpAccel_A_Y_imag = gsl_interp_accel_alloc(),
    };

    gsl_complex Sigma_A_at_p2_val = Sigma_A(p2, &accels);

    gsl_interp_accel_free(accels.interpAccel_M_X_real);
    gsl_interp_accel_free(accels.interpAccel_M_Y_real);
    gsl_interp_accel_free(accels.interpAccel_A_X_real);
    gsl_interp_accel_free(accels.interpAccel_A_Y_real);
    gsl_interp_accel_free(accels.interpAccel_M_X_imag);
    gsl_interp_accel_free(accels.interpAccel_M_Y_imag);
    gsl_interp_accel_free(accels.interpAccel_A_X_imag);
    gsl_interp_accel_free(accels.interpAccel_A_Y_imag);

    return Sigma_A_at_p2_val;
}

gsl_complex QuarkSelfEnergy::Sigma_M(gsl_complex p2, InterpAccels* interpAccels)
{
    return gsl_complex_rect(gsl_spline2d_eval(splineSigma_M_real, GSL_REAL(p2), GSL_IMAG(p2), interpAccels->interpAccel_M_X_real, interpAccels->interpAccel_M_Y_real),
                            gsl_spline2d_eval(splineSigma_M_imag, GSL_REAL(p2), GSL_IMAG(p2), interpAccels->interpAccel_M_X_imag, interpAccels->interpAccel_M_Y_imag));
}

gsl_complex QuarkSelfEnergy::Sigma_A(gsl_complex p2, InterpAccels* interpAccels)
{
    return gsl_complex_rect(gsl_spline2d_eval(splineSigma_A_real, GSL_REAL(p2), GSL_IMAG(p2), interpAccels->interpAccel_A_X_real, interpAccels->interpAccel_A_Y_real),
                            gsl_spline2d_eval(splineSigma_A_imag, GSL_REAL(p2), GSL_IMAG(p2), interpAccels->interpAccel_A_X_imag, interpAccels->interpAccel_A_Y_imag));
}

void QuarkSelfEnergy::setSigmaMAt(int imag_idx, int real_idx, gsl_complex val)
{
    gsl_spline2d_set(splineSigma_M_real, Sigma_M_real_val_grid, real_idx, imag_idx, GSL_REAL(val));
    gsl_spline2d_set(splineSigma_M_imag, Sigma_M_imag_val_grid, real_idx, imag_idx, GSL_IMAG(val));
}

void QuarkSelfEnergy::setSigmaAAt(int imag_idx, int real_idx, gsl_complex val)
{
    gsl_spline2d_set(splineSigma_A_real, Sigma_A_real_val_grid, real_idx, imag_idx, GSL_REAL(val));
    gsl_spline2d_set(splineSigma_A_imag, Sigma_A_imag_val_grid, real_idx, imag_idx, GSL_IMAG(val));
}
