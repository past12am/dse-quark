//
// Created by past12am on 04/06/24.
//

#ifndef QUARKDSE_V3_QUARKSELFENERGY_HPP
#define QUARKDSE_V3_QUARKSELFENERGY_HPP


#include <gsl/gsl_spline.h>
#include "../grid/MomentumGrid.hpp"

class QuarkSelfEnergy
{
    private:
        MomentumGrid* p2Grid;

        // Splines
        double* Sigma_A_val_grid;
        double* Sigma_M_val_grid;

        gsl_spline* splineSigma_A;
        gsl_spline* splineSigma_M;

        gsl_interp_accel* interpAccelSigma_A;
        gsl_interp_accel* interpAccelSigma_M;

    public:
        QuarkSelfEnergy(MomentumGrid* p2Grid);

        double Sigma_M(double p2);
        double Sigma_A(double p2);

        double Sigma_M(double p2, gsl_interp_accel* splineInterpAccel);
        double Sigma_A(double p2, gsl_interp_accel* splineInterpAccel);

        void allocateSplineAccel_M();
        void allocateSplineAccel_A();

        void freeSplineAccel_M();
        void freeSplineAccel_A();

        void reclacSigma_A_Spline();
        void reclacSigma_M_Spline();

        double *getSigmaAValGrid() const;
        double *getSigmaMValGrid() const;

        void setSigmaAAt(int idx, double val);
        void setSigmaMAt(int idx, double val);
};



#endif //QUARKDSE_V3_QUARKSELFENERGY_HPP
