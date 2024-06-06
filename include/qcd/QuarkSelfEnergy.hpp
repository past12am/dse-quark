//
// Created by past12am on 04/06/24.
//

#ifndef QUARKDSE_V3_QUARKSELFENERGY_HPP
#define QUARKDSE_V3_QUARKSELFENERGY_HPP


#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

#include "../helpers/InterpAcceleratorsStruct.hpp"
#include "../grid/ComplexMomentumGrid.hpp"

class QuarkSelfEnergy
{
    private:
        ComplexMomentumGrid* p2Grid;

        // Splines
        double* Sigma_A_imag_val_grid;
        double* Sigma_A_real_val_grid;
        double* Sigma_M_imag_val_grid;
        double* Sigma_M_real_val_grid;


        gsl_spline2d* splineSigma_A_imag;
        gsl_spline2d* splineSigma_M_imag;
        gsl_spline2d* splineSigma_A_real;
        gsl_spline2d* splineSigma_M_real;

        gsl_interp_accel* interpAccelX_real;
        gsl_interp_accel* interpAccelY_real;
        gsl_interp_accel* interpAccelX_imag;
        gsl_interp_accel* interpAccelY_imag;

    public:
        QuarkSelfEnergy(ComplexMomentumGrid* p2Grid);

        gsl_complex Sigma_M(gsl_complex p2);
        gsl_complex Sigma_A(gsl_complex p2);

        gsl_complex Sigma_M(gsl_complex p2, InterpAccels* interpAccels);
        gsl_complex Sigma_A(gsl_complex p2, InterpAccels* interpAccels);

        void reclacSigma_A_Spline();
        void reclacSigma_M_Spline();

        gsl_complex getSigmaAValGrid(int imag_idx, int real_idx) const;
        gsl_complex getSigmaMValGrid(int imag_idx, int real_idx) const;

        void setSigmaAAt(int imag_idx, int real_idx, gsl_complex val);
        void setSigmaMAt(int imag_idx, int real_idx, gsl_complex val);
};



#endif //QUARKDSE_V3_QUARKSELFENERGY_HPP
