//
// Created by past12am on 06/06/24.
//

#ifndef QUARKDSE_V3_INTERPACCELERATORSSTRUCT_HPP
#define QUARKDSE_V3_INTERPACCELERATORSSTRUCT_HPP

#include <gsl/gsl_spline2d.h>

struct InterpAccels {
    gsl_interp_accel* interpAccel_M_X_real;
    gsl_interp_accel* interpAccel_M_Y_real;
    gsl_interp_accel* interpAccel_A_X_real;
    gsl_interp_accel* interpAccel_A_Y_real;
    gsl_interp_accel* interpAccel_M_X_imag;
    gsl_interp_accel* interpAccel_M_Y_imag;
    gsl_interp_accel* interpAccel_A_X_imag;
    gsl_interp_accel* interpAccel_A_Y_imag;
};

#endif //QUARKDSE_V3_INTERPACCELERATORSSTRUCT_HPP
