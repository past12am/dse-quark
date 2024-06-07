//
// Created by past12am on 05/06/24.
//

#ifndef QUARKDSE_V3_COMPLEXMOMENTUMGRID_HPP
#define QUARKDSE_V3_COMPLEXMOMENTUMGRID_HPP

#include <complex>
#include <gsl/gsl_complex.h>

class ComplexMomentumGrid
{
    private:
        int num_points_real;
        int num_points_real_neg;
        int num_points_imag;

        double* p2_grid_real;
        double* p2_grid_imag;

        gsl_complex** p2_grid;

        void generateMomentumGridLogarithmic(double start, double dist_from_0, double center, double end);

    public:
        ComplexMomentumGrid(int num_points_real, int num_points_imag, int num_points_real_neg, double start, double center, double end, double dist_from_0);

        int getNumPoints();
        int getNumImagPoints();
        int getNumRealPoints();

        double* getRealGrid();
        double* getImagGrid();

        gsl_complex* getMomentumGrid(int idx_imag);

        gsl_complex momentumGridAtIdx(int idx_imag, int idx_real);
};



#endif //QUARKDSE_V3_COMPLEXMOMENTUMGRID_HPP
