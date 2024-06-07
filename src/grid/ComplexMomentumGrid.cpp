//
// Created by past12am on 05/06/24.
//

#include "../../include/grid/ComplexMomentumGrid.hpp"

#include <iostream>
#include <cmath>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

ComplexMomentumGrid::ComplexMomentumGrid(int num_points_real, int num_points_imag, int num_points_real_neg, double start, double center, double end, double dist_from_0) : num_points_real(num_points_real), num_points_imag(num_points_imag), num_points_real_neg(num_points_real_neg)
{
    p2_grid = new gsl_complex*[num_points_imag];
    for (int i = 0; i < num_points_imag; i++)
    {
        p2_grid[i] = new gsl_complex[num_points_real + num_points_real_neg];
    }

    p2_grid_real = new double[num_points_real + num_points_real_neg];
    p2_grid_imag = new double[num_points_imag];

    generateMomentumGridLogarithmic(start, dist_from_0, center, end);
}

void ComplexMomentumGrid::generateMomentumGridLogarithmic(double start, double dist_from_0, double center, double end)
{
    double complex_dist = num_points_imag > 1 ? 2.0/((double) (num_points_imag - 1)) : 1;

    if(start < 0)
    {
        // Create a linear negative real axis part
        double neg_dist = (abs(start) - dist_from_0) / (num_points_real_neg - 1);

        for (int j = 0; j < num_points_imag; j++)
        {
            double imag_part = num_points_imag > 1 ? -1 + j * complex_dist : 0;

            for (int i = 0; i < num_points_real_neg; i++)
            {
                double real_part = start + i * neg_dist;

                p2_grid[j][i] = gsl_complex_rect(real_part, imag_part);
            }
        }
    }

    // Logarithmic Distribution of grid points
    double p2_grid_exp_start = log10(dist_from_0);
    double p2_grid_exp_center = log10(center);
    double p2_grid_exp_end = log10(end);

    double dist_upper = (p2_grid_exp_end - p2_grid_exp_center) / (num_points_real/2);
    double dist_lower = (p2_grid_exp_start - p2_grid_exp_center) / (num_points_real/2);

    for (int j = 0; j < num_points_imag; j++)
    {
        for (int i = 0; i <= (num_points_real)/2; i++)
        {
            double imag_part = num_points_imag > 1 ? -1 + j * complex_dist : 0;

            p2_grid[j][num_points_real_neg + (num_points_real)/2 - i] = gsl_complex_rect(pow(10.0, p2_grid_exp_center + i * dist_lower), imag_part);
            p2_grid[j][num_points_real_neg + (num_points_real)/2 + i - 1] = gsl_complex_rect(pow(10.0, p2_grid_exp_center + i * dist_upper), imag_part);
        }

        std::cout << "p2_grid: ";
        for (int i = 0; i < num_points_real + num_points_real_neg; i++)
        {
            std::cout << GSL_REAL(p2_grid[j][i]) << ((GSL_IMAG(p2_grid[j][i]) < 0) ? " - " : " + ") << std::abs(GSL_IMAG(p2_grid[j][i]))  << "i   ;  ";
        }
        std::cout << std::endl << std::endl;
    }

    for (int i = 0; i < num_points_real + num_points_real_neg; i++)
    {
        p2_grid_real[i] = GSL_REAL(p2_grid[0][i]);
    }
    for (int j = 0; j < num_points_imag; j++)
    {
        p2_grid_imag[j] = GSL_IMAG(p2_grid[j][0]);
    }
}

gsl_complex ComplexMomentumGrid::momentumGridAtIdx(int idx_imag, int idx_real)
{
    if(idx_real < getNumRealPoints() && idx_imag < getNumImagPoints())
        return p2_grid[idx_imag][idx_real];
    else
        throw std::invalid_argument("Idx out of range for momentum grid idx_imag = " + std::to_string(idx_imag) + "; idx_real = " + std::to_string(idx_real));
}

int ComplexMomentumGrid::getNumPoints()
{
    return (num_points_real + num_points_real_neg) * num_points_imag;
}

gsl_complex* ComplexMomentumGrid::getMomentumGrid(int idx_imag)
{
    return p2_grid[idx_imag];
}

int ComplexMomentumGrid::getNumImagPoints()
{
    return num_points_imag;
}

int ComplexMomentumGrid::getNumRealPoints()
{
    return num_points_real + num_points_real_neg;
}

double *ComplexMomentumGrid::getRealGrid()
{
    return p2_grid_real;
}

double *ComplexMomentumGrid::getImagGrid()
{
    return p2_grid_imag;
}
