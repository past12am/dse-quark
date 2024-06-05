//
// Created by past12am on 04/06/24.
//

#include "../../include/grid/MomentumGrid.hpp"

#include <cmath>
#include <iostream>
#include <cassert>

MomentumGrid::MomentumGrid(int num_points, double start, double center, double end) : num_points(num_points)
{
    p2_grid = new double[num_points];
    generateMomentumGridLogarithmic(num_points, start, center, end);
}

void MomentumGrid::generateMomentumGridLogarithmic(int num_points, double start, double center, double end)
{
    // Logarithmic Distribution of grid points
    double p2_grid_exp_start = log10(start);
    double p2_grid_exp_center = log10(center);
    double p2_grid_exp_end = log10(end);

    double dist_upper = (p2_grid_exp_end - p2_grid_exp_center) / (num_points/2);
    double dist_lower = (p2_grid_exp_start - p2_grid_exp_center) / (num_points/2);

    for (int i = 0; i <= (num_points)/2; i++)
    {
        p2_grid[(num_points)/2 - i] = pow(10.0, p2_grid_exp_center + i * dist_lower);
        p2_grid[(num_points)/2 + i - 1] = pow(10.0, p2_grid_exp_center + i * dist_upper);
    }

    std::cout << "p2_grid: ";
    for (int i = 0; i < num_points; i++)
    {
        if(i < num_points - 1) assert(p2_grid[i] < p2_grid[i + 1]);

        std::cout << p2_grid[i] << " - ";
    }
    std::cout << std::endl << std::endl;
}

double MomentumGrid::momentumGridAtIdx(int idx)
{
    if(idx < num_points)
        return p2_grid[idx];
    else
        throw std::invalid_argument("Idx out of range for momentum grid");
}

int MomentumGrid::getNumPoints()
{
    return num_points;
}

double* MomentumGrid::getMomentumGrid()
{
    return p2_grid;
}


