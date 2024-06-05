//
// Created by past12am on 04/06/24.
//

#ifndef QUARKDSE_V3_MOMENTUMGRID_HPP
#define QUARKDSE_V3_MOMENTUMGRID_HPP



class MomentumGrid
{
    private:
        int num_points;
        double* p2_grid;

        void generateMomentumGridLogarithmic(int num_points, double start, double center, double end);

    public:
        MomentumGrid(int num_points, double start, double center, double end);

        int getNumPoints();
        double* getMomentumGrid();

        double momentumGridAtIdx(int idx);
};



#endif //QUARKDSE_V3_MOMENTUMGRID_HPP
