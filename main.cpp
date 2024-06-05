#include <iostream>
#include <sstream>
#include <fstream>
#include "include/dse/DSE.hpp"
#include "include/Definitions.hpp"

int main()
{
    MomentumGrid* p2Grid = new MomentumGrid(1000, 1E-3, 10, 1E3);
    DSE* dse = new DSE(1E3, p2Grid);

    clock_t start = clock();

    dse->solveDSE();

    clock_t stop = clock();
    double elapsed = (double) (stop - start) / CLOCKS_PER_SEC;
    printf("\nTime elapsed: %.5f\n", elapsed);



    // Save M(p2) and A(p2) to files
    QuarkPropagator* quarkPropagator = dse->getQuarkPropagator();

    std::ostringstream fnamestrstream_A_p2;
    fnamestrstream_A_p2 << "/home/past12am/OuzoCloud/Studium/Physik/6_Semester/SE_Bachelorarbeit/QCD_Intro/QuarkDSE/output/A_p2_m_" << QUARK_MASS << ".txt";
    std::string filename_A_p2 = fnamestrstream_A_p2.str();

    std::ostringstream fnamestrstream_M_p2;
    fnamestrstream_M_p2 << "/home/past12am/OuzoCloud/Studium/Physik/6_Semester/SE_Bachelorarbeit/QCD_Intro/QuarkDSE/output/M_p2_m_" << QUARK_MASS << ".txt";
    std::string filename_M_p2 = fnamestrstream_M_p2.str();

    std::ofstream A_p2_file;
    A_p2_file.open(filename_A_p2, std::ofstream::out | std::ios::trunc);

    std::ofstream M_p2_file;
    M_p2_file.open(filename_M_p2, std::ofstream::out | std::ios::trunc);

    A_p2_file << "p2;A" << std::endl;
    M_p2_file << "p2;M" << std::endl;

    for (int i = 0; i < p2Grid->getNumPoints(); i++)
    {
        A_p2_file << p2Grid->momentumGridAtIdx(i) << ";" << quarkPropagator->A(p2Grid->momentumGridAtIdx(i)) << std::endl;
        M_p2_file << p2Grid->momentumGridAtIdx(i) << ";" << quarkPropagator->M(p2Grid->momentumGridAtIdx(i)) << std::endl;
    }

    A_p2_file.close();
    M_p2_file.close();

    return 0;
}
