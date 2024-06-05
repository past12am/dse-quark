//
// Created by past12am on 04/06/24.
//

#include "../../include/dse/DSE.hpp"
#include "../../include/Definitions.hpp"

#include <cmath>
#include <iostream>
#include <iomanip>

DSE::DSE(double L2, MomentumGrid* p2Grid) : L2(L2), p2Grid(p2Grid)
{
    selfEnergy = new QuarkSelfEnergy(p2Grid);
    quarkPropagator = new QuarkPropagator(selfEnergy);

    gaussLegendreIntegrator = new GaussLegendre(200);
    gaussChebyshevIntegrator = new GaussChebyshev(32);  // TODO set gridpoints meaningfully
}

void DSE::solveDSE()
{
    std::cout << std::fixed;
    std::cout << std::setprecision(25);

    // Initiallize A(p2) = 1 and M(p2) = m
    //      --> Set Sigma_A and Sigma_M = 0

    for (int i = 0; i < p2Grid->getNumPoints(); i++)
    {
        selfEnergy->setSigmaAAt(i, 0);
        selfEnergy->setSigmaMAt(i, 0);
    }
    selfEnergy->reclacSigma_A_Spline();
    selfEnergy->reclacSigma_M_Spline();


    // For next iteration
    double next_Sigma_A_Vals[p2Grid->getNumPoints()];
    double next_Sigma_M_Vals[p2Grid->getNumPoints()];

    bool converged = false;
    int iteration = 1;
    do
    {
        // At each gridpoint: calculate Sigma_A and Sigma_M
        for (int i = 0; i < p2Grid->getNumPoints(); i++)
        {
            next_Sigma_A_Vals[i] = performIntegration_Sigma_A(p2Grid->momentumGridAtIdx(i));
            next_Sigma_M_Vals[i] = performIntegration_Sigma_M(p2Grid->momentumGridAtIdx(i));
        }


        // Check for Convergence
        converged = true;
        for (int i = 0; i < p2Grid->getNumPoints(); i++)
        {
            if((abs(next_Sigma_A_Vals[i] - selfEnergy->getSigmaAValGrid()[i]) > 1E-10) ||
               (abs(next_Sigma_M_Vals[i] - selfEnergy->getSigmaMValGrid()[i]) > 1E-10))
            {
                converged = false;
                break;
            }
        }


        // Update Grid and recalculate Splines
        for (int i = 0; i < p2Grid->getNumPoints(); i++)
        {
            selfEnergy->setSigmaAAt(i, next_Sigma_A_Vals[i]);
            selfEnergy->setSigmaMAt(i, next_Sigma_M_Vals[i]);
        }
        selfEnergy->reclacSigma_A_Spline();
        selfEnergy->reclacSigma_M_Spline();


        // Determine Self-Energies at the renorm point
        double Sigma_M_mu2 = selfEnergy->Sigma_M(RENORM_POINT2);
        double Sigma_A_mu2 = selfEnergy->Sigma_A(RENORM_POINT2);

        Z_2 = 1.0 - Sigma_A_mu2;


        std::cout << "Iteration " << iteration << ": Sigma_A(mu^2) = " << Sigma_A_mu2 << std::endl;
        std::cout << "Iteration " << iteration << ": Sigma_M(mu^2) = " << Sigma_M_mu2 << std::endl;
        iteration++;
    } while (!converged);
}

double DSE::selfEnergyIntegralKernelSigma_A(double p2, double q2, double z)
{
    double k2 = calc_k2(p2, q2, z);
    return quarkPropagator->sigma_v(q2) * g(k2) * F(p2, q2, k2, z);
}

double DSE::selfEnergyIntegralKernelSigma_M(double p2, double q2, double z)
{
    double k2 = calc_k2(p2, q2, z);
    return 3.0 * quarkPropagator->sigma_s(q2) * g(k2);
}

double DSE::F(double p2, double q2, double k2, double z)
{
    double q = sqrt(q2);
    double p = sqrt(p2);

    return (3.0 * q * z) / p - 2 * q2/k2 * (1 - pow(z ,2));
}

double DSE::g(double k2)
{
    return pow(Z_2, 2) * 16.0 * M_PI/3.0 * alpha(k2)/k2;
}

double DSE::alpha(double k2)
{
    double x = k2/LAMBDA;
    return M_PI * pow(ETA, 7) * pow(x, 2) * exp(-pow(ETA, 2) * x) + (2.0 * M_PI * GAMMA_m * (1.0 - exp(-k2/pow(LAMBDA_t, 2)))) /
                                                                    (log(pow(M_E, 2) - 1.0 + pow(1.0 + k2/pow(LAMBDA_QCD, 2), 2)));
}

double DSE::calc_k2(double p2, double q2, double z)
{
    double p = sqrt(p2);
    double q = sqrt(q2);
    return p2 + q2 - 2.0 * p * q * z;
}

double DSE::performIntegration_Sigma_M(double p2)
{
    std::function<double(double, double, double)> Sigma_A_integrand_function = [=, this](double q2, double p2, double z) -> double {
        return selfEnergyIntegralKernelSigma_M(p2, q2, z);
    };
    double integral_val = q2Integral(p2, Sigma_A_integrand_function, L2);
    return integral_val;
}

double DSE::performIntegration_Sigma_A(double p2)
{
    std::function<double(double, double, double)> Sigma_A_integrand_function = [=, this](double q2, double p2, double z) -> double {
        return selfEnergyIntegralKernelSigma_A(p2, q2, z);
    };
    double integral_val = q2Integral(p2, Sigma_A_integrand_function, L2);
    return integral_val;
}

double DSE::q2Integral(double p2, const std::function<double(double, double, double)> &f, double upper_cutoff)
{
    std::function<double(double)> q2Integrand = [=, this](double q2) -> double {
        return q2 * zIntegral(p2, q2, f);
    };
    return (4.0 * M_PI)/2.0 * 1.0/(pow(2.0 * M_PI, 4)) * gaussLegendreIntegrator->integrate(q2Integrand, 0, upper_cutoff);
}

double DSE::zIntegral(double p2, double q2, const std::function<double(double, double, double)> &f)
{
    std::function<double(double)> zIntegrand = [=, this](double z) -> double {
        return f(q2, p2, z);
    };
    return gaussChebyshevIntegrator->integrate_f_times_sqrt(zIntegrand);
}

QuarkPropagator *DSE::getQuarkPropagator() const
{
    return quarkPropagator;
}
