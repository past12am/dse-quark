//
// Created by past12am on 04/06/24.
//

#include "../../include/dse/DSE.hpp"
#include "../../include/Definitions.hpp"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <future>
#include <gsl/gsl_complex_math.h>

DSE::DSE(double L2, ComplexMomentumGrid* p2Grid) : L2(L2), p2Grid(p2Grid)
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

    for (int j = 0; j < p2Grid->getNumImagPoints(); j++)
    {
        for (int i = 0; i < p2Grid->getNumRealPoints(); i++)
        {
            selfEnergy->setSigmaAAt(j, i, GSL_COMPLEX_ZERO);
            selfEnergy->setSigmaMAt(j, i, GSL_COMPLEX_ZERO);
        }
    }
    selfEnergy->reclacSigma_A_Spline();
    selfEnergy->reclacSigma_M_Spline();


    // For next iteration
    gsl_complex next_Sigma_A_Vals[p2Grid->getNumPoints()];
    gsl_complex next_Sigma_M_Vals[p2Grid->getNumPoints()];

    // For Threadpool
    std::future<gsl_complex> future_newSigma_A_at_p2_val[NUM_THREADS];
    std::future<gsl_complex> future_newSigma_M_at_p2_val[NUM_THREADS];

    bool converged = false;
    int iteration = 1;
    do
    {
        // At each gridpoint: calculate Sigma_A and Sigma_M
        for (int j = 0; j < p2Grid->getNumImagPoints(); j++)
        {
            for (int i = 0; i < p2Grid->getNumRealPoints(); i++)
            {
                // Calculate Sigma_A(p2) and Sigma_M(p2) and set grid values

                int allocation_counter = 0;
                for (int threadIdx = 0; threadIdx < NUM_THREADS; threadIdx++)
                {
                    future_newSigma_A_at_p2_val[threadIdx] = std::async([this, j, i]() -> gsl_complex {
                        InterpAccels interpAccels = allocInterpAccels();

                        gsl_complex res = performIntegration_Sigma_A(p2Grid->momentumGridAtIdx(j, i), &interpAccels);

                        freeInterpAccels(interpAccels);
                        return res;
                    });

                    future_newSigma_M_at_p2_val[threadIdx] = std::async([this, j, i]() -> gsl_complex {
                        InterpAccels interpAccels = allocInterpAccels();

                        gsl_complex res = performIntegration_Sigma_M(p2Grid->momentumGridAtIdx(j, i), &interpAccels);

                        freeInterpAccels(interpAccels);

                        return res;
                    });

                    // Increase i, as we allocated it for threads
                    i++;
                    allocation_counter++;

                    // If we hit the last possible i, break
                    if(i == p2Grid->getNumRealPoints())
                        break;
                }


                // Note: that i changed by allocation_counter, thus
                //      next_Sigma_A_Vals[i - allocation_counter + threadIdx]
                for (int threadIdx = 0; threadIdx < NUM_THREADS; threadIdx++)
                {
                    next_Sigma_A_Vals[i - allocation_counter + threadIdx] = future_newSigma_A_at_p2_val[threadIdx].get();
                    next_Sigma_M_Vals[i - allocation_counter + threadIdx] = future_newSigma_M_at_p2_val[threadIdx].get();

                    if(i - allocation_counter + threadIdx == p2Grid->getNumRealPoints() - 1)
                        break;
                }

                // As i will be increased in the next step (although already has been increased) --> correct
                i--;
            }
        }


        // Check for Convergence
        converged = true;
        for (int j = 0; j < p2Grid->getNumImagPoints(); j++)
        {
            for (int i = 0; i < p2Grid->getNumRealPoints(); i++)
            {
                if ((gsl_complex_abs(gsl_complex_sub(next_Sigma_A_Vals[i], selfEnergy->getSigmaAValGrid(j, i))) > 1E-10) ||
                    (gsl_complex_abs(gsl_complex_sub(next_Sigma_M_Vals[i], selfEnergy->getSigmaMValGrid(j, i))) > 1E-10))
                {
                    converged = false;
                    break;
                }
            }
        }


        // Update Grid and recalculate Splines
        for (int j = 0; j < p2Grid->getNumImagPoints(); j++)
        {
            for (int i = 0; i < p2Grid->getNumRealPoints(); i++)
            {
                selfEnergy->setSigmaAAt(j, i, next_Sigma_A_Vals[i]);
                selfEnergy->setSigmaMAt(j, i, next_Sigma_M_Vals[i]);
            }
        }
        selfEnergy->reclacSigma_A_Spline();
        selfEnergy->reclacSigma_M_Spline();


        // Determine Self-Energies at the renorm point
        gsl_complex Sigma_M_mu2 = selfEnergy->Sigma_M(gsl_complex_rect(RENORM_POINT2, 0));
        gsl_complex Sigma_A_mu2 = selfEnergy->Sigma_A(gsl_complex_rect(RENORM_POINT2, 0));

        Z_2 = gsl_complex_add_real(gsl_complex_negative(Sigma_A_mu2), 1.0);


        std::cout << "Iteration " << iteration << ": Sigma_A(mu^2) = " << GSL_REAL(Sigma_A_mu2) << ((GSL_IMAG(Sigma_A_mu2) < 0) ? " - " : " + ") << GSL_IMAG(Sigma_A_mu2) << "i" << std::endl;
        std::cout << "Iteration " << iteration << ": Sigma_M(mu^2) = " << GSL_REAL(Sigma_M_mu2) << ((GSL_IMAG(Sigma_M_mu2) < 0) ? " - " : " + ") << GSL_IMAG(Sigma_M_mu2) << "i" << std::endl;
        iteration++;
    } while (!converged);
}

gsl_complex DSE::selfEnergyIntegralKernelSigma_A(gsl_complex p2, gsl_complex q2, double z, InterpAccels* interpAccels)
{
    gsl_complex k2 = calc_k2(p2, q2, z);
    return gsl_complex_mul(gsl_complex_mul(quarkPropagator->sigma_v(q2, interpAccels), g(k2)), F(p2, q2, k2, z));
}

gsl_complex DSE::selfEnergyIntegralKernelSigma_M(gsl_complex p2, gsl_complex q2, double z, InterpAccels* interpAccels)
{
    gsl_complex k2 = calc_k2(p2, q2, z);
    return gsl_complex_mul_real(gsl_complex_mul(quarkPropagator->sigma_s(q2, interpAccels), g(k2)), 3.0);
}

gsl_complex DSE::F(gsl_complex p2, gsl_complex q2, gsl_complex k2, double z)
{
    gsl_complex q = gsl_complex_sqrt(q2);
    gsl_complex p = gsl_complex_sqrt(p2);

    return gsl_complex_sub(gsl_complex_div(gsl_complex_mul_real(q, 3.0 * z), p), gsl_complex_mul_real(gsl_complex_div(q2, k2), 2.0 * (1 - pow(z ,2))));
}

gsl_complex DSE::g(gsl_complex k2)
{
    return gsl_complex_mul(gsl_complex_pow_real(Z_2, 2), gsl_complex_mul_real(gsl_complex_div(alpha(k2), k2), 16.0 * M_PI/3.0));
}

gsl_complex DSE::alpha(gsl_complex k2)
{
    gsl_complex x = gsl_complex_div_real(k2, LAMBDA);
    return gsl_complex_add(gsl_complex_mul(gsl_complex_mul_real(gsl_complex_pow_real(x, 2), M_PI * pow(ETA, 7)), gsl_complex_exp(gsl_complex_mul_real(x, -pow(ETA, 2)))),
                           gsl_complex_div(gsl_complex_mul_real(gsl_complex_sub(GSL_COMPLEX_ONE, gsl_complex_exp(gsl_complex_div_real(k2, -pow(LAMBDA_t, 2)))), 2.0 * M_PI * GAMMA_m),
                                              gsl_complex_log(gsl_complex_add_real(gsl_complex_pow_real(gsl_complex_add_real(gsl_complex_div_real(k2, pow(LAMBDA_QCD, 2)), 1), 2), pow(M_E, 2) - 1.0))));
}

gsl_complex DSE::calc_k2(gsl_complex p2, gsl_complex q2, double z)
{
    gsl_complex p = gsl_complex_sqrt(p2);
    gsl_complex q = gsl_complex_sqrt(q2);
    return gsl_complex_sub(gsl_complex_add(p2, q2), gsl_complex_mul_real(gsl_complex_mul(p, q), 2.0 * z));
}

gsl_complex DSE::performIntegration_Sigma_M(gsl_complex p2, InterpAccels* interpAccels)
{
    std::function<gsl_complex(gsl_complex, gsl_complex, double)> Sigma_M_integrand_function = [=, this](gsl_complex q2, gsl_complex p2, double z) -> gsl_complex {
        return selfEnergyIntegralKernelSigma_M(p2, q2, z, interpAccels);
    };
    gsl_complex integral_val = q2Integral(p2, Sigma_M_integrand_function, L2);
    return integral_val;
}

gsl_complex DSE::performIntegration_Sigma_A(gsl_complex p2, InterpAccels* interpAccels)
{
    std::function<gsl_complex(gsl_complex, gsl_complex, double)> Sigma_A_integrand_function = [=, this](gsl_complex q2, gsl_complex p2, double z) -> gsl_complex {
        return selfEnergyIntegralKernelSigma_A(p2, q2, z, interpAccels);
    };
    gsl_complex integral_val = q2Integral(p2, Sigma_A_integrand_function, L2);
    return integral_val;
}

gsl_complex DSE::q2Integral(gsl_complex p2, const std::function<gsl_complex(gsl_complex, gsl_complex, double)> &f, double upper_cutoff)
{
    std::function<gsl_complex(gsl_complex)> q2Integrand = [=, this](gsl_complex q2) -> gsl_complex {
        return gsl_complex_mul(q2, zIntegral(p2, q2, f));
    };
    return gsl_complex_mul_real(gaussLegendreIntegrator->integrate(q2Integrand, 0, upper_cutoff), (4.0 * M_PI)/2.0 * 1.0/(pow(2.0 * M_PI, 4)));
}

gsl_complex DSE::zIntegral(gsl_complex p2, gsl_complex q2, const std::function<gsl_complex(gsl_complex , gsl_complex , double)> &f)
{
    std::function<gsl_complex(double)> zIntegrand = [=, this](double z) -> gsl_complex {
        return f(q2, p2, z);
    };
    return gaussChebyshevIntegrator->integrate_f_times_sqrt(zIntegrand);
}

QuarkPropagator *DSE::getQuarkPropagator() const
{
    return quarkPropagator;
}

InterpAccels DSE::allocInterpAccels()
{
    return InterpAccels{
        .interpAccel_M_X_real = gsl_interp_accel_alloc(),
        .interpAccel_M_Y_real = gsl_interp_accel_alloc(),
        .interpAccel_A_X_real = gsl_interp_accel_alloc(),
        .interpAccel_A_Y_real = gsl_interp_accel_alloc(),
        .interpAccel_M_X_imag = gsl_interp_accel_alloc(),
        .interpAccel_M_Y_imag = gsl_interp_accel_alloc(),
        .interpAccel_A_X_imag = gsl_interp_accel_alloc(),
        .interpAccel_A_Y_imag = gsl_interp_accel_alloc(),
    };
}

void DSE::freeInterpAccels(InterpAccels interpAccels)
{
    gsl_interp_accel_free(interpAccels.interpAccel_M_X_real);
    gsl_interp_accel_free(interpAccels.interpAccel_M_Y_real);
    gsl_interp_accel_free(interpAccels.interpAccel_A_X_real);
    gsl_interp_accel_free(interpAccels.interpAccel_A_Y_real);
    gsl_interp_accel_free(interpAccels.interpAccel_M_X_imag);
    gsl_interp_accel_free(interpAccels.interpAccel_M_Y_imag);
    gsl_interp_accel_free(interpAccels.interpAccel_A_X_imag);
    gsl_interp_accel_free(interpAccels.interpAccel_A_Y_imag);
}
