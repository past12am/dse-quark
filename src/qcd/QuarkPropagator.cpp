//
// Created by past12am on 04/06/24.
//

#include "../../include/qcd/QuarkPropagator.hpp"
#include "../../include/Definitions.hpp"

#include <cmath>
#include <gsl/gsl_complex_math.h>

QuarkPropagator::QuarkPropagator(QuarkSelfEnergy* selfEnergy) : selfEnergy(selfEnergy)
{

}

gsl_complex QuarkPropagator::sigma_v(gsl_complex p2, InterpAccels* interpAccels)
{
    return gsl_complex_div(GSL_COMPLEX_ONE, (gsl_complex_mul(A(p2, interpAccels), (gsl_complex_add(p2, gsl_complex_pow_real(M(p2, interpAccels), 2))))));
}

gsl_complex QuarkPropagator::sigma_s(gsl_complex p2, InterpAccels* interpAccels)
{
    return gsl_complex_mul(sigma_v(p2, interpAccels), M(p2, interpAccels));
}

gsl_complex QuarkPropagator::M(gsl_complex p2, InterpAccels* interpAccels)
{
    return gsl_complex_mul(gsl_complex_div(GSL_COMPLEX_ONE, A(p2, interpAccels)), (gsl_complex_add_real(gsl_complex_sub(selfEnergy->Sigma_M(p2, interpAccels), selfEnergy->Sigma_M(gsl_complex_rect(RENORM_POINT2, 0), interpAccels)), QUARK_MASS)));
}

gsl_complex QuarkPropagator::A(gsl_complex p2, InterpAccels* interpAccels)
{
    return gsl_complex_add_real(gsl_complex_sub(selfEnergy->Sigma_A(p2, interpAccels), selfEnergy->Sigma_A(gsl_complex_rect(RENORM_POINT2, 0), interpAccels)), 1.0);
}


gsl_complex QuarkPropagator::sigma_v(gsl_complex p2)
{
    return gsl_complex_div(GSL_COMPLEX_ONE, (gsl_complex_mul(A(p2), (gsl_complex_add(p2, gsl_complex_pow_real(M(p2), 2))))));
}

gsl_complex QuarkPropagator::sigma_s(gsl_complex p2)
{

    return gsl_complex_mul(sigma_v(p2), M(p2));
}

gsl_complex QuarkPropagator::M(gsl_complex p2)
{
    return gsl_complex_mul(gsl_complex_div(GSL_COMPLEX_ONE, A(p2)), (gsl_complex_add_real(gsl_complex_sub(selfEnergy->Sigma_M(p2), selfEnergy->Sigma_M(gsl_complex_rect(RENORM_POINT2, 0))), QUARK_MASS)));
}

gsl_complex QuarkPropagator::A(gsl_complex p2)
{
    return gsl_complex_add_real(gsl_complex_sub(selfEnergy->Sigma_A(p2), selfEnergy->Sigma_A(gsl_complex_rect(RENORM_POINT2, 0))), 1.0);
}