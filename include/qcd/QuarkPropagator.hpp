//
// Created by past12am on 04/06/24.
//

#ifndef QUARKDSE_V3_QUARKPROPAGATOR_HPP
#define QUARKDSE_V3_QUARKPROPAGATOR_HPP



#include "QuarkSelfEnergy.hpp"



class QuarkPropagator
{
    private:
        QuarkSelfEnergy* selfEnergy;

    public:
        QuarkPropagator(QuarkSelfEnergy* selfEnergy);

        gsl_complex M(gsl_complex p2);
        gsl_complex A(gsl_complex p2);
        gsl_complex M(gsl_complex p2, InterpAccels* interpAccels);
        gsl_complex A(gsl_complex p2, InterpAccels* interpAccels);

        gsl_complex sigma_v(gsl_complex p2);
        gsl_complex sigma_s(gsl_complex p2);
        gsl_complex sigma_v(gsl_complex p2, InterpAccels* interpAccels);
        gsl_complex sigma_s(gsl_complex p2, InterpAccels* interpAccels);
};



#endif //QUARKDSE_V3_QUARKPROPAGATOR_HPP
