
#ifndef THINSHELLSOLUTION_HPP_
#define THINSHELLSOLUTION_HPP_

#include "cmath"
#include "spline.h"
#include <vector>

class ThinShellSolution
{

  private:          // private member variables/arrays
    double A0;      // Klein Gordon mass squared, KG scalr field central
                    // amplitude
    double lambda;  // phi 4 coupling in Klein gordon potential
    double sigma;   // 0.2 works with PC = 0.05 // parameter for solitonic stars
    bool solitonic; // false fro mini/lambda star. true for solitonic star
    double L;       // L, length of domain, dx.
    double BS_frequency; // frequency of the boson star
    double aspect_mass;
    double adm_mass;
    double radius;
    double compactness;

    void calculate_aspect_mass();
    void calculate_adm_mass();
    double calculate_radius();

  public:
    ThinShellSolution();

    void set_initialcondition_params(BosonStar_params_t m_params_BosonStar,
                                     Potential::params_t m_params_potential,
                                     const double max_r);
    void initialise_from_file();
    double get_BSfrequency() const;
    void output_csv();

    tk::spline ASpline, PhiSpline, fSpline, r_from_R_Spline;
};

#include "ThinShellSolution.impl.hpp"

#endif /* THINSHELLSOLUTION_HPP_ */
