/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef EFFECTIVEPOTENTIALEXTRACTION_HPP_
#define EFFECTIVEPOTENTIALEXTRACTION_HPP_

#include "SmallDataIO.hpp" // for writing data
#include "SphericalExtraction.hpp"
#include "UserVariables.hpp" //This files needs c_NUM - total number of components

/*!
   The class allows the user to extract data from the grid for the
   effective potential for null geodesics over spherical shells at specified
   radii. The values may then be written to an output file, or integrated across
   the surfaces.
*/

class EffectivePotentialExtraction : public SphericalExtraction
{
  public:
    string m_filename = "EffectivePotential";

    //! The constructor
    EffectivePotentialExtraction(const spherical_extraction_params_t &a_params,
                                 double a_dt, double a_time, bool a_first_step,
                                 double a_restart_time = 0.0)
        : SphericalExtraction(a_params, a_dt, a_time, a_first_step,
                              a_restart_time)
    {
        add_var(c_V_eff, VariableType::diagnostic);
    }

    //! Execute the query
    void execute_query(AMRInterpolator<Lagrange<4>> *a_interpolator)
    {
        // extract the values of the Weyl scalars on the spheres
        extract(a_interpolator);

        if (m_params.write_extraction)
            write_extraction(m_params.extraction_file_prefix);

        std::vector<double> integrals;

        // the integrand lambda function
        auto integrand = [](std::vector<double> effective_potential_vals,
                            double r, double, double)
        { return effective_potential_vals[0]; };

        add_integrand(integrand, integrals);

        // do the integration over the surface
        integrate();

        // write the integrals
        write_to_dat(integrals);
    }

    void write_to_dat(std::vector<double> vals)
    {

        std::vector<string> title_line(m_params.num_extraction_radii);
        string dummy_string;
        for (int i = 0; i < m_params.num_extraction_radii; i++)
        {
            dummy_string = "r = " + to_string(m_params.extraction_radii[i]);
            title_line[i] = dummy_string;
        }

        SmallDataIO potential_file(m_filename, m_dt, m_time, m_restart_time,
                                   SmallDataIO::APPEND, m_first_step);

        if (m_time > 0)
            potential_file.remove_duplicate_time_data();

        if (m_time == 0.)
        {
            potential_file.write_header_line(title_line);
        }

        potential_file.write_time_data_line(vals);
    }

    ~EffectivePotentialExtraction() {}
};

#endif /* EFFECTIVEPOTENTIALEXTRACTION_HPP_ */
