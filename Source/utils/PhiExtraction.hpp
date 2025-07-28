/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef PHIEXTRACTION_HPP_
#define PHIEXTRACTION_HPP_

#include "SphericalExtraction.hpp"

//!  The class allows extraction of the values of the scalar field on
//!  spherical shells at specified radii, and integration over those shells
/*!
   The class allows the user to extract data from the grid for the scalar field
   components over spherical shells at specified radii. The values may then be
   written to an output file, or integrated across the surfaces.
*/
class PhiExtraction : public SphericalExtraction
{
  public:
    //! The constructor
    PhiExtraction(const spherical_extraction_params_t &a_params, double a_dt,
                   double a_time, bool a_first_step,
                   double a_restart_time = 0.0)
        : SphericalExtraction(a_params, a_dt, a_time, a_first_step,
                              a_restart_time)
    {
        add_var(c_phi_Re, VariableType::evolution);
        add_var(c_phi_Im, VariableType::evolution);
    }

    //! The old constructor which assumes it is called in specificPostTimeStep
    //! so the first time step is when m_time == m_dt
    PhiExtraction(const spherical_extraction_params_t &a_params, double a_dt,
                   double a_time, double a_restart_time = 0.0)
        : PhiExtraction(a_params, a_dt, a_time, (a_dt == a_time),
                         a_restart_time)
    {
    }

    //! Execute the query
    void execute_query(AMRInterpolator<Lagrange<4>> *a_interpolator)
    {
        // extract the values of the Phi scalars on the spheres
        extract(a_interpolator);

        if (m_params.write_extraction)
            write_extraction(m_params.extraction_file_prefix);

        // now calculate and write the requested spherical harmonic modes
        std::vector<std::pair<std::vector<double>, std::vector<double>>>
            mode_integrals(m_num_modes);

        // note that this is normalised by multiplying by radius
        auto normalised_Phi_complex =
            [](std::vector<double> Phi_reim_parts, double r, double, double)
        {
            // here the std::vector<double> passed will just have
            // the real and imaginary parts of the scalar field as its
            // only components
            return std::make_pair(Phi_reim_parts[0],
                   Phi_reim_parts[1]);
        };

        // add the modes that will be integrated
        for (int imode = 0; imode < m_num_modes; ++imode)
        {
            const auto &mode = m_modes[imode];
            constexpr int es = 0;
            add_mode_integrand(es, mode.first, mode.second,
                               normalised_Phi_complex, mode_integrals[imode]);
        }

        // do the integration over the surface
        integrate();

        // write the integrals
        for (int imode = 0; imode < m_num_modes; ++imode)
        {
            const auto &mode = m_modes[imode];
            std::string integrals_filename = m_params.integral_file_prefix +
                                             std::to_string(mode.first) +
                                             std::to_string(mode.second);
            std::vector<std::vector<double>> integrals_for_writing = {
                std::move(mode_integrals[imode].first),
                std::move(mode_integrals[imode].second)};
            std::vector<std::string> labels = {"integral Re", "integral Im"};
	    write_integrals(integrals_filename, integrals_for_writing, labels);
        }
    }
};

#endif /* PHIEXTRACTION_HPP_ */
