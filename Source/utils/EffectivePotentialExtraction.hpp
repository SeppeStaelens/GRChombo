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
    spherical_extraction_params_t m_params;

    //! The constructor
    EffectivePotentialExtraction(const spherical_extraction_params_t &a_params,
                                 double a_dt, double a_time, bool a_first_step,
                                 double a_restart_time = 0.0)
        : SphericalExtraction(a_params, a_dt, a_time, a_first_step,
                              a_restart_time), m_params(a_params)
    {
        add_var(c_V_eff, VariableType::diagnostic);
	add_var(c_unit, VariableType::diagnostic);
	//std::vector<int> metric_vars = {c_chi, c_h11, c_h12, c_h13, c_h22, c_h23, c_h33};
	//add_evolution_vars(metric_vars);
	add_var(c_volume, VariableType::diagnostic);
    }

    //! Execute the query
    void execute_query(AMRInterpolator<Lagrange<4>> *a_interpolator, std::string data_path)
    {
        extract(a_interpolator);

        if (m_params.write_extraction)
            write_extraction(data_path + m_params.extraction_file_prefix);

        std::vector<double> integrals_V_eff;
	std::vector<double> integrals_unit;
	std::vector<double> integrals_volume;

        // the integrand lambda function
        //auto integrand = [](std::vector<double> vals,
        //                    double r, double theta, double phi)
        //{ 
	//	double g_th_th = 
	//	return vals[0]; 
	//};

        add_var_integrand(0, integrals_V_eff);
	add_var_integrand(1, integrals_unit);
	add_var_integrand(2, integrals_volume);
	//add_integrand(integrand, integrals_unit)

        // do the integration over the surface
        integrate();
	
	std::vector<double> values(integrals_V_eff.size());
	std::vector<double> values2(integrals_V_eff.size());
	for (int i = 0; i < integrals_V_eff.size(); i++)
	{
		values[i] = sqrt(integrals_V_eff[i] / integrals_unit[i]);
		values2[i] = sqrt(integrals_V_eff[i] / (4*M_PI))/m_params.extraction_radii[i] / sqrt(integrals_volume[i] / (4*M_PI));
	}
        // write the integrals
        write_to_dat(values, data_path, "EffectivePotential_");
	write_to_dat(values2, data_path, "EffectivePotential2_");
    }

    void write_to_dat(std::vector<double> vals, std::string data_path, std::string filename)
    {
	
        std::vector<string> title_line(m_params.num_extraction_radii);
        string dummy_string;
        for (int i = 0; i < m_params.num_extraction_radii; i++)
        {
            dummy_string = "r = " + to_string(m_params.extraction_radii[i]);
            title_line[i] = dummy_string;
        }
	
        SmallDataIO file(data_path + filename, m_dt, m_time, m_restart_time,
                                   SmallDataIO::APPEND, m_first_step);

        if (m_time > 0)
            file.remove_duplicate_time_data();

        if (m_time == 0.)
        {
            file.write_header_line(title_line);
        }

        file.write_time_data_line(vals);
    }

    ~EffectivePotentialExtraction() {;}
};

#endif /* EFFECTIVEPOTENTIALEXTRACTION_HPP_ */
