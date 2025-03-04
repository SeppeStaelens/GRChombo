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
    std::vector<double> m_light_rings;
    std::vector<double> m_potential_extrema;
    bool found_light_rings = false;

    //! The constructor
    EffectivePotentialExtraction(const spherical_extraction_params_t &a_params,
                                 double a_dt, double a_time, bool a_first_step,
                                 double a_restart_time = 0.0)
        : SphericalExtraction(a_params, a_dt, a_time, a_first_step,
                              a_restart_time),
          m_params(a_params)
    {
        add_var(c_V_eff_n, VariableType::diagnostic);
        add_var(c_V_eff_d1, VariableType::diagnostic);
        add_var(c_V_eff_d2, VariableType::diagnostic);
    }

    //! Execute the query
    void execute_query(AMRInterpolator<Lagrange<4>> *a_interpolator,
                       std::string data_path)
    {
        extract(a_interpolator);

        if (m_params.write_extraction)
            write_extraction(data_path + m_params.extraction_file_prefix);

        std::vector<double> integrals_V_eff_n;
        std::vector<double> integrals_V_eff_d1;
        std::vector<double> integrals_V_eff_d2;

        add_var_integrand(0, integrals_V_eff_n);
        add_var_integrand(1, integrals_V_eff_d1);
        add_var_integrand(2, integrals_V_eff_d2);

        // do the integration over the surface
        integrate();

        std::vector<double> V_eff_values(integrals_V_eff_n.size());
        std::vector<double> V_eff_values2(integrals_V_eff_n.size());
	std::vector<double> areal_radii(integrals_V_eff_n.size());
        for (int i = 0; i < integrals_V_eff_n.size(); i++)
        {
            V_eff_values[i] =
                sqrt(integrals_V_eff_n[i] / integrals_V_eff_d1[i]);
            V_eff_values2[i] = sqrt(integrals_V_eff_n[i] / (4 * M_PI)) /
                               m_params.extraction_radii[i] /
                               sqrt(integrals_V_eff_d2[i] / (4 * M_PI));
	    areal_radii[i] = sqrt(integrals_V_eff_d2[i] / (4 * M_PI));
        }
        // write the integrals
        write_to_dat(V_eff_values, data_path, "EffectivePotential_");
        write_to_dat(V_eff_values2, data_path, "EffectivePotential2_");

        // find the light rings
        find_light_rings(V_eff_values, areal_radii);
        if (found_light_rings)
        {
            write_light_rings_to_dat(m_light_rings, data_path, "LightRings_");
            found_light_rings = false;
            m_light_rings.clear();
	    m_potential_extrema.clear();
        }
        find_light_rings(V_eff_values2, areal_radii);
        if (found_light_rings)
        {
            write_light_rings_to_dat(m_light_rings, data_path, "LightRings2_");
        }
    }

    void write_to_dat(std::vector<double> vals, std::string data_path,
                      std::string filename)
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

    //! This function estimates the extrema of the effective potential, which
    //! correspond to the light rings
    void find_light_rings(std::vector<double> vals, std::vector<double> areal)
    {
        std::vector<double> derivatives(vals.size());
        int num_radii = vals.size();
        double dr = m_params.extraction_radii[1] - m_params.extraction_radii[0];
        if (num_radii < 3)
        {
            pout() << "Need at least 3 radii to look for light rings" << endl;
        }
        else
        {
            for (int i = 1; i < num_radii - 1; i++)
            {
                derivatives[i] = (vals[i + 1] - vals[i - 1]) / (2 * dr);
            }
            derivatives[0] = (-3 * vals[0] + 4 * vals[1] - vals[2]) / (2 * dr);
            derivatives[num_radii - 1] =
                (3 * vals[num_radii - 1] - 4 * vals[num_radii - 2] +
                 vals[num_radii - 3]) /
                (2 * dr);

	    double a, b, r_i, r_extremum, r_ip1, K, r_areal_ex;
            for (int i = 0; i < num_radii - 1; i++)
            {
                if (derivatives[i] == 0)
                {
                    found_light_rings = true;
                    m_light_rings.push_back(m_params.extraction_radii[i]);
		    m_light_rings.push_back(areal[i]);
                }
                else if (derivatives[i] * derivatives[i + 1] < 0)
                {
                    found_light_rings = true;
		    r_i = m_params.extraction_radii[i];
  		    a = (derivatives[i+1] - derivatives[i]) / (2 * dr);
		    b = derivatives[i] - 2 * a * r_i;
		    r_extremum = -b / (2*a);
		    m_light_rings.push_back(r_extremum);
		    m_potential_extrema.push_back(vals[i] + a*(r_extremum*r_extremum - r_i*r_i) + b*(r_extremum - r_i)); 
                
		    r_ip1 = m_params.extraction_radii[i+1];
		    K = (areal[i+1] - areal[i]) / (r_ip1 - r_i);
		    r_areal_ex = areal[i] + K * (r_extremum - r_i);
		    m_light_rings.push_back(r_areal_ex);
		}
            }
        }
    }

    void write_light_rings_to_dat(std::vector<double> vals,
                                  std::string data_path, std::string filename)
    {
        std::vector<string> title_line(6);
        title_line[0] = "nr_of_LRs";
        title_line[1] = "r1";
	title_line[2] = "r1_areal";
        title_line[3] = "r2";
	title_line[4] = "r2_areal";
	title_line[5] = "delta_V";

        SmallDataIO file(data_path + filename, m_dt, m_time, m_restart_time,
                         SmallDataIO::APPEND, m_first_step);

        if (m_time > 0)
            file.remove_duplicate_time_data();

        if (m_time == 0.)
        {
            file.write_header_line(title_line);
        }
        vals.insert(vals.begin(), vals.size());
	double potential_diff = 0;
	for (int i = 0; i < floor(m_potential_extrema.size() / 2); i++)
        {
            potential_diff = abs(m_potential_extrema[2*i] - m_potential_extrema[2*i+1]);
	    vals.insert(vals.end(), potential_diff);
	}
	file.write_time_data_line(vals);
    }

    ~EffectivePotentialExtraction() { ; }
};

#endif /* EFFECTIVEPOTENTIALEXTRACTION_HPP_ */
