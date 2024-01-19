/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(CONSTRAINTVIOLATIONS_HPP_)
#error "This file should only be included through ConstraintViolations.hpp"
#endif

#ifndef CONSTRAINTVIOLATIONS_IMPL_HPP_
#define CONSTRAINTVIOLATIONS_IMPL_HPP_

void ConstraintViolations::execute()
{
    compute_norms();
    write_norms();
}

std::pair<double, double> ConstraintViolations::get_norms()
{
    return std::make_pair(m_Ham_norm, m_Mom_norm);
}

void ConstraintViolations::compute_norms()
{
    m_Ham_norm = m_gr_amr->compute_norm(Interval(m_Ham_comp, m_Ham_comp),
                                      m_norm_exponent, m_dx_coarse);
    m_Mom_norm = m_gr_amr->compute_norm(m_Mom_comps, m_norm_exponent,m_dx_coarse);
}

void ConstraintViolations::write_norms() const
{
    CH_TIME("ConstraintViolations::write_norms");

    SmallDataIO file(m_filename, m_dt, m_time, m_restart_time,
                     SmallDataIO::APPEND, m_first_step);

    // remove any duplicate data if this is a restart
    // note that this only does something if this is the first timestep after
    // a restart
    file.remove_duplicate_time_data();

    // need to write headers if this is the first timestep
    if (m_first_step)
    {
        std::stringstream norm_exponent_ss;
        norm_exponent_ss << m_norm_exponent;
        std::string Ham_header_string = "L^" + norm_exponent_ss.str() + "_Ham";
        std::string Mom_header_string = "L^" + norm_exponent_ss.str() + "_Mom";
        // write headers
        file.write_header_line({Ham_header_string, Mom_header_string});
    }

    // write data
    file.write_time_data_line({m_Ham_norm, m_Mom_norm});
}

#endif /* CONSTRAINTVIOLATIONS_IMPL_HPP_ */
