/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef CONSTRAINTVIOLATIONS_HPP_
#define CONSTRAINTVIOLATIONS_HPP_

#include "GRAMR.hpp"
#include "Interval.H"
#include "SmallDataIO.hpp"
#include <sstream>

//! This class calculates the p-norm of the Hamiltonian and Momentum
//! constraints and writes them to a file
class ConstraintViolations
{
public:
    // Constructor
    ConstraintViolations(const int a_Ham_comp, const Interval a_Mom_comps,
                         GRAMR *a_gr_amr, const double a_dx_coarse,
                         const double a_dt, const double a_time,
                         const double a_restart_time,
                         const std::string a_filename,
                         const bool a_first_step,
                         const double a_norm_exponent = 2.)
        : m_Ham_comp(a_Ham_comp), m_Mom_comps(a_Mom_comps), m_gr_amr(a_gr_amr),
        m_dx_coarse(a_dx_coarse), m_dt(a_dt), m_time(a_time),
        m_restart_time(a_restart_time), m_filename(a_filename),
        m_first_step(a_first_step), m_norm_exponent(a_norm_exponent) {}

    // Calculates norms and writes to file
    void execute();

    // Returns norms
    std::pair<double, double> get_norms();

private:
    const int m_Ham_comp;
    const Interval m_Mom_comps;
    GRAMR *m_gr_amr;
    const double m_dx_coarse;
    const double m_dt;
    const double m_time;
    const double m_restart_time;
    const std::string m_filename;
    const bool m_first_step;
    const double m_norm_exponent;
    double m_Ham_norm;
    double m_Mom_norm;

    // Calculate norms
    void compute_norms();

    // Write norms to file
    void write_norms() const;
};

#include "ConstraintViolations.impl.hpp"

#endif /* CONSTRAINTVIOLATIONS_HPP_ */
