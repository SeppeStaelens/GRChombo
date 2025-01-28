/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(THINSHELLSOLUTION_HPP_)
#error "This file should only be included through ThinShellSolution.hpp"
#endif

#ifndef THINSHELLSOLUTION_IMPL_HPP_
#define THINSHELLSOLUTION_IMPL_HPP_

//#include "spline.h"
#include "FourthOrderNeville.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

ThinShellSolution::ThinShellSolution() {}

void ThinShellSolution::initialise_from_file()
{
    // Read in the data files with the solution from Fortran
    std::ifstream A_file("A.dat");
    std::ifstream m_file("m.dat");
    std::ifstream phi_file("Phi.dat");
    std::ifstream info_file("output.dat");

    if (!A_file.is_open() || !phi_file.is_open() || !m_file.is_open() ||
        !info_file.is_open())
    {
        std::cerr << "Error reading thinshell files!" << endl;
        exit(1);
    }

    // count the number of lines
    int line_count = 0;
    std::string line;
    while (std::getline(A_file, line))
    {
        line_count++;
    }
    A_file.clear();
    A_file.seekg(0, std::ios::beg);

    pout() << "The files contain " << line_count << " lines." << endl;

    // Read in the data
    std::string lineA, linePhi, linem, lineInfo;

    std::vector<double> A_vals(line_count);
    std::vector<double> phi_vals(line_count);
    std::vector<double> m_vals(line_count);
    std::vector<double> X_vals(line_count);
    std::vector<double> r_vals(line_count);

    double junk;
    int j = 0;
    while (std::getline(A_file, lineA) && j < line_count)
    {
        std::istringstream iss(lineA);
        if (iss >> r_vals[j] >> A_vals[j])
            j++;
    }
    j = 0;
    while (std::getline(m_file, linem) && j < line_count)
    {
        std::istringstream iss(linem);
        if (iss >> junk >> m_vals[j])
            j++;
    }

    double phi_offset, A_central, omega_pre_rescale, omega_file, Mass,
        compactness, r_99;

    // read in info values including omega, M, r_99, C
    while (std::getline(info_file, lineInfo))
    {
        std::istringstream iss(lineInfo);
        if (!lineInfo.empty() &&
            lineInfo[0] != '#') // ignore leading commented lines
        {
            if (!(iss >> A_central >> omega_pre_rescale >> junk >>
                  BS_frequency >> phi_offset >> Mass >> junk >> junk >>
                  compactness >> r_99 >> junk >> junk >> junk >> junk >> junk >>
                  junk >> junk))
                std::cout
                    << "WARNING: reading thinshell output.dat may have failed "
                    << endl;
        }
    }

    // read in phi-values, accounting for offset
    j = 0;
    while (std::getline(phi_file, linePhi) && j < line_count)
    {
        std::istringstream iss(linePhi);
        if (iss >> junk >> phi_vals[j])
        {
            phi_vals[j] -= phi_offset; // enforce phi(infty) = 0
            j++;
        }
    }

    X_vals[0] = 1.0;
    for (int i = 1; i < line_count; i++)
    {
        X_vals[i] = sqrt(r_vals[i] / (r_vals[i] - 2 * m_vals[i]));
    }

    pout() << "We now have central values A = " << A_vals[0]
           << ", X = " << X_vals[0] << ", phi = " << phi_vals[0]
           << ", r = " << r_vals[0] << endl;

    // // Set the derivatives to be zero at the boundaries
    // ASpline.set_boundary(tk::spline::first_deriv, 0.0, tk::spline::first_deriv,
    //                      0.0);
    // XSpline.set_boundary(tk::spline::first_deriv, 0.0, tk::spline::first_deriv,
    //                      0.0);
    // PhiSpline.set_boundary(tk::spline::first_deriv, 0.0,
    //                        tk::spline::first_deriv, 0.0);

    // // Set up the spline interpolators
    // PhiSpline.set_points(r_vals, phi_vals, tk::spline::cspline_hermite);
    // ASpline.set_points(r_vals, A_vals, tk::spline::cspline_hermite);
    // XSpline.set_points(r_vals, X_vals, tk::spline::cspline_hermite);

    // Set up the spline interpolators
    PhiSpline.set_data_points(&r_vals, &phi_vals);
    ASpline.set_data_points(&r_vals, &A_vals);
    XSpline.set_data_points(&r_vals, &X_vals);

    pout() << "TEST splines: phi, A, X at 0 are " << ASpline.interpolate(0.) << ", "
           << XSpline.interpolate(0.) << ", " << PhiSpline.interpolate(0.) << endl;

    // Asymptotic behaviour of different radii
    double max_iso_R = L;
    double max_arial_r = max_iso_R + Mass + Mass * Mass / (4 * max_iso_R);
    double f_max_arial_r = max_iso_R / max_arial_r;

    pout() << "The maximal isotropic radius, " << max_iso_R
           << ", corresponds to maximal areal radius " << max_arial_r
           << " giving f_at_max = " << f_max_arial_r << endl;

    double int_radius = 0.001, dr = 0.001, f_c = 1.;

    int len_int_array = floor(max_arial_r / dr) + 1;

    std::vector<double> f_vals(len_int_array);
    std::vector<double> iso_R_vals(len_int_array);
    std::vector<double> int_r_vals(len_int_array);

    j = 1;
    f_vals[0] = f_c;
    int_r_vals[0] = 0.;
    while (int_radius < max_arial_r)
    {
        int_r_vals[j] = int_radius;
        double k1 = (XSpline.interpolate(int_radius) - 1.) * f_c / int_radius;
        double k2 = (XSpline.interpolate(int_radius + dr / 2) - 1.) * (f_c + k1 * dr / 2) /
                    (int_radius + dr / 2);
        double k3 = (XSpline.interpolate(int_radius + dr / 2) - 1.) * (f_c + k2 * dr / 2) /
                    (int_radius + dr / 2);
        double k4 = (XSpline.interpolate(int_radius + dr) - 1.) * (f_c + k3 * dr) /
                    (int_radius + dr);
        f_c += dr * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
        f_vals[j] = f_c;
        int_radius += dr;
        j++;
    }

    double unscaled_f_at_max_r = f_vals[len_int_array - 1];
    pout() << "unscaled f at max r = " << unscaled_f_at_max_r << std::endl;
    j = 0;
    while (j < len_int_array)
    {
        f_vals[j] *= f_max_arial_r / unscaled_f_at_max_r;
        iso_R_vals[j] = f_vals[j] * int_r_vals[j];
        j++;
    }

    // r_from_R_Spline.set_points(iso_R_vals, int_r_vals,
    //                            tk::spline::cspline_hermite);
    // fSpline.set_points(int_r_vals, f_vals, tk::spline::cspline_hermite);
    r_from_R_Spline.set_data_points(&iso_R_vals, &int_r_vals);
    fSpline.set_data_points(&int_r_vals, &f_vals);

    // Calculate the aspect mass and the ADM mass at the boundary of the
    // physical domain. They shoud be similar, but not equal!
    aspect_mass = calculate_aspect_mass(max_iso_R);
    adm_mass = calculate_adm_mass(max_iso_R);
    radius = calculate_radius();
    compactness = aspect_mass / radius;

    pout() << "-----------------------------------------------" << endl;
    pout() << "Central Amplitude : " << A0 << endl;
    pout() << "ADM mass : " << adm_mass << endl;
    pout() << "Aspect mass : " << aspect_mass << endl;
    pout() << "Radius : " << radius << endl;
    pout() << "Compactness : " << compactness << endl;
    pout() << "w : " << BS_frequency << endl;
}

// Find the aspect mass
double ThinShellSolution::calculate_aspect_mass(double radius)
{
    double areal_radius = r_from_R_Spline.interpolate(radius);
    return 2. * radius * (sqrt(1. / fSpline.interpolate(areal_radius)) - 1.);
}

// Find the ADM mass
double ThinShellSolution::calculate_adm_mass(double radius)
{
    double areal_radius = r_from_R_Spline.interpolate(radius);
    double Xvalue = XSpline.interpolate(areal_radius);
    return -areal_radius * (1. / Xvalue - 1.) / fSpline.interpolate(areal_radius);
}

double ThinShellSolution::calculate_radius(double dx)
{
    double mass_99 = 0.99 * aspect_mass;
    double lower_radius = 0.;
    double enclosed_mass = 0.;
    while (enclosed_mass < mass_99)
    {
        lower_radius += dx;
        enclosed_mass = calculate_aspect_mass(lower_radius);
    }
    pout() << "radius = " << lower_radius << endl;

    return lower_radius;
}

double ThinShellSolution::get_BSfrequency() const { return BS_frequency; }

void ThinShellSolution::set_initialcondition_params(
    BosonStar_params_t m_params_BosonStar,
    Potential::params_t m_params_potential, const double max_r)
{

    A0 = m_params_BosonStar.central_amplitude_CSF;
    lambda = m_params_potential.phi4_coeff;
    solitonic = m_params_potential.solitonic;
    sigma = m_params_potential.sigma_soliton;
    L = max_r * 1.05; // just to make sure the function domain is slightly
                      // larger than the required cube
}

#endif /* THINSHELLSOLUTION_IMPL_HPP_ */
