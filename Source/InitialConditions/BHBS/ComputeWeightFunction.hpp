#ifndef COMPUTEWEIGHTFUNCTION_HPP_
#define COMPUTEWEIGHTFUNCTION_HPP_

// Chombo includes
#include "IntVect.H"

// #include "simd.hpp"
#include "VarsTools.hpp"
#include <array>
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DebuggingTools.hpp"
#include "BosonStarParams.hpp"
#include "WeightFunction.hpp"
#include "simd.hpp"

class ComputeWeightFunction
{
    protected:

    BosonStar_params_t m_params_BosonStar;
    BlackHole_params_t m_params_BlackHole;
    double m_dx;

    public:

    template <class data_t> struct weightfunc_t
    {
        data_t weight1;
        data_t weight2;
    };

    ComputeWeightFunction(BosonStar_params_t a_params_BosonStar,
                   BlackHole_params_t a_params_BlackHole, 
                   const double a_dx)
                    : m_dx(a_dx), 
                    m_params_BosonStar(a_params_BosonStar), 
                    m_params_BlackHole(a_params_BlackHole)
    {   
    }

    //! Function to compute the value of the initial vars on the grid
    template <class data_t>
    void compute(Cell<data_t> current_cell) const
    {
        weightfunc_t<data_t> out;

        Coordinates<double> coords(current_cell, m_dx,
        m_params_BosonStar.star_centre);
	
	    double separation = m_params_BosonStar.binary_separation;
    	double impact_parameter = m_params_BosonStar.BS_impact_parameter;
	    double q = m_params_BosonStar.mass_ratio;
        double rapidity = m_params_BosonStar.BS_rapidity;
        double rapidity2 = m_params_BlackHole.BH_rapidity;
        double radius_width1 = m_params_BosonStar.radius_width1;
        double radius_width2 = m_params_BosonStar.radius_width2;

	double x_p2 = (separation) * cosh(rapidity);
        double z_p2 = 0.; //set /tilde{t} to zero
        double y_p2 = -impact_parameter;

        WeightFunction weightfunction(separation, x_p2, y_p2, z_p2, m_params_BlackHole.weight_function_order);

        double profile_func1 = weightfunction.profile_chi((coords.x-q*separation/(q+1))*cosh(rapidity), coords.y+q*impact_parameter/(q+1.), coords.z);

        current_cell.store_vars(profile_func1, c_profile1);
    
    }

    
};

#endif /* COMPUTEWEIGHTFUNCTION_HPP_ */
