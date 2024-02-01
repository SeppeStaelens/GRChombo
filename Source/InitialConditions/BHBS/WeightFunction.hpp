#ifndef WEIGHTFUNCTION_HPP_
#define WEIGHTFUNCTION_HPP_

#include <array>

class WeightFunction
{
    public:

        double m_separation;
        double m_x_BS, m_y_BS, m_z_BS;
        int m_order;

        WeightFunction(double separation, double x_BS, double y_BS, double z_BS , int n){
            m_separation = separation;
            m_x_BS = x_BS;
            m_y_BS = y_BS;
            m_z_BS = z_BS;
            m_order = n;
        }

        double profile_chi(double coord_x, double coord_y, double coord_z)
        {
            double dist_BS = sqrt(m_x_BS*m_x_BS + m_y_BS*m_y_BS + m_z_BS*m_z_BS);
            double numer = m_separation * pow(m_separation - dist_BS, m_order);
            double denom = pow(m_separation, m_order + 1) + pow(dist_BS, m_order + 1);
            return numer / denom;
        }

        double profile_chi_2(double coord_x, double coord_y, double coord_z)
        {
            
        }
};

#endif /* WEIGHTFUNCTION_HPP_ */
