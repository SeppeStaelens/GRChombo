/**
 * @file FourthOrderNeville.hpp
 * @author Seppe Staelens
 * @brief Implementation of Neville's algorithm for fourth order interpolation of a point.
 * @version 0.1
 * @date 2025-01-27
 * @copyright Copyright (c) 2025
 */

#ifndef FOURTHORDERNEVILLE_HPP_
#define FOURTHORDERNEVILLE_HPP_

/**
 * @brief Class for the fourth-order accurate implementation of Neville's algorithm in 1 dimension.
 */
class FourthOrderNeville
{
private:
    //! Pointer to the x values of the data points.
    std::vector<double> *m_x = nullptr;
    //! Pointer to the y values of the data points.
    std::vector<double> *m_y = nullptr;
    //! The number of data points.
    int m_n;

public:
    /**
     * @brief Construct a new Fourth Order Neville object.
     */
    FourthOrderNeville() {}

    /**
     * @brief Set the data points for the interpolation.
     * @param x x values of the data points.
     * @param y y values of the data points.
     */
    void set_data_points(std::vector<double> *x, std::vector<double> *y)
    {
        if (x == nullptr || y == nullptr)
        {
            throw std::invalid_argument("Input vectors must not be null.");
        }
        if (x->size() != y->size())
        {
            throw std::invalid_argument("Input vectors must be of the same size.");
        }
        m_x = x;
        m_y = y;
        m_n = x->size();
    }

    /**
     * @brief Interpolates the value of a point using Neville's algorithm.
     * @param x0 The x value of the point to interpolate.
     * @return The interpolated value of the point.
     */
    double interpolate(double x0) const
    {
        int i_low = 0, i_high = m_n - 1;
        while (i_high - i_low > 1)
        {
            int i_mid = (i_low + i_high) / 2;
            if ((*m_x)[i_mid] < x0)
            {
                i_low = i_mid;
            }
            else
            {
                i_high = i_mid;
            }
        }

        if (i_low <= 1)
        {
            i_low = 2;
        }
        else if (i_high == m_n - 1)
        {
            i_low = m_n - 3;
        }

        int i0 = i_low - 2;
        int i1 = i_low - 1;
        int i2 = i_low;
        int i3 = i_low + 1;
        int i4 = i_low + 2;

        double x00 = (*m_x)[i0];
        double x11 = (*m_x)[i1];
        double x22 = (*m_x)[i2];
        double x33 = (*m_x)[i3];
        double x44 = (*m_x)[i4];

        double y00 = (*m_y)[i0];
        double y11 = (*m_y)[i1];
        double y22 = (*m_y)[i2];
        double y33 = (*m_y)[i3];
        double y44 = (*m_y)[i4];

        double y01 = ((x0 - x11) * y00 + (x00 - x0) * y11) / (x00 - x11);
        double y12 = ((x0 - x22) * y11 + (x11 - x0) * y22) / (x11 - x22);
        double y23 = ((x0 - x33) * y22 + (x22 - x0) * y33) / (x22 - x33);
        double y34 = ((x0 - x44) * y33 + (x33 - x0) * y44) / (x33 - x44);
        double y02 = ((x0 - x22) * y01 + (x00 - x0) * y12) / (x00 - x22);
        double y13 = ((x0 - x33) * y12 + (x11 - x0) * y23) / (x11 - x33);
        double y24 = ((x0 - x44) * y23 + (x22 - x0) * y34) / (x22 - x44);
        double y03 = ((x0 - x33) * y02 + (x00 - x0) * y13) / (x00 - x33);
        double y14 = ((x0 - x44) * y13 + (x11 - x0) * y24) / (x11 - x44);
        double y04 = ((x0 - x44) * y03 + (x00 - x0) * y14) / (x00 - x44);

        return y04;
    }
};

#endif /* FOURTHORDERNEVILLE_HPP_ */
