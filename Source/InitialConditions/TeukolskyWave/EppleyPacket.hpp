
#ifndef EPPLEYPACKET_HPP_
#define EPPLEYPACKET_HPP_

#include "cmath"
#include "EppleyPacketParams.hpp"
#include <vector>

struct EppleyPacketDerivs
{
    double F0, F1, F2, F3, F4;
};

//! Base EppleyPacket class
class EppleyPacket
{
  private:          // private member variables/arrays
    double m_amplitude; //!< amplitude of the packet
    double m_sigma;      //!< width of the packet
    double m_radial_offset; //!< offset for radial coordinate to not center the wave on the center

  public:
    //! F function and its first four derivatives, where x = r \pm t
    EppleyPacketDerivs get_F_derivs(double x) const;

    EppleyPacket(EppleyPacket_params_t a_params):
        m_amplitude(a_params.amplitude),
        m_sigma(a_params.sigma),
        m_radial_offset(a_params.radial_offset)
    {};
};

struct EvenEppleyPacketCoefficients
{
    double A, B, C;
};
//! Base class for the even parity Eppley packets
class EvenEppleyPacket : EppleyPacket
{
  public:
    EvenEppleyPacketCoefficients get_ABC(double r) const;

    EvenEppleyPacket(EppleyPacket_params_t m_params) : EppleyPacket(m_params) {};
};

struct OddEppleyPacketCoefficients
{
    double K, L;
};
//! Base class for the odd parity Eppley packets
class OddEppleyPacket : EppleyPacket
{
  public:
    OddEppleyPacketCoefficients get_KL(double r) const;

    OddEppleyPacket(EppleyPacket_params_t m_params) : EppleyPacket(m_params) {};
};

struct EppleyPacketMetricComponents
{
    double gxx, gxy, gxz, gyy, gyz, gzz;
};
//! Specific Eppley Packet classes
//! m = 0 and m = 2 are the only ones implemented so far for the even parity case
class EppleyPacketM0 : public EvenEppleyPacket
{
  public:
    EppleyPacketM0(EppleyPacket_params_t m_params) : EvenEppleyPacket(m_params) {}

    EppleyPacketMetricComponents get_metric_components(double x, double y, double z) const;
};

class EppleyPacketM2 : public EvenEppleyPacket
{
  public:
    EppleyPacketM2(EppleyPacket_params_t m_params) : EvenEppleyPacket(m_params) {}

    EppleyPacketMetricComponents get_metric_components(double x, double y, double z) const;
};

//! m = 2 class for the odd parity case
class OddEppleyPacketM2 : public OddEppleyPacket
{
  public:
    OddEppleyPacketM2(EppleyPacket_params_t m_params) : OddEppleyPacket(m_params) {}

    EppleyPacketMetricComponents get_metric_components(double x, double y, double z) const;
};

#include "EppleyPacket.impl.hpp"

#endif /* EPPLEYPACKET_HPP_ */
