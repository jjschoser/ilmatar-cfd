#include "EquationOfState.H"

#include <cassert>
#ifdef DEBUG
    #include <iostream>
#endif

IdealGas::IdealGas(const REAL gamma) : m_gamma(gamma)
{
    assert(gamma > 1.0);
}

REAL IdealGas::getGamma() const
{
    return m_gamma;
}

REAL IdealGas::getPressure(const REAL rho, const REAL e) const
{
    #ifdef DEBUG
        if(rho <= 0.0 || e <= 0.0)
        {
            std::cout << "Invalid EoS input rho = " << rho << ", e = " << e << std::endl;
        }
    #endif
    return (m_gamma - 1.0) * rho * e;
}

REAL IdealGas::getSpecificInternalEnergy(const REAL rho, const REAL p) const
{
    #ifdef DEBUG
        if(rho <= 0.0 || p <= 0.0)
        {
            std::cout << "Invalid EoS input rho = " << rho << ", p = " << p << std::endl;
        }
    #endif
    return p / ((m_gamma - 1.0) * rho);
}

REAL IdealGas::getSoundSpeed(const REAL rho, const REAL p) const
{
    #ifdef DEBUG
        if(rho <= 0.0 || p <= 0.0)
        {
            std::cout << "Invalid EoS input rho = " << rho << ", p = " << p << std::endl;
        }
    #endif
    return std::sqrt(m_gamma * p / rho);
}