#include "EquationOfState.H"

#ifdef DEBUG
    #include <cassert>
#endif

IdealGas::IdealGas(const REAL gamma) : m_gamma(gamma)
{
    #ifdef DEBUG
        assert(gamma > 1.0);
    #endif
}

REAL IdealGas::getGamma() const
{
    return m_gamma;
}

REAL IdealGas::getPressure(const REAL rho, const REAL e) const
{
    #ifdef DEBUG
        assert(rho > 0.0 && e > 0.0);
    #endif
    return (m_gamma - 1.0) * rho * e;
}

REAL IdealGas::getSpecificInternalEnergy(const REAL rho, const REAL p) const
{
    #ifdef DEBUG
        assert(rho > 0.0 && p > 0.0);
    #endif
    return p / ((m_gamma - 1.0) * rho);
}

REAL IdealGas::getSoundSpeed(const REAL rho, const REAL p) const
{
    #ifdef DEBUG
        assert(rho > 0.0 && p > 0.0);
    #endif
    return std::sqrt(m_gamma * p / rho);
}