///////////////////////////////////////////////////////////////////////////////
//
// File: AdamsBashforthTimeIntegrationSchemes.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2018 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Combined header file for all Adams Bashforth based time
// integration schemes.
//
///////////////////////////////////////////////////////////////////////////////

// Note : If adding a new integrator be sure to register the
// integrator with the Time Integration Scheme Facatory in
// SchemeInitializor.cpp.

#ifndef NEKTAR_LIB_UTILITIES_TIME_INTEGRATION_AB_TIME_INTEGRATION_SCHEME
#define NEKTAR_LIB_UTILITIES_TIME_INTEGRATION_AB_TIME_INTEGRATION_SCHEME

#define LUE LIB_UTILITIES_EXPORT

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/TimeIntegration/RungeKuttaTimeIntegrationSchemes.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationSchemeGLM.h>

namespace Nektar::LibUtilities
{

///////////////////////////////////////////////////////////////////////////////
// Adams Bashforth Order N

class AdamsBashforthTimeIntegrationScheme : public TimeIntegrationSchemeGLM
{
public:
    AdamsBashforthTimeIntegrationScheme(std::string variant, size_t order,
                                        std::vector<NekDouble> freeParams)
        : TimeIntegrationSchemeGLM(variant, order, freeParams)
    {
        // Currently up to 4th order is implemented.
        ASSERTL1(1 <= order && order <= 4,
                 "AdamsBashforth Time integration scheme bad order (1-4): " +
                     std::to_string(order));

        m_integration_phases = TimeIntegrationAlgorithmGLMVector(order);

        for (size_t n = 0; n < order; ++n)
        {
            m_integration_phases[n] = TimeIntegrationAlgorithmGLMSharedPtr(
                new TimeIntegrationAlgorithmGLM(this));
        }

        // Next to last phase
        if (order > 1)
        {
            AdamsBashforthTimeIntegrationScheme::SetupSchemeData(
                m_integration_phases[order - 2], order - 1);
        }

        // Last phase
        AdamsBashforthTimeIntegrationScheme::SetupSchemeData(
            m_integration_phases[order - 1], order);

        // Initial phases
        switch (order)
        {
            case 1:
                // No intial phases.
                break;

            case 2:
                // Intial phase set above
                break;

            case 3:
                // Order 2
                RungeKuttaTimeIntegrationScheme::SetupSchemeData(
                    m_integration_phases[0], "", 2, std::vector<NekDouble>());
                break;

            case 4:
                // SSP Order 3
                RungeKuttaTimeIntegrationScheme::SetupSchemeData(
                    m_integration_phases[0], "SSP", 3,
                    std::vector<NekDouble>());
                // SSP Order 3
                RungeKuttaTimeIntegrationScheme::SetupSchemeData(
                    m_integration_phases[1], "SSP", 3,
                    std::vector<NekDouble>());
                break;

            default:
                ASSERTL1(false,
                         "AdamsBashforth Time integration scheme bad order: " +
                             std::to_string(order));
        }
    }

    ~AdamsBashforthTimeIntegrationScheme() override
    {
    }

    static TimeIntegrationSchemeSharedPtr create(
        std::string variant, size_t order, std::vector<NekDouble> freeParams)
    {
        TimeIntegrationSchemeSharedPtr p = MemoryManager<
            AdamsBashforthTimeIntegrationScheme>::AllocateSharedPtr(variant,
                                                                    order,
                                                                    freeParams);

        return p;
    }

    static std::string className;

    LUE static void SetupSchemeData(TimeIntegrationAlgorithmGLMSharedPtr &phase,
                                    size_t order)
    {
        // clang-format off
        constexpr NekDouble coefficients[5][4] =
            { {      0.,       0.,      0.,      0. },
              // 1st Order
              {      1.,       0.,      0.,      0. },
              // 2nd Order
              {  3./ 2.,  -1./ 2.,      0.,      0. },
              // 3rd Order
              { 23./12., -16./12.,  5./12.,      0. },
              // 4th Order
              { 55./24., -59./24., 37./24., -9./24.} };
        // clang-format on

        phase->m_schemeType = eExplicit;
        phase->m_order      = order;
        phase->m_name =
            std::string("AdamsBashforthOrder" + std::to_string(phase->m_order));

        phase->m_numsteps  = phase->m_order;
        phase->m_numstages = 1;

        phase->m_A = Array<OneD, Array<TwoD, NekDouble>>(1);
        phase->m_B = Array<OneD, Array<TwoD, NekDouble>>(1);

        phase->m_A[0] =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numstages, 0.0);
        phase->m_B[0] =
            Array<TwoD, NekDouble>(phase->m_numsteps, phase->m_numstages, 0.0);
        phase->m_U =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numsteps, 0.0);
        phase->m_V =
            Array<TwoD, NekDouble>(phase->m_numsteps, phase->m_numsteps, 0.0);

        // Coefficients

        // B Coefficient for first row first column
        phase->m_B[0][0][0] = coefficients[phase->m_order][0];

        // B evaluation value shuffling second row first column
        if (phase->m_order > 1)
        {
            phase->m_B[0][1][0] = 1.0;
        }

        // U Curent time step evaluation first row first column
        phase->m_U[0][0] = 1.0;
        phase->m_V[0][0] = 1.0;

        // V Coefficients for first row additional columns
        for (size_t n = 1; n < phase->m_order; ++n)
        {
            phase->m_V[0][n] = coefficients[phase->m_order][n];
        }

        // V evaluation value shuffling row n column n-1
        for (size_t n = 2; n < phase->m_order; ++n)
        {
            phase->m_V[n][n - 1] = 1.0;
        }

        phase->m_numMultiStepValues         = 1;
        phase->m_numMultiStepImplicitDerivs = 0;
        phase->m_numMultiStepExplicitDerivs = phase->m_order - 1;
        phase->m_timeLevelOffset    = Array<OneD, size_t>(phase->m_numsteps);
        phase->m_timeLevelOffset[0] = 0;

        // For order > 1 derivatives are needed.
        for (size_t n = 1; n < phase->m_order; ++n)
        {
            phase->m_timeLevelOffset[n] = n;
        }

        phase->CheckAndVerify();
    }

protected:
    LUE std::string v_GetName() const override
    {
        return std::string("AdamsBashforth");
    }

    LUE NekDouble v_GetTimeStability() const override
    {
        if (GetOrder() == 1)
        {
            return 2.0;
        }
        else if (GetOrder() == 2)
        {
            return 1.0;
        }
        else if (GetOrder() == 3)
        {
            return 0.545454545454545;
        }
        else if (GetOrder() == 4)
        {
            return 0.3;
        }
        else
        {
            return 2.0;
        }
    }

}; // end class AdamsBashforthTimeIntegrationScheme

////////////////////////////////////////////////////////////////////////////////
// Backwards compatibility
class AdamsBashforthOrder1TimeIntegrationScheme
    : public AdamsBashforthTimeIntegrationScheme
{
public:
    AdamsBashforthOrder1TimeIntegrationScheme(std::string variant, size_t order,
                                              std::vector<NekDouble> freeParams)
        : AdamsBashforthTimeIntegrationScheme("", 1, freeParams)
    {
        boost::ignore_unused(variant, order);
    }

    static TimeIntegrationSchemeSharedPtr create(
        [[maybe_unused]] std::string variant, [[maybe_unused]] size_t order,
        std::vector<NekDouble> freeParams)
    {
        TimeIntegrationSchemeSharedPtr p = MemoryManager<
            AdamsBashforthTimeIntegrationScheme>::AllocateSharedPtr("", 1,
                                                                    freeParams);
        return p;
    }

    static std::string className;

protected:
    static std::string TimeIntegrationMethodLookupId;

}; // end class AdamsBashforthOrder1TimeIntegrationScheme

class AdamsBashforthOrder2TimeIntegrationScheme
    : public AdamsBashforthTimeIntegrationScheme
{
public:
    AdamsBashforthOrder2TimeIntegrationScheme(std::string variant, size_t order,
                                              std::vector<NekDouble> freeParams)
        : AdamsBashforthTimeIntegrationScheme("", 2, freeParams)
    {
        boost::ignore_unused(variant, order);
    }

    static TimeIntegrationSchemeSharedPtr create(
        [[maybe_unused]] std::string variant, [[maybe_unused]] size_t order,
        std::vector<NekDouble> freeParams)
    {
        TimeIntegrationSchemeSharedPtr p = MemoryManager<
            AdamsBashforthTimeIntegrationScheme>::AllocateSharedPtr("", 2,
                                                                    freeParams);
        return p;
    }

    static std::string className;

protected:
    static std::string TimeIntegrationMethodLookupId;

}; // end class AdamsBashforthOrder2TimeIntegrationScheme

class AdamsBashforthOrder3TimeIntegrationScheme
    : public AdamsBashforthTimeIntegrationScheme
{
public:
    AdamsBashforthOrder3TimeIntegrationScheme(std::string variant, size_t order,
                                              std::vector<NekDouble> freeParams)
        : AdamsBashforthTimeIntegrationScheme("", 3, freeParams)
    {
        boost::ignore_unused(variant, order);
    }

    static TimeIntegrationSchemeSharedPtr create(
        [[maybe_unused]] std::string variant, [[maybe_unused]] size_t order,
        std::vector<NekDouble> freeParams)
    {
        TimeIntegrationSchemeSharedPtr p = MemoryManager<
            AdamsBashforthTimeIntegrationScheme>::AllocateSharedPtr("", 3,
                                                                    freeParams);
        return p;
    }

    static std::string className;

protected:
    static std::string TimeIntegrationMethodLookupId;

}; // end class AdamsBashforthOrder3TimeIntegrationScheme

class AdamsBashforthOrder4TimeIntegrationScheme
    : public AdamsBashforthTimeIntegrationScheme
{
public:
    AdamsBashforthOrder4TimeIntegrationScheme(std::string variant, size_t order,
                                              std::vector<NekDouble> freeParams)
        : AdamsBashforthTimeIntegrationScheme("", 4, freeParams)
    {
        boost::ignore_unused(variant, order);
    }

    static TimeIntegrationSchemeSharedPtr create(
        [[maybe_unused]] std::string variant, [[maybe_unused]] size_t order,
        std::vector<NekDouble> freeParams)
    {
        TimeIntegrationSchemeSharedPtr p = MemoryManager<
            AdamsBashforthTimeIntegrationScheme>::AllocateSharedPtr("", 4,
                                                                    freeParams);
        return p;
    }

    static std::string className;

protected:
    static std::string TimeIntegrationMethodLookupId;

}; // end class AdamsBashforthOrder4TimeIntegrationScheme

} // namespace Nektar::LibUtilities

#endif
