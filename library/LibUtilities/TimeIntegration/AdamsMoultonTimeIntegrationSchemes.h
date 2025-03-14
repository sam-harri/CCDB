///////////////////////////////////////////////////////////////////////////////
//
// File: AdamsMoultonTimeIntegrationSchemes.h
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
// Description: Combined header file for all Adams Moulton based time
// integration schemes.
//
///////////////////////////////////////////////////////////////////////////////

// Note : If adding a new integrator be sure to register the
// integrator with the Time Integration Scheme Facatory in
// SchemeInitializor.cpp.

#ifndef NEKTAR_LIB_UTILITIES_TIME_INTEGRATION_AM_TIME_INTEGRATION_SCHEME
#define NEKTAR_LIB_UTILITIES_TIME_INTEGRATION_AM_TIME_INTEGRATION_SCHEME

#define LUE LIB_UTILITIES_EXPORT

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/TimeIntegration/DIRKTimeIntegrationSchemes.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationSchemeGLM.h>

namespace Nektar::LibUtilities
{

///////////////////////////////////////////////////////////////////////////////
// Adams Moulton Order N

class AdamsMoultonTimeIntegrationScheme : public TimeIntegrationSchemeGLM
{
public:
    AdamsMoultonTimeIntegrationScheme(std::string variant, size_t order,
                                      std::vector<NekDouble> freeParams)
        : TimeIntegrationSchemeGLM(variant, order, freeParams)
    {
        // Currently up to 4th order is implemented.
        ASSERTL1(1 <= order && order <= 4,
                 "AdamsMoulton Time integration scheme bad order (1-4): " +
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
            AdamsMoultonTimeIntegrationScheme::SetupSchemeData(
                m_integration_phases[order - 2], order - 1);
        }

        // Last phase
        AdamsMoultonTimeIntegrationScheme::SetupSchemeData(
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
                DIRKTimeIntegrationScheme::SetupSchemeData(
                    m_integration_phases[0], 2);
                break;

            case 4:
                // Order 3
                DIRKTimeIntegrationScheme::SetupSchemeData(
                    m_integration_phases[0], 3);
                // Order 3
                DIRKTimeIntegrationScheme::SetupSchemeData(
                    m_integration_phases[1], 3);
                break;

            default:
                ASSERTL1(false,
                         "AdamsMoulton Time integration scheme bad order: " +
                             std::to_string(order));
        }
    }

    ~AdamsMoultonTimeIntegrationScheme() override
    {
    }

    static TimeIntegrationSchemeSharedPtr create(
        std::string variant, size_t order, std::vector<NekDouble> freeParams)
    {
        TimeIntegrationSchemeSharedPtr p =
            MemoryManager<AdamsMoultonTimeIntegrationScheme>::AllocateSharedPtr(
                variant, order, freeParams);

        return p;
    }

    static std::string className;

    LUE static void SetupSchemeData(TimeIntegrationAlgorithmGLMSharedPtr &phase,
                                    size_t order)
    {
        // clang-format off
        constexpr NekDouble coefficients[5][4] =
            { {      0.,       0.,      0.,     0. },
              // 1st Order
              {      1.,       0.,      0.,     0. },
              // 2nd Order
              {  1./ 2.,   1./ 2.,      0.,     0. },
              // 3rd Order
              {  5./12.,   8./12., -1./12.,     0. },
              // 4th Order
              {  9./24.,  19./24., -5./24.,  1./24. } };
        // clang-format on

        phase->m_schemeType = eDiagonallyImplicit;
        phase->m_order      = order;
        phase->m_name =
            std::string("AdamsMoultonOrder" + std::to_string(phase->m_order));

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

        // A/B Coefficient for first row first column
        phase->m_A[0][0][0] = coefficients[phase->m_order][0];
        phase->m_B[0][0][0] = coefficients[phase->m_order][0];

        // B evaluation value shuffling second row first column
        if (phase->m_order > 1)
        {
            phase->m_B[0][1][0] = 1.0;
        }

        // U/V Coefficient for first row first column
        phase->m_U[0][0] = 1.0;
        phase->m_V[0][0] = 1.0;

        // U/V Coefficients for first row additional columns
        for (size_t n = 1; n < phase->m_order; ++n)
        {
            phase->m_U[0][n] = coefficients[phase->m_order][n];
            phase->m_V[0][n] = coefficients[phase->m_order][n];
        }

        // V evaluation value shuffling row n column n-1
        for (size_t n = 2; n < phase->m_order; ++n)
        {
            phase->m_V[n][n - 1] = 1.0;
        }

        phase->m_numMultiStepValues         = 1;
        phase->m_numMultiStepImplicitDerivs = phase->m_order - 1;
        phase->m_numMultiStepExplicitDerivs = 0;
        phase->m_timeLevelOffset    = Array<OneD, size_t>(phase->m_numsteps);
        phase->m_timeLevelOffset[0] = 0;

        // For order > 1 derivatives are needed.
        for (size_t n = 1; n < phase->m_order; ++n)
        {
            phase->m_timeLevelOffset[n] = n - 1;
        }

        phase->CheckAndVerify();
    }

protected:
    LUE std::string v_GetName() const override
    {
        return std::string("AdamsMoulton");
    }

    LUE NekDouble v_GetTimeStability() const override
    {
        return 1.0;
    }

}; // end class AdamsMoultonTimeIntegrationScheme

////////////////////////////////////////////////////////////////////////////////
// Backwards compatibility
class AdamsMoultonOrder1TimeIntegrationScheme
    : public AdamsMoultonTimeIntegrationScheme
{
public:
    AdamsMoultonOrder1TimeIntegrationScheme(std::string variant, size_t order,
                                            std::vector<NekDouble> freeParams)
        : AdamsMoultonTimeIntegrationScheme("", 1, freeParams)
    {
        boost::ignore_unused(variant, order);
    }

    static TimeIntegrationSchemeSharedPtr create(
        [[maybe_unused]] std::string variant, [[maybe_unused]] size_t order,
        std::vector<NekDouble> freeParams)
    {
        TimeIntegrationSchemeSharedPtr p =
            MemoryManager<AdamsMoultonTimeIntegrationScheme>::AllocateSharedPtr(
                "", 1, freeParams);

        return p;
    }

    static std::string className;

protected:
    static std::string TimeIntegrationMethodLookupId;

}; // end class AdamsMoultonOrder1TimeIntegrationScheme

class AdamsMoultonOrder2TimeIntegrationScheme
    : public AdamsMoultonTimeIntegrationScheme
{
public:
    AdamsMoultonOrder2TimeIntegrationScheme(std::string variant, size_t order,
                                            std::vector<NekDouble> freeParams)
        : AdamsMoultonTimeIntegrationScheme("", 2, freeParams)
    {
        boost::ignore_unused(variant, order);
    }

    static TimeIntegrationSchemeSharedPtr create(
        [[maybe_unused]] std::string variant, [[maybe_unused]] size_t order,
        std::vector<NekDouble> freeParams)
    {
        TimeIntegrationSchemeSharedPtr p =
            MemoryManager<AdamsMoultonTimeIntegrationScheme>::AllocateSharedPtr(
                "", 2, freeParams);

        return p;
    }

    static std::string className;

protected:
    static std::string TimeIntegrationMethodLookupId;

}; // end class AdamsMoultonOrder2TimeIntegrationScheme

class AdamsMoultonOrder3TimeIntegrationScheme
    : public AdamsMoultonTimeIntegrationScheme
{
public:
    AdamsMoultonOrder3TimeIntegrationScheme(std::string variant, size_t order,
                                            std::vector<NekDouble> freeParams)
        : AdamsMoultonTimeIntegrationScheme("", 3, freeParams)
    {
        boost::ignore_unused(variant, order);
    }

    static TimeIntegrationSchemeSharedPtr create(
        [[maybe_unused]] std::string variant, [[maybe_unused]] size_t order,
        std::vector<NekDouble> freeParams)
    {
        TimeIntegrationSchemeSharedPtr p =
            MemoryManager<AdamsMoultonTimeIntegrationScheme>::AllocateSharedPtr(
                "", 3, freeParams);

        return p;
    }

    static std::string className;

protected:
    static std::string TimeIntegrationMethodLookupId;

}; // end class AdamsMoultonOrder3TimeIntegrationScheme

class AdamsMoultonOrder4TimeIntegrationScheme
    : public AdamsMoultonTimeIntegrationScheme
{
public:
    AdamsMoultonOrder4TimeIntegrationScheme(std::string variant, size_t order,
                                            std::vector<NekDouble> freeParams)
        : AdamsMoultonTimeIntegrationScheme("", 4, freeParams)
    {
        boost::ignore_unused(variant, order);
    }

    static TimeIntegrationSchemeSharedPtr create(
        [[maybe_unused]] std::string variant, [[maybe_unused]] size_t order,
        std::vector<NekDouble> freeParams)
    {
        TimeIntegrationSchemeSharedPtr p =
            MemoryManager<AdamsMoultonTimeIntegrationScheme>::AllocateSharedPtr(
                "", 4, freeParams);

        return p;
    }

    static std::string className;

protected:
    static std::string TimeIntegrationMethodLookupId;

}; // end class AdamsMoultonOrder4TimeIntegrationScheme

} // namespace Nektar::LibUtilities

#endif
