///////////////////////////////////////////////////////////////////////////////
//
// File: BDFImplicitTimeIntegrationSchemes.h
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
// Description: Combined header file for all BDF Implicit based time integration
// schemes.
//
///////////////////////////////////////////////////////////////////////////////

// Note : If adding a new integrator be sure to register the
// integrator with the Time Integration Scheme Facatory in
// SchemeInitializor.cpp.

#ifndef NEKTAR_LIB_UTILITIES_TIME_INTEGRATION_BDF_TIME_INTEGRATION_SCHEME
#define NEKTAR_LIB_UTILITIES_TIME_INTEGRATION_BDF_TIME_INTEGRATION_SCHEME

#define LUE LIB_UTILITIES_EXPORT

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/TimeIntegration/DIRKTimeIntegrationSchemes.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationSchemeGLM.h>

namespace Nektar::LibUtilities
{

///////////////////////////////////////////////////////////////////////////////
// BDF Implicit Order N

class BDFImplicitTimeIntegrationScheme : public TimeIntegrationSchemeGLM
{
public:
    BDFImplicitTimeIntegrationScheme(std::string variant, size_t order,
                                     std::vector<NekDouble> freeParams)
        : TimeIntegrationSchemeGLM(variant, order, freeParams)
    {
        // Currently up to 4th order is implemented.
        ASSERTL1(1 <= order && order <= 4,
                 "BDFImplicit Time integration scheme bad order (1-4): " +
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
            BDFImplicitTimeIntegrationScheme::SetupSchemeData(
                m_integration_phases[order - 2], order - 1);
        }

        // Last phase
        BDFImplicitTimeIntegrationScheme::SetupSchemeData(
            m_integration_phases[order - 1], order);

        // Initial phases
        switch (order)
        {
            case 1:
                // No intial phase.
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
                         "BDFImplicit Time integration scheme bad order: " +
                             std::to_string(order));
        }
    }

    ~BDFImplicitTimeIntegrationScheme() override
    {
    }

    static TimeIntegrationSchemeSharedPtr create(
        std::string variant, size_t order, std::vector<NekDouble> freeParams)
    {
        TimeIntegrationSchemeSharedPtr p =
            MemoryManager<BDFImplicitTimeIntegrationScheme>::AllocateSharedPtr(
                variant, order, freeParams);

        return p;
    }

    static std::string className;

    LUE static void SetupSchemeData(TimeIntegrationAlgorithmGLMSharedPtr &phase,
                                    size_t order)
    {
        constexpr NekDouble ABcoefficients[5] = {0.,
                                                 1.,         // 1st Order
                                                 2. / 3.,    // 2nd Order
                                                 6. / 11.,   // 3rd Order
                                                 12. / 25.}; // 4th Order

        // clang-format off
        constexpr NekDouble UVcoefficients[5][4] =
            { {       0.,      0.,     0.,       0. },
              // 1st Order
              {       1.,      0.,     0.,       0. },
              // 2nd Order
              {  4./ 3.,  -1./ 3.,     0.,       0. },
              // 3rd Order
              { 18./11.,  -9./11.,  2./11.,      0. },
              // 4th Order
              { 48./25., -36./25., 16./25., -3./25. } };
        // clang-format on

        phase->m_schemeType = eDiagonallyImplicit;
        phase->m_order      = order;
        phase->m_name =
            std::string("BDFImplicitOrder" + std::to_string(phase->m_order));

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
        phase->m_A[0][0][0] = ABcoefficients[phase->m_order];
        phase->m_B[0][0][0] = ABcoefficients[phase->m_order];

        // U/V Coefficients for first row additional columns
        for (size_t n = 0; n < phase->m_order; ++n)
        {
            phase->m_U[0][n] = UVcoefficients[phase->m_order][n];
            phase->m_V[0][n] = UVcoefficients[phase->m_order][n];
        }

        // V evaluation value shuffling row n column n-1
        for (size_t n = 1; n < phase->m_order; ++n)
        {
            phase->m_V[n][n - 1] = 1.0;
        }

        phase->m_numMultiStepValues         = phase->m_order;
        phase->m_numMultiStepImplicitDerivs = 0;
        phase->m_numMultiStepExplicitDerivs = 0;
        phase->m_timeLevelOffset = Array<OneD, size_t>(phase->m_numsteps);
        for (size_t n = 0; n < phase->m_order; ++n)
        {
            phase->m_timeLevelOffset[n] = n;
        }

        phase->CheckAndVerify();
    }

protected:
    LUE std::string v_GetName() const override
    {
        return std::string("BDFImplicit");
    }

    LUE NekDouble v_GetTimeStability() const override
    {
        return 1.0;
    }

}; // end class BDFImplicitTimeIntegrator

////////////////////////////////////////////////////////////////////////////////
// Backwards compatibility
class BDFImplicitOrder1TimeIntegrationScheme
    : public BDFImplicitTimeIntegrationScheme
{
public:
    BDFImplicitOrder1TimeIntegrationScheme(std::string variant, size_t order,
                                           std::vector<NekDouble> freeParams)
        : BDFImplicitTimeIntegrationScheme("", 1, freeParams)
    {
        boost::ignore_unused(variant, order);
    }

    static TimeIntegrationSchemeSharedPtr create(
        [[maybe_unused]] std::string variant, [[maybe_unused]] size_t order,
        std::vector<NekDouble> freeParams)
    {
        TimeIntegrationSchemeSharedPtr p =
            MemoryManager<BDFImplicitTimeIntegrationScheme>::AllocateSharedPtr(
                "", 1, freeParams);
        return p;
    }

    static std::string className;

protected:
    static std::string TimeIntegrationMethodLookupId;

}; // end class BDFImplicitOrder1TimeIntegrationScheme

class BDFImplicitOrder2TimeIntegrationScheme
    : public BDFImplicitTimeIntegrationScheme
{
public:
    BDFImplicitOrder2TimeIntegrationScheme(std::string variant, size_t order,
                                           std::vector<NekDouble> freeParams)
        : BDFImplicitTimeIntegrationScheme("", 2, freeParams)
    {
        boost::ignore_unused(variant, order);
    }

    static TimeIntegrationSchemeSharedPtr create(
        [[maybe_unused]] std::string variant, [[maybe_unused]] size_t order,
        std::vector<NekDouble> freeParams)
    {
        TimeIntegrationSchemeSharedPtr p =
            MemoryManager<BDFImplicitTimeIntegrationScheme>::AllocateSharedPtr(
                "", 2, freeParams);
        return p;
    }

    static std::string className;

protected:
    static std::string TimeIntegrationMethodLookupId;

}; // end class BDFImplicitOrder2TimeIntegrationScheme

class BDFImplicitOrder3TimeIntegrationScheme
    : public BDFImplicitTimeIntegrationScheme
{
public:
    BDFImplicitOrder3TimeIntegrationScheme(std::string variant, size_t order,
                                           std::vector<NekDouble> freeParams)
        : BDFImplicitTimeIntegrationScheme("", 3, freeParams)
    {
        boost::ignore_unused(variant, order);
    }

    static TimeIntegrationSchemeSharedPtr create(
        [[maybe_unused]] std::string variant, [[maybe_unused]] size_t order,
        std::vector<NekDouble> freeParams)
    {
        TimeIntegrationSchemeSharedPtr p =
            MemoryManager<BDFImplicitTimeIntegrationScheme>::AllocateSharedPtr(
                "", 3, freeParams);
        return p;
    }

    static std::string className;

protected:
    static std::string TimeIntegrationMethodLookupId;

}; // end class BDFImplicitOrder3TimeIntegrationScheme

class BDFImplicitOrder4TimeIntegrationScheme
    : public BDFImplicitTimeIntegrationScheme
{
public:
    BDFImplicitOrder4TimeIntegrationScheme(std::string variant, size_t order,
                                           std::vector<NekDouble> freeParams)
        : BDFImplicitTimeIntegrationScheme("", 4, freeParams)
    {
        boost::ignore_unused(variant, order);
    }

    static TimeIntegrationSchemeSharedPtr create(
        [[maybe_unused]] std::string variant, [[maybe_unused]] size_t order,
        std::vector<NekDouble> freeParams)
    {
        TimeIntegrationSchemeSharedPtr p =
            MemoryManager<BDFImplicitTimeIntegrationScheme>::AllocateSharedPtr(
                "", 4, freeParams);
        return p;
    }

    static std::string className;

protected:
    static std::string TimeIntegrationMethodLookupId;

}; // end class BDFImplicitOrder4TimeIntegrationScheme

} // namespace Nektar::LibUtilities

#endif
