///////////////////////////////////////////////////////////////////////////////
//
// File: IMEXTimeIntegrationSchemes.h
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
// Description: Combined header file for all basic IMEX time integration
// schemes.
//
///////////////////////////////////////////////////////////////////////////////

// Note : If adding a new integrator be sure to register the
// integrator with the Time Integration Scheme Facatory in
// SchemeInitializor.cpp.

#ifndef NEKTAR_LIB_UTILITIES_TIME_INTEGRATION_IMEX_TIME_INTEGRATION_SCHEME
#define NEKTAR_LIB_UTILITIES_TIME_INTEGRATION_IMEX_TIME_INTEGRATION_SCHEME

#define LUE LIB_UTILITIES_EXPORT

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/TimeIntegration/IMEXGearTimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/IMEXdirkTimeIntegrationSchemes.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationSchemeGLM.h>

namespace Nektar::LibUtilities
{

///////////////////////////////////////////////////////////////////////////////
// IMEX Order N

class IMEXTimeIntegrationScheme : public TimeIntegrationSchemeGLM
{
public:
    IMEXTimeIntegrationScheme(std::string variant, size_t order,
                              std::vector<NekDouble> freeParams)
        : TimeIntegrationSchemeGLM(variant, order, freeParams)
    {
        if (variant == "dirk" || variant == "DIRK")
        {
            ASSERTL1(freeParams.size() == 2,
                     "IMEX DIRK Time integration scheme invalid number "
                     "of free parameters, expected two "
                     "<implicit stages, explicit stages>, received " +
                         std::to_string(freeParams.size()));

            size_t s     = freeParams[0];
            size_t sigma = freeParams[1];

            m_integration_phases    = TimeIntegrationAlgorithmGLMVector(1);
            m_integration_phases[0] = TimeIntegrationAlgorithmGLMSharedPtr(
                new TimeIntegrationAlgorithmGLM(this));

            if (order == 1 && s == 1 && sigma == 1)
            {
                // This phase is Forward Backward Euler which has two steps.
                IMEXdirkTimeIntegrationScheme::SetupSchemeData_1_1_1(
                    m_integration_phases[0]);
            }
            else
            {
                IMEXdirkTimeIntegrationScheme::SetupSchemeData(
                    m_integration_phases[0], order, freeParams);
            }
        }
        else if (variant == "Gear")
        {
            m_integration_phases    = TimeIntegrationAlgorithmGLMVector(2);
            m_integration_phases[0] = TimeIntegrationAlgorithmGLMSharedPtr(
                new TimeIntegrationAlgorithmGLM(this));
            m_integration_phases[1] = TimeIntegrationAlgorithmGLMSharedPtr(
                new TimeIntegrationAlgorithmGLM(this));

            IMEXdirkTimeIntegrationScheme::SetupSchemeData(
                m_integration_phases[0], 2, std::vector<NekDouble>{2, 2});
            IMEXGearTimeIntegrationScheme::SetupSchemeData(
                m_integration_phases[1]);
        }
        else if (variant == "")
        {
            // Currently up to 4th order is implemented.
            ASSERTL1(1 <= order && order <= 4,
                     "IMEX Time integration scheme bad order (1-4): " +
                         std::to_string(order));

            m_integration_phases = TimeIntegrationAlgorithmGLMVector(order);

            for (size_t n = 0; n < order; ++n)
            {
                m_integration_phases[n] = TimeIntegrationAlgorithmGLMSharedPtr(
                    new TimeIntegrationAlgorithmGLM(this));
            }

            // Last phase
            IMEXTimeIntegrationScheme::SetupSchemeData(
                m_integration_phases[order - 1], order);

            // Initial phases
            switch (order)
            {
                case 1:
                    // No intial phase.
                    break;

                case 2:
                    IMEXTimeIntegrationScheme::SetupSchemeData(
                        m_integration_phases[0], 1);
                    break;

                case 3:
                    IMEXdirkTimeIntegrationScheme::SetupSchemeData(
                        m_integration_phases[0], 3,
                        std::vector<NekDouble>{3, 4});
                    IMEXdirkTimeIntegrationScheme::SetupSchemeData(
                        m_integration_phases[1], 3,
                        std::vector<NekDouble>{3, 4});
                    break;

                case 4:
                    IMEXdirkTimeIntegrationScheme::SetupSchemeData(
                        m_integration_phases[0], 3,
                        std::vector<NekDouble>{3, 4});
                    IMEXdirkTimeIntegrationScheme::SetupSchemeData(
                        m_integration_phases[1], 3,
                        std::vector<NekDouble>{3, 4});
                    IMEXdirkTimeIntegrationScheme::SetupSchemeData(
                        m_integration_phases[2], 3,
                        std::vector<NekDouble>{3, 4});
                    break;

                default:
                    ASSERTL1(false, "IMEX Time integration scheme bad order: " +
                                        std::to_string(order));
            }
        }
        else
        {
            ASSERTL1(false, "IMEX Time integration scheme bad variant: " +
                                variant + ". Must be blank, 'dirk' or 'Gear'");
        }
    }

    ~IMEXTimeIntegrationScheme() override
    {
    }

    static TimeIntegrationSchemeSharedPtr create(
        std::string variant, size_t order, std::vector<NekDouble> freeParams)
    {
        TimeIntegrationSchemeSharedPtr p =
            MemoryManager<IMEXTimeIntegrationScheme>::AllocateSharedPtr(
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

        // Nsteps = 2 * order

        // clang-format off
        constexpr NekDouble UVcoefficients[5][8] =
            { {         0.,    0.,     0.,        0.,
                        0.,    0.,     0.,        0. },
              // 1st Order
              {         1.,    1.,     0.,        0.,
                        0.,    0.,     0.,        0. },
              // 2nd Order
              {  4./ 3.,  -1./ 3.,  4./3.,   -2./ 3.,
                 0.,           0.,     0.,        0. },
              // 3rd Order
              {  18./11.,  -9./11.,  2./11.,  18./11.,
                -18./11.,   6./11.,      0.,       0. },
              // 4th Order
              { 48./25., -36./25., 16./25.,  -3./25.,
                48./25., -72./25., 48./25., -12./25. } };
        // clang-format on

        phase->m_schemeType = eIMEX;
        phase->m_order      = order;
        phase->m_name =
            std::string("IMEXOrder" + std::to_string(phase->m_order));

        phase->m_numsteps  = 2 * phase->m_order;
        phase->m_numstages = 1;

        phase->m_A = Array<OneD, Array<TwoD, NekDouble>>(2);
        phase->m_B = Array<OneD, Array<TwoD, NekDouble>>(2);

        phase->m_A[0] =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numstages, 0.0);
        phase->m_B[0] =
            Array<TwoD, NekDouble>(phase->m_numsteps, phase->m_numstages, 0.0);
        phase->m_A[1] =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numstages, 0.0);
        phase->m_B[1] =
            Array<TwoD, NekDouble>(phase->m_numsteps, phase->m_numstages, 0.0);
        phase->m_U =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numsteps, 0.0);
        phase->m_V =
            Array<TwoD, NekDouble>(phase->m_numsteps, phase->m_numsteps, 0.0);

        // Coefficients
        phase->m_B[1][phase->m_order][0] = 1.0;

        phase->m_A[0][0][0] = ABcoefficients[phase->m_order];
        phase->m_B[0][0][0] = ABcoefficients[phase->m_order];

        for (size_t n = 0; n < 2 * phase->m_order; ++n)
        {
            phase->m_U[0][n] = UVcoefficients[phase->m_order][n];
            phase->m_V[0][n] = UVcoefficients[phase->m_order][n];
        }

        // V evaluation value shuffling row n column n-1
        for (size_t n = 1; n < 2 * phase->m_order; ++n)
        {
            if (n != phase->m_order)
            {
                phase->m_V[n][n - 1] = 1.0; // constant 1
            }
        }

        phase->m_numMultiStepValues         = phase->m_order;
        phase->m_numMultiStepImplicitDerivs = 0;
        phase->m_numMultiStepExplicitDerivs = phase->m_order;

        phase->m_timeLevelOffset = Array<OneD, size_t>(phase->m_numsteps);

        // Values and derivatives needed.
        for (size_t n = 0; n < phase->m_order; ++n)
        {
            phase->m_timeLevelOffset[n]                  = n;
            phase->m_timeLevelOffset[phase->m_order + n] = n;
        }

        phase->CheckAndVerify();
    }

protected:
    LUE std::string v_GetFullName() const override
    {
        return m_integration_phases.back()->m_name;
    }

    LUE std::string v_GetName() const override
    {
        return std::string("IMEX");
    }

    LUE NekDouble v_GetTimeStability() const override
    {
        return 1.0;
    }

}; // end class IMEXTimeIntegrationScheme

////////////////////////////////////////////////////////////////////////////////
// Backwards compatibility
class IMEXOrder1TimeIntegrationScheme : public IMEXTimeIntegrationScheme
{
public:
    IMEXOrder1TimeIntegrationScheme(std::string variant, size_t order,
                                    std::vector<NekDouble> freeParams)
        : IMEXTimeIntegrationScheme("", 1, freeParams)
    {
        boost::ignore_unused(variant, order);
    }

    static TimeIntegrationSchemeSharedPtr create(
        [[maybe_unused]] std::string variant, [[maybe_unused]] size_t order,
        std::vector<NekDouble> freeParams)
    {
        TimeIntegrationSchemeSharedPtr p =
            MemoryManager<IMEXTimeIntegrationScheme>::AllocateSharedPtr(
                "", 1, freeParams);
        return p;
    }

    static std::string className;

protected:
    static std::string TimeIntegrationMethodLookupId;

}; // end class IMEXOrder1TimeIntegrationScheme

class IMEXOrder2TimeIntegrationScheme : public IMEXTimeIntegrationScheme
{
public:
    IMEXOrder2TimeIntegrationScheme(std::string variant, size_t order,
                                    std::vector<NekDouble> freeParams)
        : IMEXTimeIntegrationScheme("", 2, freeParams)
    {
        boost::ignore_unused(variant, order);
    }

    static TimeIntegrationSchemeSharedPtr create(
        [[maybe_unused]] std::string variant, [[maybe_unused]] size_t order,
        std::vector<NekDouble> freeParams)
    {
        TimeIntegrationSchemeSharedPtr p =
            MemoryManager<IMEXTimeIntegrationScheme>::AllocateSharedPtr(
                "", 2, freeParams);
        return p;
    }

    static std::string className;

protected:
    static std::string TimeIntegrationMethodLookupId;

}; // end class IMEXOrder2TimeIntegrationScheme

class IMEXOrder3TimeIntegrationScheme : public IMEXTimeIntegrationScheme
{
public:
    IMEXOrder3TimeIntegrationScheme(std::string variant, size_t order,
                                    std::vector<NekDouble> freeParams)
        : IMEXTimeIntegrationScheme("", 3, freeParams)
    {
        boost::ignore_unused(variant, order);
    }

    static TimeIntegrationSchemeSharedPtr create(
        [[maybe_unused]] std::string variant, [[maybe_unused]] size_t order,
        std::vector<NekDouble> freeParams)
    {
        TimeIntegrationSchemeSharedPtr p =
            MemoryManager<IMEXTimeIntegrationScheme>::AllocateSharedPtr(
                "", 3, freeParams);
        return p;
    }

    static std::string className;

protected:
    static std::string TimeIntegrationMethodLookupId;

}; // end class IMEXOrder3TimeIntegrationScheme

class IMEXOrder4TimeIntegrationScheme : public IMEXTimeIntegrationScheme
{
public:
    IMEXOrder4TimeIntegrationScheme(std::string variant, size_t order,
                                    std::vector<NekDouble> freeParams)
        : IMEXTimeIntegrationScheme("", 4, freeParams)
    {
        boost::ignore_unused(variant, order);
    }

    static TimeIntegrationSchemeSharedPtr create(
        [[maybe_unused]] std::string variant, [[maybe_unused]] size_t order,
        std::vector<NekDouble> freeParams)
    {

        TimeIntegrationSchemeSharedPtr p =
            MemoryManager<IMEXTimeIntegrationScheme>::AllocateSharedPtr(
                "", 4, freeParams);
        return p;
    }

    static std::string className;

protected:
    static std::string TimeIntegrationMethodLookupId;

}; // end class IMEXOrder4TimeIntegrationScheme

} // namespace Nektar::LibUtilities

#endif
