///////////////////////////////////////////////////////////////////////////////
//
// File: TimeIntegrationScheme.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
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
// Description: Header file of time integration scheme base class,
// this class is the parent class for the TimeIntegrationSchemeGLM,
// TimeIntegrationSchemeFIT, TimeIntegrationSchemeGEM, and
// TimeIntegrationSchemeSDC classes.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_TIME_INTEGRATION_TIME_INTEGRATION_SCHEME
#define NEKTAR_LIB_UTILITIES_TIME_INTEGRATION_TIME_INTEGRATION_SCHEME

#include <string>

#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>

#include <LibUtilities/TimeIntegration/TimeIntegrationSchemeOperators.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationTypes.hpp>

#define LUE LIB_UTILITIES_EXPORT

namespace Nektar::LibUtilities
{

/// Datatype of the NekFactory used to instantiate classes derived from the
/// EquationSystem class.
typedef NekFactory<std::string, TimeIntegrationScheme, std::string, size_t,
                   std::vector<NekDouble>>
    TimeIntegrationSchemeFactory;

// Allows a code to create a TimeIntegrator. Usually used like this:
//
//    LibUtilities::TimeIntegrationSchemeSharedPtr timeIntegrationScheme =
//      LibUtilities::GetTimeIntegrationSchemeFactory().CreateInstance(
//                  "IMEX", "dirk", 3, std::vector<size_t>{2,3} );
LUE TimeIntegrationSchemeFactory &GetTimeIntegrationSchemeFactory();

/**
 * @brief Base class for time integration schemes.
 */
class TimeIntegrationScheme
{
public:
    // Access methods
    LUE std::string GetFullName() const
    {
        return v_GetFullName();
    }
    LUE std::string GetName() const
    {
        return v_GetName();
    }
    LUE std::string GetVariant() const
    {
        return v_GetVariant();
    }
    LUE size_t GetOrder() const
    {
        return v_GetOrder();
    }
    LUE std::vector<NekDouble> GetFreeParams()
    {
        return v_GetFreeParams();
    }
    LUE TimeIntegrationSchemeType GetIntegrationSchemeType()
    {
        return v_GetIntegrationSchemeType();
    }
    LUE NekDouble GetTimeStability() const
    {
        return v_GetTimeStability();
    }
    LUE size_t GetNumIntegrationPhases()
    {
        return v_GetNumIntegrationPhases();
    }

    /**
     * \brief Gets the solution vector of the ODE
     */
    LUE const TripleArray &GetSolutionVector() const
    {
        return v_GetSolutionVector();
    }

    LUE TripleArray &UpdateSolutionVector()
    {
        return v_UpdateSolutionVector();
    }

    /**
     * \brief Sets the solution vector of the ODE
     */
    LUE void SetSolutionVector(const size_t Offset, const DoubleArray &y)
    {
        v_SetSolutionVector(Offset, y);
    }

    // The worker methods
    /**
     * \brief Explicit integration of an ODE.
     *
     * This function explicitely perfroms a single integration step of the ODE
     * system:
     * \f[
     * \frac{d\boldsymbol{y}}{dt}=\boldsymbol{f}(t,\boldsymbol{y})
     * \f]
     *
     * \param timestep The size of the timestep, i.e. \f$\Delta t\f$.
     * \param f an object of the class FuncType, where FuncType should have a
     * method FuncType::ODEforcing
     *       to evaluate the right hand side \f$f(t,\boldsymbol{y})\f$ of the
     * ODE.
     * \param y on input: the vectors \f$\boldsymbol{y}^{[n-1]}\f$ and
     * \f$t^{[n-1]}\f$ (which corresponds to the
     *    solution at the old time level)
     * \param y on output:  the vectors \f$\boldsymbol{y}^{[n]}\f$ and
     * \f$t^{[n]}\f$ (which corresponds to the
     *    solution at the old new level)
     * \return The actual solution \f$\boldsymbol{y}^{n}\f$ at the new time
     * level
     *    (which in fact is also embedded in the argument y).
     */
    LUE void InitializeScheme(const NekDouble deltaT, ConstDoubleArray &y_0,
                              const NekDouble time,
                              const TimeIntegrationSchemeOperators &op)
    {
        v_InitializeScheme(deltaT, y_0, time, op);
    }

    LUE ConstDoubleArray &TimeIntegrate(const size_t timestep,
                                        const NekDouble delta_t)
    {
        return v_TimeIntegrate(timestep, delta_t);
    }

    LUE void print(std::ostream &os) const
    {
        v_print(os);
    }
    LUE void printFull(std::ostream &os) const
    {
        v_printFull(os);
    }

    // Friend classes
    LUE friend std::ostream &operator<<(std::ostream &os,
                                        const TimeIntegrationScheme &rhs);
    LUE friend std::ostream &operator<<(
        std::ostream &os, const TimeIntegrationSchemeSharedPtr &rhs);

protected:
    LUE virtual std::string v_GetFullName() const;
    LUE virtual std::string v_GetName() const                  = 0;
    LUE virtual std::string v_GetVariant() const               = 0;
    LUE virtual size_t v_GetOrder() const                      = 0;
    LUE virtual std::vector<NekDouble> v_GetFreeParams() const = 0;
    LUE virtual TimeIntegrationSchemeType v_GetIntegrationSchemeType()
        const                                                  = 0;
    LUE virtual NekDouble v_GetTimeStability() const           = 0;
    LUE virtual size_t v_GetNumIntegrationPhases() const       = 0;
    LUE virtual const TripleArray &v_GetSolutionVector() const = 0;
    LUE virtual TripleArray &v_UpdateSolutionVector()          = 0;
    LUE virtual void v_SetSolutionVector(const size_t Offset,
                                         const DoubleArray &y) = 0;
    LUE virtual void v_InitializeScheme(
        const NekDouble deltaT, ConstDoubleArray &y_0, const NekDouble time,
        const TimeIntegrationSchemeOperators &op)                          = 0;
    LUE virtual ConstDoubleArray &v_TimeIntegrate(const size_t timestep,
                                                  const NekDouble delta_t) = 0;
    LUE virtual void v_print(std::ostream &os) const                       = 0;
    LUE virtual void v_printFull(std::ostream &os) const                   = 0;

    // These methods should never be used directly, only used by child classes.
    LUE TimeIntegrationScheme(std::string variant, size_t order,
                              std::vector<NekDouble> freeParams);

    LUE TimeIntegrationScheme(const TimeIntegrationScheme &in) = delete;
    virtual ~TimeIntegrationScheme()                           = default;

}; // end class TimeIntegrationScheme

LUE std::ostream &operator<<(std::ostream &os,
                             const TimeIntegrationScheme &rhs);
LUE std::ostream &operator<<(std::ostream &os,
                             const TimeIntegrationSchemeSharedPtr &rhs);

} // namespace Nektar::LibUtilities

#endif
