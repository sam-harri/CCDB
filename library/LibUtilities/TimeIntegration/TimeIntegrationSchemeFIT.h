///////////////////////////////////////////////////////////////////////////////
//
// File: TimeIntegrationSchemeFIT.h
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
// Description: Header file of time integration scheme FIT base class
//
///////////////////////////////////////////////////////////////////////////////

// Note: The file is named TimeIntegrationSchemeFIT to parallel the
// TimeIntegrationSchemeGLM file but the class is named
// FractionalInTimeIntegrationScheme so keep with the factory naming
// convention.

#ifndef NEKTAR_LIB_UTILITIES_TIME_INTEGRATION_TIME_INTEGRATION_SCHEME_FIT
#define NEKTAR_LIB_UTILITIES_TIME_INTEGRATION_TIME_INTEGRATION_SCHEME_FIT

#define LUE LIB_UTILITIES_EXPORT

#include <string>

#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>

namespace Nektar::LibUtilities
{

///////////////////////////////////////////////////////////////////////////////
/// Class for fractional-in-time integration.
class FractionalInTimeIntegrationScheme : public TimeIntegrationScheme
{
public:
    /// Constructor
    FractionalInTimeIntegrationScheme(std::string variant, size_t order,
                                      std::vector<NekDouble> freeParams);

    /// Destructor
    ~FractionalInTimeIntegrationScheme() override
    {
    }

    /// Creator
    static TimeIntegrationSchemeSharedPtr create(
        std::string variant, size_t order, std::vector<NekDouble> freeParams)
    {
        TimeIntegrationSchemeSharedPtr p =
            MemoryManager<FractionalInTimeIntegrationScheme>::AllocateSharedPtr(
                variant, order, freeParams);

        return p;
    }

    static std::string className;

    // Friend classes
    LUE friend std::ostream &operator<<(
        std::ostream &os, const FractionalInTimeIntegrationScheme &rhs);
    LUE friend std::ostream &operator<<(
        std::ostream &os,
        const FractionalInTimeIntegrationSchemeSharedPtr &rhs);

protected:
    // Access methods from the base class that are virtual
    LUE std::string v_GetName() const override
    {
        return m_name;
    }

    LUE std::string v_GetVariant() const override
    {
        return m_variant;
    }

    LUE size_t v_GetOrder() const override
    {
        return m_order;
    }

    LUE std::vector<NekDouble> v_GetFreeParams() const override
    {
        return m_freeParams;
    }

    LUE TimeIntegrationSchemeType v_GetIntegrationSchemeType() const override
    {
        return m_schemeType;
    }

    LUE NekDouble v_GetTimeStability() const override
    {
        return 1.0;
    }

    LUE size_t v_GetNumIntegrationPhases() const override
    {
        return 1;
    }

    /**
     * \brief Gets the solution vector of the ODE
     */
    const TripleArray &v_GetSolutionVector() const override
    {
        return m_u;
    }
    TripleArray &v_UpdateSolutionVector() override
    {
        return m_u;
    }

    /**
     * \brief Sets the solution vector of the ODE
     */
    void v_SetSolutionVector(const size_t Offset, const DoubleArray &y) override
    {
        m_u[Offset] = y;
    }

    // The worker methods from the base class that are virtual
    LUE void v_InitializeScheme(
        const NekDouble deltaT, ConstDoubleArray &y_0, const NekDouble time,
        const TimeIntegrationSchemeOperators &op) override;

    LUE ConstDoubleArray &v_TimeIntegrate(const size_t timestep,
                                          const NekDouble delta_t) override;

    LUE void v_print(std::ostream &os) const override;
    LUE void v_printFull(std::ostream &os) const override;

    struct Instance
    {
        size_t base;

        size_t index;         // Index of this instance
        bool active;          // Used to determine if active
        size_t activecounter; // counter used to flip active bit
        size_t activebase;

        // Major storage for auxilliary ODE solutions.
        // Storage for values of y currently used to update u
        ComplexTripleArray stage_y;
        std::pair<size_t, size_t>
            stage_ind; // Time-step counters indicating the
                       // interval ymain is associated with

        // Staging allocation
        bool stage_active;
        size_t stage_ccounter;
        size_t stage_cbase; // This base is halved after the first cycle
        size_t stage_fcounter;
        size_t stage_fbase; // This base is halved after the first cycle

        // Ceiling stash allocation
        size_t cstash_counter; // Counter used to determine
                               // when to stash
        size_t cstash_base;    // base for counter
        ComplexTripleArray cstash_y;
        std::pair<size_t, size_t> cstash_ind; // ind(1) is never used:
                                              // it always matches main.ind(1)

        // Ceiling sandbox allocation
        bool csandbox_active; // Flag to determine when
                              // stash 2 is utilized
        size_t csandbox_counter;
        ComplexTripleArray csandbox_y;
        std::pair<size_t, size_t> csandbox_ind;

        // Floor stash
        size_t fstash_base;
        ComplexTripleArray fstash_y;
        std::pair<size_t, size_t> fstash_ind;

        // Floor sandbox
        bool fsandbox_active;
        size_t fsandbox_activebase;
        size_t fsandbox_stashincrement;
        ComplexTripleArray fsandbox_y;
        std::pair<size_t, size_t> fsandbox_ind;

        // Talbot quadrature rule
        ComplexSingleArray z;
        ComplexSingleArray w;

        TripleArray As;

        ComplexSingleArray E;
        ComplexDoubleArray Eh;
        ComplexDoubleArray AtEh;
    };

    inline size_t modIncrement(const size_t counter, const size_t base) const;

    inline size_t computeL(const size_t base, const size_t m) const;

    inline size_t computeQML(const size_t base, const size_t m);

    inline size_t computeTaus(const size_t base, const size_t m);

    void talbotQuadrature(const size_t nQuadPts, const NekDouble mu,
                          const NekDouble nu, const NekDouble sigma,
                          ComplexSingleArray &lamb,
                          ComplexSingleArray &w) const;

    void integralClassInitialize(const size_t index, Instance &instance) const;

    void updateStage(const size_t timeStep, Instance &instance);

    void finalIncrement(const size_t timeStep);

    void integralContribution(const size_t timeStep, const size_t tauml,
                              const Instance &instance);

    void timeAdvance(const size_t timeStep, Instance &instance,
                     ComplexTripleArray &y);

    void advanceSandbox(const size_t timeStep, Instance &instance);

    // Variables common to all schemes.
    std::string m_name;
    std::string m_variant;
    size_t m_order{0};
    std::vector<NekDouble> m_freeParams;

    TimeIntegrationSchemeType m_schemeType{eFractionalInTime};
    TimeIntegrationSchemeOperators m_op;

    // Varaibles and methods specific to FIT integration schemes.
    NekDouble m_deltaT{0};
    NekDouble m_T{0};       // Finial time
    size_t m_maxTimeSteps;  // Number of time steps.
    NekDouble m_alpha{0.3}; // Value for exp integration.
    size_t m_base{4};       // "Base" of the algorithm.
    size_t m_nQuadPts{20};  // Number of Talbot quadrature rule points
    NekDouble m_sigma{0};
    NekDouble m_mu0{8};
    NekDouble m_nu{0.6};

    size_t m_nvars{0};   // Number of variables in the integration scheme.
    size_t m_npoints{0}; // Number of points    in the integration scheme.

    size_t m_Lmax{0}; // Maxium number of integral groups.
    Array<OneD, Instance> m_integral_classes;

    // Demarcation integers
    Array<OneD, size_t> m_qml;
    // Demarcation interval markers
    Array<OneD, size_t> m_taus;

    // Storage of the initial values.
    DoubleArray m_u0;
    // Storage of the next solution from the final increment.
    DoubleArray m_uNext;
    // Storage for the integral contribution.
    ComplexDoubleArray m_uInt;
    // Storage for the exponential factor in the integral contribution.
    ComplexSingleArray m_expFactor;

    // Storage of previous states and associated timesteps.
    TripleArray m_u;

    // Storage for the stage derivative as the data will be re-used to
    // update the solution.
    DoubleArray m_F;

    // J
    SingleArray m_J;

    // Ahat array one for each order.
    TripleArray m_Ahats;

    // Multiply the last Ahat array, transposed by J
    SingleArray m_AhattJ;

}; // end class FractionalInTimeIntegrator

LUE std::ostream &operator<<(std::ostream &os,
                             const FractionalInTimeIntegrationScheme &rhs);
LUE std::ostream &operator<<(
    std::ostream &os, const FractionalInTimeIntegrationSchemeSharedPtr &rhs);

} // namespace Nektar::LibUtilities

#endif
