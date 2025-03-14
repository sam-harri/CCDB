///////////////////////////////////////////////////////////////////////////////
//
// File: UnsteadyAdvectionDiffusion.h
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
// Description: Unsteady advection-diffusion solve routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_UNSTEADYADVECTIONDIFFUSION_H
#define NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_UNSTEADYADVECTIONDIFFUSION_H

#include <ADRSolver/EquationSystems/UnsteadyAdvection.h>
#include <SolverUtils/Diffusion/Diffusion.h>

namespace Nektar
{
class UnsteadyAdvectionDiffusion : public UnsteadyAdvection
{
public:
    friend class MemoryManager<UnsteadyAdvectionDiffusion>;

    /// Creates an instance of this class
    static SolverUtils::EquationSystemSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const SpatialDomains::MeshGraphSharedPtr &pGraph)
    {
        SolverUtils::EquationSystemSharedPtr p =
            MemoryManager<UnsteadyAdvectionDiffusion>::AllocateSharedPtr(
                pSession, pGraph);
        p->InitObject();
        return p;
    }

    /// Name of class
    static std::string className;

    /// Destructor
    ~UnsteadyAdvectionDiffusion() override = default;

    void v_ALEInitObject(
        int spaceDim,
        Array<OneD, MultiRegions::ExpListSharedPtr> &fields) override;

protected:
    bool m_subSteppingScheme;

    // Use Spectral Vanishing Viscosity
    bool m_useSpecVanVisc;
    // cut off ratio from which to start decayhing modes
    NekDouble m_sVVCutoffRatio;
    // Diffusion coefficient of SVV modes
    NekDouble m_sVVDiffCoeff;

    SolverUtils::DiffusionSharedPtr m_diffusion;

    /// Session reader
    UnsteadyAdvectionDiffusion(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const SpatialDomains::MeshGraphSharedPtr &pGraph);

    /// Evaluate the flux at each solution point for the diffusion part
    void GetFluxVectorDiff(
        const Array<OneD, Array<OneD, NekDouble>> &inarray,
        const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &qfield,
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &viscousTensor);

    /// Compute the RHS
    void DoOdeRhs(const Array<OneD, const Array<OneD, NekDouble>> &inarray,
                  Array<OneD, Array<OneD, NekDouble>> &outarray,
                  const NekDouble time);

    /// Solve implicitly the diffusion term
    void DoImplicitSolve(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, NekDouble time,
        NekDouble lambda);

    /// Initialise the object
    void v_InitObject(bool DeclareFields = true) override;

    /// Print Summary
    void v_GenerateSummary(SolverUtils::SummaryList &s) override;

    /// PreIntegration step for substepping.
    bool v_PreIntegrate(int step) override;

    void v_ExtraFldOutput(std::vector<Array<OneD, NekDouble>> &fieldcoeffs,
                          std::vector<std::string> &variables) override;

    // SubsStepping methods -> Probably could be set up in separate class
    void SubStepAdvance(int nstep, NekDouble time);

    NekDouble GetSubstepTimeStep();

    void SetUpSubSteppingTimeIntegration(
        const LibUtilities::TimeIntegrationSchemeSharedPtr &IntegrationScheme);

    void SubStepAdvection(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time);

    void SubStepProjection(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time);

    void AddAdvectionPenaltyFlux(
        const Array<OneD, const Array<OneD, NekDouble>> &velfield,
        const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        Array<OneD, Array<OneD, NekDouble>> &Outarray);

    Array<OneD, NekDouble> GetMaxStdVelocity(
        const Array<OneD, Array<OneD, NekDouble>> inarray);

    LibUtilities::TimeIntegrationSchemeSharedPtr m_subStepIntegrationScheme;
    LibUtilities::TimeIntegrationSchemeOperators m_subStepIntegrationOps;

    int m_intSteps;
    int m_minsubsteps;

private:
    NekDouble m_epsilon;
};
} // namespace Nektar

#endif
