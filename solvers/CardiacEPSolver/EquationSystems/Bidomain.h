///////////////////////////////////////////////////////////////////////////////
//
// File: Bidomain.h
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
// Description: Bidomain cardiac electrophysiology homogenised model.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_BIDOMAIN_H
#define NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_BIDOMAIN_H

#include <CardiacEPSolver/CellModels/CellModel.h>
#include <SolverUtils/UnsteadySystem.h>

using namespace Nektar::SolverUtils;

namespace Nektar
{

/// A model for cardiac conduction.
class Bidomain : public UnsteadySystem
{
public:
    friend class MemoryManager<Bidomain>;

    /// Creates an instance of this class
    static EquationSystemSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const SpatialDomains::MeshGraphSharedPtr &pGraph)
    {
        EquationSystemSharedPtr p =
            MemoryManager<Bidomain>::AllocateSharedPtr(pSession, pGraph);
        p->InitObject();
        return p;
    }

    /// Name of class
    static std::string className;

    /// Desctructor
    ~Bidomain() override;

protected:
    /// Constructor
    Bidomain(const LibUtilities::SessionReaderSharedPtr &pSession,
             const SpatialDomains::MeshGraphSharedPtr &pGraph);

    void v_InitObject(bool DeclareField = true) override;

    /// Solve for the diffusion term.
    void DoImplicitSolve(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, NekDouble time,
        NekDouble lambda);

    /// Computes the reaction terms \f$f(u,v)\f$ and \f$g(u,v)\f$.
    void DoOdeRhs(const Array<OneD, const Array<OneD, NekDouble>> &inarray,
                  Array<OneD, Array<OneD, NekDouble>> &outarray,
                  const NekDouble time);

    /// Sets a custom initial condition.
    void v_SetInitialConditions(NekDouble initialtime,
                                bool dumpInitialConditions,
                                const int domain) override;

    /// Prints a summary of the model parameters.
    void v_GenerateSummary(SummaryList &s) override;

private:
    /// Cell model.
    CellModelSharedPtr m_cell;

    NekDouble m_chi, m_capMembrane, m_sigmaix, m_sigmaiy, m_sigmaiz, m_sigmaex,
        m_sigmaey, m_sigmaez;

    StdRegions::VarCoeffMap m_vardiffi;
    StdRegions::VarCoeffMap m_vardiffie;

    Array<OneD, Array<OneD, NekDouble>> tmp1;
    Array<OneD, Array<OneD, NekDouble>> tmp2;
    Array<OneD, Array<OneD, NekDouble>> tmp3;

    /// Stimulus current
    NekDouble m_stimDuration;
};

} // namespace Nektar

#endif
