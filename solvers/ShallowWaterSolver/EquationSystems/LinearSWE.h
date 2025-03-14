///////////////////////////////////////////////////////////////////////////////
//
// File: LinearSWE.h
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
// Description: Linear Shallow water equations in primitive variables
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_SHALLOWWATERSOLVER_EQUATIONSYSTEMS_LINEARSWE_H
#define NEKTAR_SOLVERS_SHALLOWWATERSOLVER_EQUATIONSYSTEMS_LINEARSWE_H

#include <ShallowWaterSolver/EquationSystems/ShallowWaterSystem.h>

namespace Nektar
{

class LinearSWE : public ShallowWaterSystem
{
public:
    friend class MemoryManager<LinearSWE>;

    /// Creates an instance of this class
    static SolverUtils::EquationSystemSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const SpatialDomains::MeshGraphSharedPtr &pGraph)
    {
        SolverUtils::EquationSystemSharedPtr p =
            MemoryManager<LinearSWE>::AllocateSharedPtr(pSession, pGraph);
        p->InitObject();
        return p;
    }

    /// Name of class
    static std::string className;

    ~LinearSWE() override = default;

protected:
    LinearSWE(const LibUtilities::SessionReaderSharedPtr &pSession,
              const SpatialDomains::MeshGraphSharedPtr &pGraph);

    void v_InitObject(bool DeclareFields = true) override;

    void v_DoOdeRhs(const Array<OneD, const Array<OneD, NekDouble>> &inarray,
                    Array<OneD, Array<OneD, NekDouble>> &outarray,
                    const NekDouble time) override;

    void v_GenerateSummary(SolverUtils::SummaryList &s) override;

    void GetFluxVector(
        const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &flux);

    void GetVelocityVector(
        const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        Array<OneD, Array<OneD, NekDouble>> &velocity);

    void CopyBoundaryTrace(const Array<OneD, const NekDouble> &Fwd,
                           Array<OneD, NekDouble> &Bwd);

    const Array<OneD, NekDouble> &GetDepthFwd()
    {
        return m_dFwd;
    }

    const Array<OneD, NekDouble> &GetDepthBwd()
    {
        return m_dBwd;
    }

private:
    /// Still water depth traces
    Array<OneD, NekDouble> m_dFwd;
    Array<OneD, NekDouble> m_dBwd;
};

} // namespace Nektar

#endif
