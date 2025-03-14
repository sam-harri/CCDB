///////////////////////////////////////////////////////////////////////////////
//
// File: EulerImplicitCFE.h
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
// Description: EulerImplicit equations in conservative variables without
// artificial diffusion
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_EULERIMPLICITCFE_H
#define NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_EULERIMPLICITCFE_H

#include <CompressibleFlowSolver/EquationSystems/CompressibleFlowSystemImplicit.h>

namespace Nektar
{

class EulerImplicitCFE : public CFSImplicit
{
public:
    friend class MemoryManager<EulerImplicitCFE>;

    /// Creates an instance of this class.
    static SolverUtils::EquationSystemSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const SpatialDomains::MeshGraphSharedPtr &pGraph)
    {
        SolverUtils::EquationSystemSharedPtr p =
            MemoryManager<EulerImplicitCFE>::AllocateSharedPtr(pSession,
                                                               pGraph);
        p->InitObject();
        return p;
    }

    /// Name of class.
    static std::string className;

    ~EulerImplicitCFE() override = default;

protected:
    EulerImplicitCFE(const LibUtilities::SessionReaderSharedPtr &pSession,
                     const SpatialDomains::MeshGraphSharedPtr &pGraph);

    void v_InitObject(bool DeclareFields = true) override;

    void v_DoDiffusion(
        [[maybe_unused]] const Array<OneD, Array<OneD, NekDouble>> &inarray,
        [[maybe_unused]] Array<OneD, Array<OneD, NekDouble>> &outarray,
        [[maybe_unused]] const Array<OneD, Array<OneD, NekDouble>> &pFwd,
        [[maybe_unused]] const Array<OneD, Array<OneD, NekDouble>> &pBwd) final
    {
        NEKERROR(ErrorUtil::efatal,
                 "v_DoDiffusion is not implemented for implicit solvers");
    }

    bool v_SupportsShockCaptType(const std::string type) const final;
};
} // namespace Nektar
#endif
