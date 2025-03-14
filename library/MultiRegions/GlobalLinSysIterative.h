///////////////////////////////////////////////////////////////////////////////
//
// File: GlobalLinSysIterative.h
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
// Description: GlobalLinSysIterative header
//
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_LIB_MULTIREGIONS_GLOBALLINSYSITERATIVE_H
#define NEKTAR_LIB_MULTIREGIONS_GLOBALLINSYSITERATIVE_H

#include <LibUtilities/LinearAlgebra/NekLinSysIterGMRES.h>
#include <MultiRegions/GlobalLinSys.h>
#include <MultiRegions/MultiRegionsDeclspec.h>
#include <MultiRegions/Preconditioner.h>

namespace Nektar::MultiRegions
{
// Forward declarations
class ExpList;

/// A global linear system.
class GlobalLinSysIterative : virtual public GlobalLinSys
{
public:
    /// Constructor for full direct matrix solve.
    MULTI_REGIONS_EXPORT GlobalLinSysIterative(
        const GlobalLinSysKey &pKey, const std::weak_ptr<ExpList> &pExpList,
        const std::shared_ptr<AssemblyMap> &pLocToGloMap);

    MULTI_REGIONS_EXPORT ~GlobalLinSysIterative() override;

    void DoMatrixMultiply(const Array<OneD, NekDouble> &pInput,
                          Array<OneD, NekDouble> &pOutput)
    {
        v_DoMatrixMultiply(pInput, pOutput);
    }

protected:
    /// Global to universal unique map
    Array<OneD, int> m_map;

    /// dot product of rhs to normalise stopping criterion
    NekDouble m_rhs_magnitude;

    /// cnt to how many times rhs_magnitude is called
    NekDouble m_rhs_mag_sm;

    PreconditionerSharedPtr m_precon;

    std::string m_precontype;

    int m_totalIterations;

    /// Whether to apply projection technique
    bool m_useProjection;

    /// Root if parallel
    bool m_root;

    /// Iterative solver: Conjugate Gradient, GMRES
    std::string m_linSysIterSolver;

    /// Storage for solutions to previous linear problems
    std::vector<Array<OneD, NekDouble>> m_prevLinSol;
    std::vector<Array<OneD, NekDouble>> m_prevBasis;
    DNekMatSharedPtr m_coeffMatrix;
    Array<OneD, NekDouble> m_coeffMatrixFactor;
    Array<OneD, int> m_ipivot;
    int m_numSuccessiveRHS;
    bool m_isAconjugate;
    std::string m_matrixType;
    bool m_isNonSymmetricLinSys;
    int m_numPrevSols;
    bool m_isAbsoluteTolerance;

    LibUtilities::NekSysOperators m_NekSysOp;

    LibUtilities::NekLinSysIterSharedPtr m_linsol;

    static std::string IteratSolverlookupIds[];
    static std::string IteratSolverdef;

    /// projection technique
    void DoProjection(const int pNumRows,
                      const Array<OneD, const NekDouble> &pInput,
                      Array<OneD, NekDouble> &pOutput, const int pNumDir,
                      const bool isAconjugate);

    virtual void v_UniqueMap() = 0;

    virtual void v_DoMatrixMultiply(const Array<OneD, NekDouble> &pInput,
                                    Array<OneD, NekDouble> &pOutput) = 0;

private:
    void UpdateKnownSolutions(const int pGlobalBndDofs,
                              const Array<OneD, const NekDouble> &pSolution,
                              const int pNumDirBndDofs,
                              const bool isAconjugate);

    int ResetKnownSolutionsToLatestOne();

    void DoPreconditionerFlag(const Array<OneD, NekDouble> &pInput,
                              Array<OneD, NekDouble> &pOutput,
                              const bool &isLocal = false)
    {
        m_precon->DoPreconditioner(pInput, pOutput, isLocal);
    }

    void DoMatrixMultiplyFlag(const Array<OneD, NekDouble> &pInput,
                              Array<OneD, NekDouble> &pOutput,
                              [[maybe_unused]] const bool &controlFlag)
    {
        v_DoMatrixMultiply(pInput, pOutput);
    }

    void DoAssembleLocFlag(const Array<OneD, NekDouble> &pInput,
                           Array<OneD, NekDouble> &pOutput, const bool &ZeroDir)
    {
        m_precon->DoAssembleLoc(pInput, pOutput, ZeroDir);
    }
};
} // namespace Nektar::MultiRegions

#endif
