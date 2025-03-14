///////////////////////////////////////////////////////////////////////////////
//
// File: PreconditionerDiagonal.cpp
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
// Description: PreconditionerDiagonal definition
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/GlobalLinSysIterative.h>
#include <MultiRegions/GlobalLinSysIterativeStaticCond.h>
#include <MultiRegions/GlobalMatrixKey.h>
#include <MultiRegions/PreconditionerDiagonal.h>
#include <cmath>

using namespace std;

namespace Nektar::MultiRegions
{
/**
 * Registers the class with the Factory.
 */
string PreconditionerDiagonal::className =
    GetPreconFactory().RegisterCreatorFunction(
        "Diagonal", PreconditionerDiagonal::create, "Diagonal Preconditioning");
/**
 * @class Preconditioner
 *
 * This class implements diagonal preconditioning for the conjugate
 * gradient matrix solver.
 */

PreconditionerDiagonal::PreconditionerDiagonal(
    const std::shared_ptr<GlobalLinSys> &plinsys,
    const AssemblyMapSharedPtr &pLocToGloMap)
    : Preconditioner(plinsys, pLocToGloMap)
{
}

void PreconditionerDiagonal::v_InitObject()
{
}

void PreconditionerDiagonal::v_BuildPreconditioner()
{
    GlobalSysSolnType solvertype = m_locToGloMap.lock()->GetGlobalSysSolnType();
    if (solvertype == eIterativeFull)
    {
        DiagonalPreconditionerSum();
    }
    else if (solvertype == eIterativeStaticCond ||
             solvertype == eIterativeMultiLevelStaticCond ||
             solvertype == ePETScStaticCond ||
             solvertype == ePETScMultiLevelStaticCond)
    {
        StaticCondDiagonalPreconditionerSum();
    }
    else
    {
        ASSERTL0(0, "Unsupported solver type");
    }
}

/**
 * Diagonal preconditioner computed by summing the relevant elements of
 * the local matrix system.
 */
void PreconditionerDiagonal::DiagonalPreconditionerSum()
{
    std::shared_ptr<MultiRegions::ExpList> expList =
        ((m_linsys.lock())->GetLocMat()).lock();

    auto asmMap = m_locToGloMap.lock();
    int nGlobal = asmMap->GetNumGlobalCoeffs();
    int nDir    = asmMap->GetNumGlobalDirBndCoeffs();
    int nInt    = nGlobal - nDir;
    int nElmt   = expList->GetNumElmts();

    Array<OneD, NekDouble> vOutput(nGlobal, 0.0);
    for (int n = 0, cnt = 0; n < nElmt; ++n)
    {
        auto loc_mat = (m_linsys.lock())->GetBlock(n);
        int loc_row  = loc_mat->GetRows();
        for (int i = 0; i < loc_row; ++i)
        {
            int gid1 = asmMap->GetLocalToGlobalMap(cnt + i);
            vOutput[gid1] += (*loc_mat)(i, i);
        }
        cnt += loc_row;
    }

    // Assemble diagonal contributions across processes
    asmMap->UniversalAssemble(vOutput);

    m_diagonals = Array<OneD, NekDouble>(nInt);
    Vmath::Sdiv(nInt, 1.0, &vOutput[nDir], 1, &m_diagonals[0], 1);
}

/**
 * Diagonal preconditioner defined as the inverse of the main
 * diagonal of the Schur complement
 *
 */
void PreconditionerDiagonal::StaticCondDiagonalPreconditionerSum()
{
    auto asmMap = m_locToGloMap.lock();

    int nGlobalBnd = asmMap->GetNumGlobalBndCoeffs();
    int nDirBnd    = asmMap->GetNumGlobalDirBndCoeffs();
    int rows       = nGlobalBnd - nDirBnd;

    Array<OneD, NekDouble> vOutput(nGlobalBnd, 0.0);

    // Extract diagonal contributions
    Array<OneD, NekDouble> diagonals = AssembleStaticCondGlobalDiagonals();
    for (unsigned int i = 0; i < rows; ++i)
    {
        vOutput[nDirBnd + i] = diagonals[i];
    }

    // Assemble diagonal contributions across processes
    asmMap->UniversalAssembleBnd(vOutput);

    m_diagonals = Array<OneD, NekDouble>(rows);
    Vmath::Sdiv(rows, 1.0, &vOutput[nDirBnd], 1, &m_diagonals[0], 1);
}

/**
 *
 */
void PreconditionerDiagonal::v_DoPreconditioner(
    const Array<OneD, NekDouble> &pInput, Array<OneD, NekDouble> &pOutput,
    const bool &IsLocal)
{
    auto asmMap = m_locToGloMap.lock();

    GlobalSysSolnType solvertype = asmMap->GetGlobalSysSolnType();

    bool isFull = solvertype == eIterativeFull ? true : false;

    int nGlobal = (isFull) ? asmMap->GetNumGlobalCoeffs()
                           : asmMap->GetNumGlobalBndCoeffs();
    int nDir    = asmMap->GetNumGlobalDirBndCoeffs();
    int nNonDir = nGlobal - nDir;

    if (IsLocal)
    {
        Array<OneD, NekDouble> wk(nGlobal);
        (isFull) ? asmMap->Assemble(pInput, wk)
                 : asmMap->AssembleBnd(pInput, wk);
        Vmath::Vmul(nNonDir, wk.data() + nDir, 1, m_diagonals.data(), 1,
                    wk.data() + nDir, 1);
        Vmath::Zero(nDir, wk, 1);
        (isFull) ? asmMap->GlobalToLocal(wk, pOutput)
                 : asmMap->GlobalToLocalBnd(wk, pOutput);
    }
    else
    {
        Vmath::Vmul(nNonDir, &pInput[0], 1, &m_diagonals[0], 1, &pOutput[0], 1);
    }
}

string PreconditionerNull::className =
    GetPreconFactory().RegisterCreatorFunction(
        "Null", PreconditionerNull::create, "No Preconditioning");

/**
 * @class Null Preconditioner
 *
 * This class implements no preconditioning for the conjugate
 * gradient matrix solver.
 */
PreconditionerNull::PreconditionerNull(
    const std::shared_ptr<GlobalLinSys> &plinsys,
    const AssemblyMapSharedPtr &pLocToGloMap)
    : Preconditioner(plinsys, pLocToGloMap)
{
}

/**
 *
 */
void PreconditionerNull::v_InitObject()
{
}

/**
 *
 */
void PreconditionerNull::v_BuildPreconditioner()
{
}

/**
 *
 */
void PreconditionerNull::v_DoPreconditioner(
    const Array<OneD, NekDouble> &pInput, Array<OneD, NekDouble> &pOutput,
    const bool &isLocal)

{
    auto asmMap = m_locToGloMap.lock();

    GlobalSysSolnType solvertype = asmMap->GetGlobalSysSolnType();

    bool isFull = solvertype == eIterativeFull ? true : false;

    int nGlobal = (isFull) ? asmMap->GetNumGlobalCoeffs()
                           : asmMap->GetNumGlobalBndCoeffs();
    int nDir    = asmMap->GetNumGlobalDirBndCoeffs();

    if (isLocal)
    {
        Array<OneD, NekDouble> wk(nGlobal);
        (isFull) ? asmMap->Assemble(pInput, wk)
                 : asmMap->AssembleBnd(pInput, wk);
        Vmath::Zero(nDir, wk, 1);
        (isFull) ? asmMap->GlobalToLocal(wk, pOutput)
                 : asmMap->GlobalToLocalBnd(wk, pOutput);
    }
    else
    {
        Vmath::Vcopy(pInput.size(), pInput, 1, pOutput, 1);
    }
}

/**
 * Registers the class with the Factory.
 */
string PreconditionerJacobi::className =
    GetPreconFactory().RegisterCreatorFunction(
        "Jacobi", PreconditionerJacobi::create, "Jacobi Preconditioning");
/**
 * @class PreconditionerJacobi
 *
 * This class implements jacobi preconditioning  building on diagonal version
 * above
 */

PreconditionerJacobi::PreconditionerJacobi(
    const std::shared_ptr<GlobalLinSys> &plinsys,
    const AssemblyMapSharedPtr &pLocToGloMap)
    : PreconditionerDiagonal(plinsys, pLocToGloMap)
{
}

void PreconditionerJacobi::v_InitObject()
{
}

void PreconditionerJacobi::v_BuildPreconditioner()
{
    PreconditionerDiagonal::v_BuildPreconditioner();

    auto expList = ((m_linsys.lock())->GetLocMat()).lock();
    std::shared_ptr<LibUtilities::SessionReader> session =
        expList->GetSession();

    std::string var = m_locToGloMap.lock()->GetVariable();

    if (session->DefinesGlobalSysSolnInfo(var, "JacobiIterations"))
    {
        m_niter = boost::lexical_cast<int>(
            session->GetGlobalSysSolnInfo(var, "JacobiIterations").c_str());
    }
    else
    {
        m_niter = 3;
    }
}

/**
 *
 */
void PreconditionerJacobi::v_DoPreconditioner(
    const Array<OneD, NekDouble> &pInput, Array<OneD, NekDouble> &pOutput,
    const bool &IsLocal)
{
    auto asmMap = m_locToGloMap.lock();

    GlobalSysSolnType solvertype = asmMap->GetGlobalSysSolnType();

    int nGlobal = solvertype == eIterativeFull
                      ? asmMap->GetNumGlobalCoeffs()
                      : asmMap->GetNumGlobalBndCoeffs();
    int nDir    = asmMap->GetNumGlobalDirBndCoeffs();
    int nNonDir = nGlobal - nDir;
    Array<OneD, NekDouble> wk(nGlobal);

    if (IsLocal)
    {
        int nLocal = solvertype == eIterativeFull
                         ? asmMap->GetNumLocalCoeffs()
                         : asmMap->GetNumLocalBndCoeffs();

        Array<OneD, NekDouble> wk1(nLocal);
        asmMap->Assemble(pInput, wk);
        Vmath::Vmul(nNonDir, wk.data() + nDir, 1, m_diagonals.data(), 1,
                    wk.data() + nDir, 1);
        Vmath::Zero(nDir, wk, 1);

        for (int n = 1; n < m_niter; ++n)
        {
            asmMap->GlobalToLocal(wk, pOutput);

            // do Ax operator
            std::dynamic_pointer_cast<GlobalLinSysIterative>(m_linsys.lock())
                ->DoMatrixMultiply(pOutput, wk1);

            Vmath::Vsub(nLocal, pInput, 1, wk1, 1, wk1, 1);

            asmMap->Assemble(wk1, pOutput);
            Vmath::Vvtvp(nNonDir, pOutput.data() + nDir, 1, m_diagonals.data(),
                         1, wk.data() + nDir, 1, wk.data() + nDir, 1);
        }

        asmMap->GlobalToLocal(wk, pOutput);
    }
    else
    {
        Array<OneD, NekDouble> wk1(nGlobal);
        Vmath::Vmul(nNonDir, pInput.data(), 1, m_diagonals.data(), 1,
                    wk.data() + nDir, 1);
        Vmath::Zero(nDir, wk, 1);

        for (int n = 1; n < m_niter; ++n)
        {
            // do Ax operator
            std::dynamic_pointer_cast<GlobalLinSysIterative>(m_linsys.lock())
                ->DoMatrixMultiply(wk, wk1);

            // b - Ax
            Vmath::Vsub(nNonDir, pInput.data(), 1, wk1.data() + nDir, 1,
                        wk1.data() + nDir, 1);

            // new sol = 1/diag (b-Ax) + old sol
            Vmath::Vvtvp(nNonDir, wk1.data() + nDir, 1, m_diagonals.data(), 1,
                         wk.data() + nDir, 1, wk.data() + nDir, 1);
        }

        Vmath::Vcopy(nNonDir, wk.data() + nDir, 1, pOutput.data(), 1);
    }
}

} // namespace Nektar::MultiRegions
