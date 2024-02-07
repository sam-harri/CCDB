///////////////////////////////////////////////////////////////////////////////
//
// File: NekNonlinSysNewton.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
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
// Description: NekNonlinSysNewton definition
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/LinearAlgebra/NekNonlinSysNewton.h>

using namespace std;

namespace Nektar::LibUtilities
{
/**
 * @class  NekNonlinSysNewton
 *
 * Solves a nonlinear system using iterative methods.
 */
string NekNonlinSysNewton::className =
    LibUtilities::GetNekNonlinSysFactory().RegisterCreatorFunction(
        "Newton", NekNonlinSysNewton::create, "NekNonlinSysNewton solver.");

NekNonlinSysNewton::NekNonlinSysNewton(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const LibUtilities::CommSharedPtr &vRowComm, const int nscale,
    const NekSysKey &pKey)
    : NekNonlinSys(pSession, vRowComm, nscale, pKey)
{
}

void NekNonlinSysNewton::v_InitObject()
{
    NekSys::v_InitObject();
    m_Residual = Array<OneD, NekDouble>(m_SysDimen, 0.0);
    m_DeltSltn = Array<OneD, NekDouble>(m_SysDimen, 0.0);
}

void NekNonlinSysNewton::v_SetSysOperators(const NekSysOperators &in)
{
    NekSys::v_SetSysOperators(in);
    m_linsol->SetSysOperators(in);
}

/**
 *
 **/
int NekNonlinSysNewton::v_SolveSystem(
    const int nGlobal, const Array<OneD, const NekDouble> &pInput,
    Array<OneD, NekDouble> &pOutput, const int nDir, const NekDouble tol,
    [[maybe_unused]] const NekDouble factor)
{
    ASSERTL0(nDir == 0, "nDir != 0 not tested");
    ASSERTL0(nGlobal == m_SysDimen, "nGlobal != m_SysDimen");

    int ntotal = nGlobal - nDir;

    m_SourceVec     = pInput;
    m_Solution      = pOutput;
    m_NtotLinSysIts = 0;
    m_converged     = false;

    Vmath::Vcopy(ntotal, pInput, 1, m_Solution, 1);

    NekDouble resnormOld = 0.0;
    int NttlNonlinIte    = 0;
    for (; NttlNonlinIte < m_maxiter; ++NttlNonlinIte)
    {
        m_operator.DoNekSysResEval(m_Solution, m_Residual);

        m_converged = v_ConvergenceCheck(NttlNonlinIte, m_Residual, tol);
        if (m_converged)
        {
            break;
        }

        NekDouble LinSysRelativeIteTol =
            CalcInexactNewtonForcing(NttlNonlinIte, resnormOld, m_SysResNorm);
        resnormOld = m_SysResNorm;
        m_linsol->setRhsMagnitude(m_SysResNorm);
        int ntmpLinSysIts = m_linsol->SolveSystem(
            ntotal, m_Residual, m_DeltSltn, 0, LinSysRelativeIteTol);
        m_NtotLinSysIts += ntmpLinSysIts;

        Vmath::Vsub(ntotal, m_Solution, 1, m_DeltSltn, 1, m_Solution, 1);
    }

    if ((!m_converged || m_verbose) && m_root && m_FlagWarnings)
    {
        int nwidthcolm = 11;

        WARNINGL0(m_converged,
                  "     # Nonlinear solver not converge in DoImplicitSolve");
        cout << right << scientific << setw(nwidthcolm)
             << setprecision(nwidthcolm - 6)
             << "     * Newton-Its converged (RES=" << sqrt(m_SysResNorm)
             << " Res/(DtRHS): " << sqrt(m_SysResNorm / m_SysResNorm0)
             << " with " << setw(3) << NttlNonlinIte << " Non-Its)" << endl;
    }

    return NttlNonlinIte;
}

bool NekNonlinSysNewton::v_ConvergenceCheck(
    const int nIteration, const Array<OneD, const NekDouble> &Residual,
    const NekDouble tol)
{
    m_SysResNorm = Vmath::Dot(Residual.size(), Residual, Residual);
    m_rowComm->AllReduce(m_SysResNorm, Nektar::LibUtilities::ReduceSum);

    if (nIteration == 0)
    {
        m_SysResNorm0 = m_SysResNorm;
    }

    NekDouble resratio = m_SysResNorm / m_SysResNorm0;
    NekDouble restol   = m_NonlinIterTolRelativeL2;

    return resratio < restol * restol || m_SysResNorm < tol * tol;
}

NekDouble NekNonlinSysNewton::CalcInexactNewtonForcing(
    const int &nIteration, const NekDouble &resnormOld,
    const NekDouble &resnorm)
{
    if (nIteration == 0 || !m_InexactNewtonForcing)
    {
        return m_LinSysRelativeTolInNonlin;
    }
    else
    {
        NekDouble tmpForc =
            m_forcingGamma * pow((resnorm / resnormOld), m_forcingAlpha);
        return max(min(m_LinSysRelativeTolInNonlin, tmpForc), 1.0E-6);
    }
}

} // namespace Nektar::LibUtilities
