///////////////////////////////////////////////////////////////////////////////
//
// File: NekNonlinSys.cpp
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
// Description: NekNonlinSys definition
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/LinearAlgebra/NekNonlinSys.h>

namespace Nektar::LibUtilities
{
/**
 * @class  NekNonlinSys
 *
 * Solves a nonlinear system using iterative methods.
 */
NekNonlinSysFactory &GetNekNonlinSysFactory()
{
    static NekNonlinSysFactory instance;
    return instance;
}

NekNonlinSys::NekNonlinSys(const LibUtilities::SessionReaderSharedPtr &pSession,
                           const LibUtilities::CommSharedPtr &vRowComm,
                           const int nDimen, const NekSysKey &pKey)
    : NekSys(pSession, vRowComm, nDimen, pKey)
{
    m_maxiter                   = pKey.m_NekNonlinSysMaxIterations;
    m_NonlinIterTolRelativeL2   = pKey.m_NonlinIterTolRelativeL2;
    m_LinSysRelativeTolInNonlin = pKey.m_LinSysRelativeTolInNonlin;
    m_LinSysIterSolverType      = pKey.m_LinSysIterSolverTypeInNonlin;

    ASSERTL0(LibUtilities::GetNekLinSysIterFactory().ModuleExists(
                 m_LinSysIterSolverType),
             "NekLinSysIter '" + m_LinSysIterSolverType +
                 "' is not defined.\n");

    m_linsol = LibUtilities::GetNekLinSysIterFactory().CreateInstance(
        m_LinSysIterSolverType, pSession, m_rowComm, m_SysDimen, pKey);
    m_linsol->SetFlagWarnings(false);
}

void NekNonlinSys::v_InitObject()
{
    NekSys::v_InitObject();
}

} // namespace Nektar::LibUtilities
