///////////////////////////////////////////////////////////////////////////////
//
// File: GlobalLinSysXxtFull.cpp
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
// Description: GlobalLinSysXxtFull definition
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Communication/Xxt.hpp>
#include <MultiRegions/ExpList.h>
#include <MultiRegions/GlobalLinSysXxtFull.h>

using namespace std;

namespace Nektar::MultiRegions
{
/**
 * @class GlobalLinSysXxtFull
 */

/**
 * Registers the class with the Factory.
 */
string GlobalLinSysXxtFull::className =
    GetGlobalLinSysFactory().RegisterCreatorFunction(
        "XxtFull", GlobalLinSysXxtFull::create, "Xxt Full Matrix.");

/// Constructor for full direct matrix solve.
GlobalLinSysXxtFull::GlobalLinSysXxtFull(
    const GlobalLinSysKey &pLinSysKey, const std::weak_ptr<ExpList> &pExp,
    const std::shared_ptr<AssemblyMap> &pLocToGloMap)
    : GlobalLinSys(pLinSysKey, pExp, pLocToGloMap),
      GlobalLinSysXxt(pLinSysKey, pExp, pLocToGloMap)
{

    ASSERTL1(m_linSysKey.GetGlobalSysSolnType() == eXxtFullMatrix,
             "This routine should only be used when using a Full XXT"
             " matrix solve");

    AssembleMatrixArrays(pLocToGloMap);
}

/**
 * Construct the local matrix row index, column index and value index
 * arrays and initialize the XXT data structure with this information.
 * @param   locToGloMap Local to global mapping information.
 */
void GlobalLinSysXxtFull::AssembleMatrixArrays(
    const std::shared_ptr<AssemblyMap> &pLocToGloMap)
{
    ExpListSharedPtr vExp = m_expList.lock();
    unsigned int nElmt    = vExp->GetNumElmts();
    DNekScalMatSharedPtr loc_mat;
    unsigned int iCount    = 0;
    unsigned int rCount    = 0;
    unsigned int nRows     = 0;
    unsigned int nEntries  = 0;
    unsigned int numDirBnd = pLocToGloMap->GetNumGlobalDirBndCoeffs();
    unsigned int nLocal    = pLocToGloMap->GetNumLocalCoeffs();
    const Array<OneD, NekDouble> &vMapSign =
        pLocToGloMap->GetLocalToGlobalSign();
    bool doSign    = pLocToGloMap->GetSignChange();
    unsigned int i = 0, j = 0, k = 0, n = 0;
    int gid1;
    Array<OneD, unsigned int> vSizes(nElmt);

    // First construct a map of the number of local DOFs in each block
    // and the number of matrix entries for each block

    // Dimension of matrix is just the linear vertex space
    if ((m_linSysKey.GetMatrixType() == StdRegions::ePreconLinearSpace) ||
        (m_linSysKey.GetMatrixType() == StdRegions::ePreconLinearSpaceMass))
    {
        for (n = 0; n < nElmt; ++n)
        {
            vSizes[n] = vExp->GetExp(n)->GetNverts();
            nEntries += vSizes[n] * vSizes[n];
        }
    }
    else
    {
        for (n = 0; n < nElmt; ++n)
        {
            vSizes[n] = vExp->GetExp(n)->GetNcoeffs();
            nEntries += vSizes[n] * vSizes[n];
        }
    }

    // Set up i-index, j-index and value arrays
    m_Ai = Array<OneD, unsigned int>(nEntries);
    m_Aj = Array<OneD, unsigned int>(nEntries);
    m_Ar = Array<OneD, double>(nEntries, 0.0);

    // Set up the universal ID array for XXT
    Array<OneD, unsigned long> vId(nLocal);

    // Loop over each elemental block, extract matrix indices and value
    // and set the universal ID array
    for (n = iCount = 0; n < nElmt; ++n)
    {
        loc_mat = GetBlock(n);
        nRows   = loc_mat->GetRows();

        for (i = 0; i < nRows; ++i)
        {
            gid1 = pLocToGloMap->GetLocalToGlobalMap(iCount + i);
            for (j = 0; j < nRows; ++j)
            {
                k       = rCount + i * vSizes[n] + j;
                m_Ai[k] = iCount + i;
                m_Aj[k] = iCount + j;
                m_Ar[k] = (*loc_mat)(i, j);
                if (doSign)
                {
                    m_Ar[k] *= vMapSign[iCount + i] * vMapSign[iCount + j];
                }
            }

            // Dirichlet DOFs are not included in the solve, so we set
            // these to the special XXT id=0.
            if (gid1 < numDirBnd)
            {
                vId[iCount + i] = 0;
            }
            else
            {
                vId[iCount + i] = pLocToGloMap->GetGlobalToUniversalMap(gid1);
            }
        }
        iCount += vSizes[n];
        rCount += vSizes[n] * vSizes[n];
    }

    // Set up XXT and output some stats
    LibUtilities::CommSharedPtr vComm = pLocToGloMap->GetComm();
    m_crsData = Xxt::Init(nLocal, vId, m_Ai, m_Aj, m_Ar, vComm);
    if (m_verbose)
    {
        Xxt::nektar_crs_stats(m_crsData);
    }
}

/**
 * Solve the linear system using a full global matrix system.
 */
void GlobalLinSysXxtFull::v_Solve(
    const Array<OneD, const NekDouble> &pLocInput,
    Array<OneD, NekDouble> &pLocOutput,
    const AssemblyMapSharedPtr &pLocToGloMap,
    const Array<OneD, const NekDouble> &pDirForcing)
{
    bool dirForcCalculated = (bool)pDirForcing.size();
    int nDirDofs           = pLocToGloMap->GetNumGlobalDirBndCoeffs();
    int nLocDofs           = pLocToGloMap->GetNumLocalCoeffs();

    if (nDirDofs)
    {
        std::shared_ptr<MultiRegions::ExpList> expList = m_expList.lock();
        Array<OneD, NekDouble> rhs(nLocDofs);

        // Calculate the Dirichlet forcing
        if (dirForcCalculated)
        {
            // Assume pDirForcing is in local space
            ASSERTL0(
                pDirForcing.size() >= nLocDofs,
                "DirForcing is not of sufficient size. Is it in local space?");
            Vmath::Vsub(nLocDofs, pLocInput, 1, pDirForcing, 1, rhs, 1);
        }
        else
        {
            // Calculate initial condition and Dirichlet forcing and subtract it
            // from the rhs
            expList->GeneralMatrixOp(m_linSysKey, pLocOutput, rhs);

            // Iterate over all the elements computing Robin BCs where
            // necessary
            for (auto &r : m_robinBCInfo) // add robin mass matrix
            {
                RobinBCInfoSharedPtr rBC;
                Array<OneD, NekDouble> rhsloc;

                int n      = r.first;
                int offset = expList->GetCoeff_Offset(n);

                LocalRegions::ExpansionSharedPtr vExp = expList->GetExp(n);
                // Add local matrix contribution
                for (rBC = r.second; rBC; rBC = rBC->next)
                {
                    vExp->AddRobinTraceContribution(
                        rBC->m_robinID, rBC->m_robinPrimitiveCoeffs,
                        pLocOutput + offset, rhsloc = rhs + offset);
                }
            }
            Vmath::Vsub(nLocDofs, pLocInput, 1, rhs, 1, rhs, 1);
        }

        Array<OneD, NekDouble> diff(nLocDofs);

        // Solve for perturbation from initial guess in pOutput
        SolveLinearSystem(nLocDofs, rhs, diff, pLocToGloMap);

        // Add back initial and boundary condition
        Vmath::Vadd(nLocDofs, diff, 1, pLocOutput, 1, pLocOutput, 1);
    }
    else
    {
        SolveLinearSystem(nLocDofs, pLocInput, pLocOutput, pLocToGloMap);
    }
}

/// Solve the linear system for given input and output vectors.
void GlobalLinSysXxtFull::v_SolveLinearSystem(
    [[maybe_unused]] const int pNumRows,
    const Array<OneD, const NekDouble> &pInput, Array<OneD, NekDouble> &pOutput,
    const AssemblyMapSharedPtr &pLocToGloMap,
    [[maybe_unused]] const int pNumDir)
{
    int nLocal = pNumRows;

    Vmath::Zero(nLocal, pOutput, 1);

    // Set Output into correct sign
    if (pLocToGloMap->GetSignChange())
    {
        Array<OneD, NekDouble> vlocal(nLocal);
        Vmath::Vmul(nLocal, pLocToGloMap->GetLocalToGlobalSign(), 1, pInput, 1,
                    vlocal, 1);

        Xxt::Solve(pOutput, m_crsData, vlocal);

        Vmath::Vmul(nLocal, pLocToGloMap->GetLocalToGlobalSign(), 1, pOutput, 1,
                    pOutput, 1);
    }
    else
    {
        Xxt::Solve(pOutput, m_crsData, pInput);
    }
}

} // namespace Nektar::MultiRegions
