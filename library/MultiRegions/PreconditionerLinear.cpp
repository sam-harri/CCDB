///////////////////////////////////////////////////////////////////////////////
//
// File: PreconditionerLinear.cpp
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
// Description: PreconditionerLinear definition
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/GlobalLinSys.h>
#include <MultiRegions/GlobalLinSysIterativeStaticCond.h>
#include <MultiRegions/GlobalLinSysXxtFull.h>
#include <MultiRegions/GlobalMatrixKey.h>
#include <MultiRegions/PreconditionerLinear.h>

#ifdef NEKTAR_USING_PETSC
#include <MultiRegions/GlobalLinSysPETScFull.h>
#endif

#include <LocalRegions/MatrixKey.h>
#include <cmath>

using namespace std;

namespace Nektar::MultiRegions
{
/**
 * Registers the class with the Factory.
 */

string PreconditionerLinear::className1 =
    GetPreconFactory().RegisterCreatorFunction(
        "FullLinearSpace", PreconditionerLinear::create,
        "Full Linear space inverse Preconditioning");

std::string PreconditionerLinear::solveType =
    LibUtilities::SessionReader::RegisterDefaultSolverInfo("LinearPreconSolver",
                                                           "Xxt");
std::string PreconditionerLinear::solveTypeIds[] = {
    LibUtilities::SessionReader::RegisterEnumValue(
        "LinearPreconSolver", "PETSc", MultiRegions::eLinearPreconPETSc),
    LibUtilities::SessionReader::RegisterEnumValue(
        "LinearPreconSolver", "Xxt", MultiRegions::eLinearPreconXxt)};

/**
 * @class PreconditionerLinear
 *
 * This class implements preconditioning for the conjugate
 * gradient matrix solver.
 */

PreconditionerLinear::PreconditionerLinear(
    const std::shared_ptr<GlobalLinSys> &plinsys,
    const AssemblyMapSharedPtr &pLocToGloMap)
    : Preconditioner(plinsys, pLocToGloMap)
{
}

void PreconditionerLinear::v_InitObject()
{
}

void PreconditionerLinear::v_BuildPreconditioner()
{
    GlobalSysSolnType sType = m_locToGloMap.lock()->GetGlobalSysSolnType();
    ASSERTL0(sType == eIterativeStaticCond || sType == ePETScStaticCond,
             "This type of preconditioning is not implemented "
             "for this solver");

    std::shared_ptr<MultiRegions::ExpList> expList =
        ((m_linsys.lock())->GetLocMat()).lock();

    LinearPreconSolver solveType =
        expList->GetSession()->GetSolverInfoAsEnum<LinearPreconSolver>(
            "LinearPreconSolver");

    GlobalSysSolnType linSolveType;

    switch (solveType)
    {
        case eLinearPreconPETSc:
        {
            linSolveType = ePETScFullMatrix;
#ifndef NEKTAR_USING_PETSC
            NEKERROR(ErrorUtil::efatal, "Nektar++ has not been compiled with "
                                        "PETSc support.");
#endif
            break;
        }
        case eLinearPreconXxt:
        default:
        {
            linSolveType = eXxtFullMatrix;
            break;
        }
    }

    m_vertLocToGloMap =
        m_locToGloMap.lock()->LinearSpaceMap(*expList, linSolveType);

    SetupInvMult(m_vertLocToGloMap);

    // Generate linear solve system.
    StdRegions::MatrixType mType =
        m_linsys.lock()->GetKey().GetMatrixType() == StdRegions::eMass
            ? StdRegions::ePreconLinearSpaceMass
            : StdRegions::ePreconLinearSpace;

    GlobalLinSysKey preconKey(mType, m_vertLocToGloMap,
                              (m_linsys.lock())->GetKey().GetConstFactors());

    switch (solveType)
    {
        case eLinearPreconXxt:
        {
            m_vertLinsys =
                MemoryManager<GlobalLinSysXxtFull>::AllocateSharedPtr(
                    preconKey, expList, m_vertLocToGloMap);
            break;
        }
        case eLinearPreconPETSc:
        {
#ifdef NEKTAR_USING_PETSC
            m_vertLinsys =
                MemoryManager<GlobalLinSysPETScFull>::AllocateSharedPtr(
                    preconKey, expList, m_vertLocToGloMap);
#else
            ASSERTL0(false, "Nektar++ has not been compiled with "
                            "PETSc support.");
#endif
        }
    }
}

/**
 *
 */
void PreconditionerLinear::v_DoPreconditioner(
    const Array<OneD, NekDouble> &pInput, Array<OneD, NekDouble> &pOutput,
    const bool &isLocal)
{
    ASSERTL0(isLocal == false, "PreconditionerLinear is only set up for Global"
                               " iteratives sovles");
    v_DoPreconditionerWithNonVertOutput(pInput, pOutput, NullNekDouble1DArray,
                                        NullNekDouble1DArray);
}

/**
 *
 */
void PreconditionerLinear::v_DoPreconditionerWithNonVertOutput(
    const Array<OneD, NekDouble> &pInput, Array<OneD, NekDouble> &pOutput,
    const Array<OneD, NekDouble> &pNonVertOutput,
    Array<OneD, NekDouble> &pVertForce)
{
    GlobalSysSolnType solvertype = m_locToGloMap.lock()->GetGlobalSysSolnType();
    switch (solvertype)
    {
        case MultiRegions::eIterativeStaticCond:
        case MultiRegions::ePETScStaticCond:
        {
            int i, val;
            int nloc = m_vertLocToGloMap->GetNumLocalCoeffs();
            // mapping from full space to vertices
            Array<OneD, int> LocToGloBnd =
                m_vertLocToGloMap->GetLocalToGlobalBndMap();

            // Global to local for linear solver (different from above)
            Array<OneD, int> LocToGlo =
                m_vertLocToGloMap->GetLocalToGlobalMap();

            // number of Dir coeffs in from full problem
            int nDirFull = m_locToGloMap.lock()->GetNumGlobalDirBndCoeffs();

#if 0 
            int nglo = m_vertLocToGloMap->GetNumGlobalCoeffs();
            Array<OneD, NekDouble> In(nglo, 0.0);
            Array<OneD, NekDouble> Out(nglo, 0.0);

            // Gather rhs to local linear space.
            for (i = 0; i < nloc; ++i)
            {
                val = LocToGloBnd[i];
                if (val >= nDirFull)
                {
                    In[LocToGlo[i]] = pInput[val - nDirFull];
                }
            }
#else
            Array<OneD, NekDouble> In(nloc, 0.0);
            Array<OneD, NekDouble> Out(nloc, 0.0);

            // Gather rhs to local linear space.
            for (i = 0; i < nloc; ++i)
            {
                val = LocToGloBnd[i];
                if (val >= nDirFull)
                {
                    In[i] = pInput[val - nDirFull] * m_invMult[i];
                }
            }
#endif

            // Do solve without enforcing any boundary conditions.
            m_vertLinsys->SolveLinearSystem(
                m_vertLocToGloMap->GetNumGlobalCoeffs(), In, Out,
                m_vertLocToGloMap,
                m_vertLocToGloMap->GetNumGlobalDirBndCoeffs());

            if (pNonVertOutput != NullNekDouble1DArray)
            {
                ASSERTL1(pNonVertOutput.size() >= pOutput.size(),
                         "Non Vert output is not of sufficient length");
                Vmath::Vcopy(pOutput.size(), pNonVertOutput, 1, pOutput, 1);
            }
            else
            {
                // Copy input to output as a unit preconditioner on
                // any other value
                Vmath::Vcopy(pInput.size(), pInput, 1, pOutput, 1);
            }

#if 0 
            if (pVertForce != NullNekDouble1DArray)
            {
                Vmath::Zero(pVertForce.size(), pVertForce, 1);
                // Scatter back soln from linear solve
                for (i = 0; i < nloc; ++i)
                {
                    val = LocToGloBnd[i];
                    if (val >= nDirFull)
                    {
                        pOutput[val - nDirFull] = Out[LocToGlo[i]];
                        // copy vertex forcing into this vector
                        pVertForce[val - nDirFull] = In[LocToGlo[i]];
                    }
                }
            }
            else
            {
                // Scatter back soln from linear solve
                for (i = 0; i < nloc; ++i)
                {
                    val = LocToGloBnd[i];
                    if (val >= nDirFull)
                    {
                        pOutput[val - nDirFull] = Out[LocToGlo[i]];
                    }
                }
            }
#else
            if (pVertForce != NullNekDouble1DArray)
            {
                Vmath::Zero(pVertForce.size(), pVertForce, 1);
                // Scatter back soln from linear solve
                for (i = 0; i < nloc; ++i)
                {
                    val = LocToGloBnd[i];
                    if (val >= nDirFull)
                    {
                        pOutput[val - nDirFull] = Out[i];
                        // copy vertex forcing into this vector
                        pVertForce[val - nDirFull] = In[i];
                    }
                }
            }
            else
            {
                // Scatter back soln from linear solve
                for (i = 0; i < nloc; ++i)
                {
                    val = LocToGloBnd[i];
                    if (val >= nDirFull)
                    {
                        pOutput[val - nDirFull] = Out[i];
                    }
                }
            }
#endif
        }
        break;
        default:
            ASSERTL0(0, "Unsupported solver type");
            break;
    }
}

/**
 * Create the inverse multiplicity map.
 * @param   locToGloMap Local to global mapping information.
 */
void PreconditionerLinear::SetupInvMult(
    const std::shared_ptr<AssemblyMap> &pLocToGloMap)
{
    const Array<OneD, const int> &vMap = pLocToGloMap->GetLocalToGlobalMap();
    unsigned int nGlo                  = pLocToGloMap->GetNumGlobalCoeffs();
    unsigned int nEntries              = pLocToGloMap->GetNumLocalCoeffs();
    unsigned int i;

    // Count the multiplicity of each global DOF on this process
    Array<OneD, NekDouble> vCounts(nGlo, 0.0);
    for (i = 0; i < nEntries; ++i)
    {
        vCounts[vMap[i]] += 1.0;
    }

    // Get universal multiplicity by globally assembling counts
    pLocToGloMap->UniversalAssemble(vCounts);

    // Construct a map of 1/multiplicity for use in XXT solve
    m_invMult = Array<OneD, NekDouble>(nEntries);
    for (i = 0; i < nEntries; ++i)
    {
        m_invMult[i] = 1.0 / vCounts[vMap[i]];
    }
}

} // namespace Nektar::MultiRegions
