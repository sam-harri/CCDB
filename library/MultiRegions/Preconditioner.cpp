///////////////////////////////////////////////////////////////////////////////
//
// File: Preconditioner.cpp
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
// Description: Preconditioner definition
//
///////////////////////////////////////////////////////////////////////////////

#include <LocalRegions/MatrixKey.h>
#include <MultiRegions/GlobalMatrixKey.h>
#include <MultiRegions/Preconditioner.h>
#include <cmath>

namespace Nektar::MultiRegions
{

// register default solver value as diagonal
std::string Preconditioner::def =
    LibUtilities::SessionReader::RegisterDefaultSolverInfo("Preconditioner",
                                                           "Diagonal");
/**
 * @class Preconditioner
 *
 * This class implements preconditioning for the conjugate
 * gradient matrix solver.
 */

Preconditioner::Preconditioner(const std::shared_ptr<GlobalLinSys> &plinsys,
                               const AssemblyMapSharedPtr &pLocToGloMap)
    : m_linsys(plinsys), m_preconType(pLocToGloMap->GetPreconType()),
      m_locToGloMap(pLocToGloMap)
{
}

/**
 *
 */
PreconFactory &GetPreconFactory()
{
    static PreconFactory instance;
    return instance;
}

void Preconditioner::v_InitObject()
{
    NEKERROR(ErrorUtil::efatal, "Method does not exist");
}

/**
 * \brief Apply a preconditioner to the conjugate gradient method
 */
void Preconditioner::v_DoPreconditioner(
    [[maybe_unused]] const Array<OneD, NekDouble> &pInput,
    [[maybe_unused]] Array<OneD, NekDouble> &pOutput,
    [[maybe_unused]] const bool &IsLocal)
{
    NEKERROR(ErrorUtil::efatal, "Method does not exist");
}

/**
 * \brief Apply an assembly and scatter back to lcoal array
 */
void Preconditioner::DoAssembleLoc(const Array<OneD, NekDouble> &pInput,
                                   Array<OneD, NekDouble> &pOutput,
                                   const bool &ZeroDir)
{
    auto assMap                  = m_locToGloMap.lock();
    GlobalSysSolnType solvertype = assMap->GetGlobalSysSolnType();

    if (solvertype == eIterativeFull)
    {
        assMap->Assemble(pInput, pOutput);
        if (ZeroDir)
        {
            int nDir = assMap->GetNumGlobalDirBndCoeffs();
            Vmath::Zero(nDir, pOutput, 1);
        }
        assMap->GlobalToLocal(pOutput, pOutput);
    }
    else // bnd version.
    {
        assMap->AssembleBnd(pInput, pOutput);
        if (ZeroDir)
        {
            int nDir = assMap->GetNumGlobalDirBndCoeffs();
            Vmath::Zero(nDir, pOutput, 1);
        }
        assMap->GlobalToLocalBnd(pOutput, pOutput);
    }
}

/**
 * \brief Apply a preconditioner to the conjugate gradient method with
 * an output for non-vertex degrees of freedom.
 */
void Preconditioner::v_DoPreconditionerWithNonVertOutput(
    [[maybe_unused]] const Array<OneD, NekDouble> &pInput,
    [[maybe_unused]] Array<OneD, NekDouble> &pOutput,
    [[maybe_unused]] const Array<OneD, NekDouble> &pNonVertOutput,
    [[maybe_unused]] Array<OneD, NekDouble> &pVertForce)
{
    NEKERROR(ErrorUtil::efatal, "Method does not exist");
}

/**
 * \brief Transform from original basis to low energy basis
 */
void Preconditioner::v_DoTransformBasisToLowEnergy(
    [[maybe_unused]] Array<OneD, NekDouble> &pInOut)
{
}

/**
 * \brief Transform from low energy coeffs to orignal basis
 */
void Preconditioner::v_DoTransformCoeffsFromLowEnergy(
    [[maybe_unused]] Array<OneD, NekDouble> &pInOut)
{
}

/**
 * \brief Multiply by the block inverse transformation matrix
 */
void Preconditioner::v_DoTransformCoeffsToLowEnergy(
    [[maybe_unused]] const Array<OneD, NekDouble> &pInput,
    [[maybe_unused]] Array<OneD, NekDouble> &pOutput)
{
}

/**
 * \brief Multiply by the block transposed inverse transformation matrix
 */
void Preconditioner::v_DoTransformBasisFromLowEnergy(
    [[maybe_unused]] const Array<OneD, NekDouble> &pInput,
    [[maybe_unused]] Array<OneD, NekDouble> &pOutput)
{
    NEKERROR(ErrorUtil::efatal, "Method does not exist");
}

void Preconditioner::v_BuildPreconditioner()
{
}

/**
 * \brief Get block elemental transposed transformation matrix
 * \f$\mathbf{R}^{T}\f$
 */
DNekScalMatSharedPtr Preconditioner::v_TransformedSchurCompl(
    [[maybe_unused]] int offset, [[maybe_unused]] int bnd_offset,
    const std::shared_ptr<DNekScalMat> &loc_mat)
{
    return loc_mat;
}

/**
 * @brief Performs global assembly of diagonal entries to global Schur
 * complement matrix.
 */
Array<OneD, NekDouble> Preconditioner::AssembleStaticCondGlobalDiagonals()
{
    auto asmMap = m_locToGloMap.lock();

    int nGlobalBnd = asmMap->GetNumGlobalBndCoeffs();
    int nDirBnd    = asmMap->GetNumGlobalDirBndCoeffs();
    int rows       = nGlobalBnd - nDirBnd;

    DNekScalBlkMatSharedPtr loc_mat;
    DNekScalMatSharedPtr bnd_mat;
    int sign1, sign2, gid1, gid2, i, j, n, cnt;
    Array<OneD, NekDouble> diagonals(rows, 0.0);

    // Extract diagonal contributions of globally assembled
    // schur complement matrix
    for (cnt = n = 0; n < m_linsys.lock()->GetNumBlocks(); ++n)
    {
        // Get statically condensed local matrix.
        loc_mat = (m_linsys.lock())->GetStaticCondBlock(n);

        // Extract boundary block.
        bnd_mat = loc_mat->GetBlock(0, 0);

        // Offset by number of rows.
        int bnd_row = bnd_mat->GetRows();

        for (i = 0; i < bnd_row; ++i)
        {
            gid1  = asmMap->GetLocalToGlobalBndMap(cnt + i) - nDirBnd;
            sign1 = asmMap->GetLocalToGlobalBndSign(cnt + i);

            if (gid1 < 0)
            {
                continue;
            }

            for (j = 0; j < bnd_row; ++j)
            {
                gid2  = asmMap->GetLocalToGlobalBndMap(cnt + j) - nDirBnd;
                sign2 = asmMap->GetLocalToGlobalBndSign(cnt + j);

                if (gid2 == gid1)
                {
                    diagonals[gid1] += sign1 * sign2 * (*bnd_mat)(i, j);
                }
            }
        }
        cnt += bnd_row;
    }

    return diagonals;
}
} // namespace Nektar::MultiRegions
