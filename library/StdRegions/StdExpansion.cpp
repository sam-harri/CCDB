///////////////////////////////////////////////////////////////////////////////
//
// File: StdExpansion.cpp
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
// Description: Definition of methods in class StdExpansion which is
// the base class to all expansion shapes
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Foundations/ManagerAccess.h> // for BasisManager, etc
#include <StdRegions/StdExpansion.h>

namespace Nektar::StdRegions
{
/** \brief Default constructor */
StdExpansion::StdExpansion()
{
}

StdExpansion::StdExpansion(const int numcoeffs, const int numbases,
                           const LibUtilities::BasisKey &Ba,
                           const LibUtilities::BasisKey &Bb,
                           const LibUtilities::BasisKey &Bc)
    : m_base(numbases), m_elmt_id(0), m_ncoeffs(numcoeffs),
      m_stdMatrixManager(std::bind(&StdExpansion::CreateStdMatrix, this,
                                   std::placeholders::_1),
                         std::string("StdExpansionStdMatrix")),
      m_stdStaticCondMatrixManager(
          std::bind(&StdExpansion::CreateStdStaticCondMatrix, this,
                    std::placeholders::_1),
          std::string("StdExpansionStdStaticCondMatrix"))
{
    switch (m_base.size())
    {
        case 3:
            ASSERTL2(Bc != LibUtilities::NullBasisKey,
                     "NULL Basis attempting to be used.");
            m_base[2] = LibUtilities::BasisManager()[Bc];
            /* Falls through. */
        case 2:
            ASSERTL2(Bb != LibUtilities::NullBasisKey,
                     "NULL Basis attempting to be used.");
            m_base[1] = LibUtilities::BasisManager()[Bb];
            /* Falls through. */
        case 1:
            ASSERTL2(Ba != LibUtilities::NullBasisKey,
                     "NULL Basis attempting to be used.");
            m_base[0] = LibUtilities::BasisManager()[Ba];
            break;
        default:
            break;
            // ASSERTL0(false, "numbases incorrectly specified");
    };

} // end constructor

StdExpansion::StdExpansion(const StdExpansion &T)
    : std::enable_shared_from_this<StdExpansion>(T), m_base(T.m_base),
      m_elmt_id(T.m_elmt_id), m_ncoeffs(T.m_ncoeffs),
      m_stdMatrixManager(T.m_stdMatrixManager),
      m_stdStaticCondMatrixManager(T.m_stdStaticCondMatrixManager)
{
}

// Destructor
StdExpansion::~StdExpansion()
{
}

NekDouble StdExpansion::Linf(const Array<OneD, const NekDouble> &phys,
                             const Array<OneD, const NekDouble> &sol)
{
    NekDouble val;
    int ntot = GetTotPoints();
    Array<OneD, NekDouble> wsp(ntot);

    if (sol == NullNekDouble1DArray)
    {
        Vmath::Vabs(ntot, phys, 1, wsp, 1);
    }
    else
    {
        Vmath::Vsub(ntot, sol, 1, phys, 1, wsp, 1);
        Vmath::Vabs(ntot, wsp, 1, wsp, 1);
    }

    val = Vmath::Vamax(ntot, wsp, 1);
    return val;
}

NekDouble StdExpansion::L2(const Array<OneD, const NekDouble> &phys,
                           const Array<OneD, const NekDouble> &sol)
{
    NekDouble val;
    int ntot = GetTotPoints();
    Array<OneD, NekDouble> wsp(ntot);

    if (sol.size() == 0)
    {
        Vmath::Vmul(ntot, phys, 1, phys, 1, wsp, 1);
    }
    else
    {
        Vmath::Vsub(ntot, sol, 1, phys, 1, wsp, 1);
        Vmath::Vmul(ntot, wsp, 1, wsp, 1, wsp, 1);
    }

    val = v_Integral(wsp);

    return (val < 0.0) ? 0.0 : sqrt(val);
}

NekDouble StdExpansion::H1(const Array<OneD, const NekDouble> &phys,
                           const Array<OneD, const NekDouble> &sol)
{
    int i;
    NekDouble val;
    int ntot    = GetTotPoints();
    int coordim = v_GetCoordim();
    Array<OneD, NekDouble> wsp(3 * ntot);
    Array<OneD, NekDouble> wsp_deriv = wsp + ntot;
    Array<OneD, NekDouble> sum       = wsp_deriv + ntot;

    if (sol == NullNekDouble1DArray)
    {
        Vmath::Vcopy(ntot, phys, 1, wsp, 1);
        Vmath::Vmul(ntot, phys, 1, phys, 1, sum, 1);
    }
    else
    {
        Vmath::Vsub(ntot, sol, 1, phys, 1, wsp, 1);
        Vmath::Vmul(ntot, wsp, 1, wsp, 1, sum, 1);
    }

    for (i = 0; i < coordim; ++i)
    {
        v_PhysDeriv(i, wsp, wsp_deriv);
        Vmath::Vvtvp(ntot, wsp_deriv, 1, wsp_deriv, 1, sum, 1, sum, 1);
    }

    val = sqrt(v_Integral(sum));

    return val;
}

DNekBlkMatSharedPtr StdExpansion::CreateStdStaticCondMatrix(
    const StdMatrixKey &mkey)
{
    DNekBlkMatSharedPtr returnval;

    DNekMatSharedPtr mat = GetStdMatrix(mkey);
    int nbdry = NumBndryCoeffs(); // also checks to see if this is a boundary
                                  // interior decomposed expansion
    int nint = m_ncoeffs - nbdry;
    DNekMatSharedPtr A =
        MemoryManager<DNekMat>::AllocateSharedPtr(nbdry, nbdry);
    DNekMatSharedPtr B = MemoryManager<DNekMat>::AllocateSharedPtr(nbdry, nint);
    DNekMatSharedPtr C = MemoryManager<DNekMat>::AllocateSharedPtr(nint, nbdry);
    DNekMatSharedPtr D = MemoryManager<DNekMat>::AllocateSharedPtr(nint, nint);

    int i, j;

    Array<OneD, unsigned int> bmap(nbdry);
    Array<OneD, unsigned int> imap(nint);
    GetBoundaryMap(bmap);
    GetInteriorMap(imap);

    for (i = 0; i < nbdry; ++i)
    {
        for (j = 0; j < nbdry; ++j)
        {
            (*A)(i, j) = (*mat)(bmap[i], bmap[j]);
        }

        for (j = 0; j < nint; ++j)
        {
            (*B)(i, j) = (*mat)(bmap[i], imap[j]);
        }
    }

    for (i = 0; i < nint; ++i)
    {
        for (j = 0; j < nbdry; ++j)
        {
            (*C)(i, j) = (*mat)(imap[i], bmap[j]);
        }

        for (j = 0; j < nint; ++j)
        {
            (*D)(i, j) = (*mat)(imap[i], imap[j]);
        }
    }

    // Calculate static condensed system
    if (nint)
    {
        D->Invert();
        (*B) = (*B) * (*D);
        (*A) = (*A) - (*B) * (*C);
    }

    // set up block matrix system
    Array<OneD, unsigned int> exp_size(2);
    exp_size[0] = nbdry;
    exp_size[1] = nint;
    returnval =
        MemoryManager<DNekBlkMat>::AllocateSharedPtr(exp_size, exp_size);

    returnval->SetBlock(0, 0, A);
    returnval->SetBlock(0, 1, B);
    returnval->SetBlock(1, 0, C);
    returnval->SetBlock(1, 1, D);

    return returnval;
}

DNekMatSharedPtr StdExpansion::CreateGeneralMatrix(const StdMatrixKey &mkey)
{
    int i;
    DNekMatSharedPtr returnval;

    switch (mkey.GetMatrixType())
    {
        case eInvMass:
        {
            StdMatrixKey masskey(eMass, mkey.GetShapeType(), *this,
                                 NullConstFactorMap, NullVarCoeffMap,
                                 mkey.GetNodalPointsType());
            DNekMatSharedPtr mmat = GetStdMatrix(masskey);

            returnval = MemoryManager<DNekMat>::AllocateSharedPtr(
                *mmat); // Populate standard mass matrix.
            returnval->Invert();
        }
        break;
        case eInvNBasisTrans:
        {
            StdMatrixKey tmpkey(eNBasisTrans, mkey.GetShapeType(), *this,
                                NullConstFactorMap, NullVarCoeffMap,
                                mkey.GetNodalPointsType());
            DNekMatSharedPtr tmpmat = GetStdMatrix(tmpkey);
            returnval               = MemoryManager<DNekMat>::AllocateSharedPtr(
                *tmpmat); // Populate  matrix.
            returnval->Invert();
        }
        break;
        case eBwdMat:
        {
            int nq = GetTotPoints();
            Array<OneD, NekDouble> tmpin(m_ncoeffs);
            Array<OneD, NekDouble> tmpout(nq);

            returnval =
                MemoryManager<DNekMat>::AllocateSharedPtr(m_ncoeffs, nq);
            Array<OneD, NekDouble> Bwd_data = returnval->GetPtr();

            StdRegions::StdMatrixKey matkey(StdRegions::eBwdTrans,
                                            this->DetShapeType(), *this);
            DNekMatSharedPtr MatBwdTrans         = GetStdMatrix(matkey);
            Array<OneD, NekDouble> BwdTrans_data = MatBwdTrans->GetPtr();

            for (i = 0; i < m_ncoeffs; ++i)
            {
                Array<OneD, NekDouble> tmpinn = BwdTrans_data + nq * i;
                tmpout                        = Bwd_data + i;

                Vmath::Vcopy(nq, tmpinn, 1, tmpout, m_ncoeffs);
            }
        }
        break;
        case eBwdTrans:
        {
            int nq = GetTotPoints();
            Array<OneD, NekDouble> tmpin(m_ncoeffs);
            Array<OneD, NekDouble> tmpout(nq);

            returnval =
                MemoryManager<DNekMat>::AllocateSharedPtr(nq, m_ncoeffs);

            for (i = 0; i < m_ncoeffs; ++i)
            {
                Vmath::Zero(m_ncoeffs, tmpin, 1);
                tmpin[i] = 1.0;

                BwdTrans_SumFac(tmpin, tmpout);

                Vmath::Vcopy(nq, tmpout.data(), 1,
                             returnval->GetRawPtr() + i * nq, 1);
            }
        }
        break;
        case eIProductWRTBase:
        {
            int nq = GetTotPoints();
            Array<OneD, NekDouble> tmpin(nq);
            Array<OneD, NekDouble> tmpout(m_ncoeffs);

            returnval =
                MemoryManager<DNekMat>::AllocateSharedPtr(m_ncoeffs, nq);

            for (i = 0; i < nq; ++i)
            {
                Vmath::Zero(nq, tmpin, 1);
                tmpin[i] = 1.0;

                IProductWRTBase_SumFac(tmpin, tmpout);

                Vmath::Vcopy(m_ncoeffs, tmpout.data(), 1,
                             returnval->GetRawPtr() + i * m_ncoeffs, 1);
            }
        }
        break;
        case eIProductWRTDerivBase0:
        {
            int nq = GetTotPoints();
            Array<OneD, NekDouble> tmpin(nq);
            Array<OneD, NekDouble> tmpout(m_ncoeffs);

            returnval =
                MemoryManager<DNekMat>::AllocateSharedPtr(nq, m_ncoeffs);

            for (i = 0; i < nq; ++i)
            {
                Vmath::Zero(nq, tmpin, 1);
                tmpin[i] = 1.0;

                IProductWRTDerivBase_SumFac(0, tmpin, tmpout);

                Vmath::Vcopy(m_ncoeffs, tmpout.data(), 1,
                             returnval->GetRawPtr() + i * m_ncoeffs, 1);
            }
        }
        break;
        case eIProductWRTDerivBase1:
        {
            int nq = GetTotPoints();
            Array<OneD, NekDouble> tmpin(nq);
            Array<OneD, NekDouble> tmpout(m_ncoeffs);

            returnval =
                MemoryManager<DNekMat>::AllocateSharedPtr(nq, m_ncoeffs);

            for (i = 0; i < nq; ++i)
            {
                Vmath::Zero(nq, tmpin, 1);
                tmpin[i] = 1.0;

                IProductWRTDerivBase_SumFac(1, tmpin, tmpout);

                Vmath::Vcopy(m_ncoeffs, tmpout.data(), 1,
                             returnval->GetRawPtr() + i * m_ncoeffs, 1);
            }
        }
        break;
        case eIProductWRTDerivBase2:
        {
            int nq = GetTotPoints();
            Array<OneD, NekDouble> tmpin(nq);
            Array<OneD, NekDouble> tmpout(m_ncoeffs);

            returnval =
                MemoryManager<DNekMat>::AllocateSharedPtr(nq, m_ncoeffs);

            for (i = 0; i < nq; ++i)
            {
                Vmath::Zero(nq, tmpin, 1);
                tmpin[i] = 1.0;

                IProductWRTDerivBase_SumFac(2, tmpin, tmpout);

                Vmath::Vcopy(m_ncoeffs, tmpout.data(), 1,
                             returnval->GetRawPtr() + i * m_ncoeffs, 1);
            }
        }
        break;
        case eDerivBase0:
        {
            int nq = GetTotPoints();
            returnval =
                MemoryManager<DNekMat>::AllocateSharedPtr(m_ncoeffs, nq);
            GenStdMatBwdDeriv(0, returnval);
        }
        break;
        case eDerivBase1:
        {
            int nq = GetTotPoints();
            returnval =
                MemoryManager<DNekMat>::AllocateSharedPtr(m_ncoeffs, nq);
            GenStdMatBwdDeriv(1, returnval);
        }
        break;
        case eDerivBase2:
        {
            int nq = GetTotPoints();
            returnval =
                MemoryManager<DNekMat>::AllocateSharedPtr(m_ncoeffs, nq);
            GenStdMatBwdDeriv(2, returnval);
        }
        break;
        case eEquiSpacedToCoeffs:
        {
            // check to see if equispaced basis
            int nummodes    = m_base[0]->GetNumModes();
            bool equispaced = true;
            for (i = 1; i < m_base.size(); ++i)
            {
                if (m_base[i]->GetNumModes() != nummodes)
                {
                    equispaced = false;
                }
            }

            ASSERTL0(equispaced,
                     "Currently need to have same num modes in all "
                     "directionmodes to use EquiSpacedToCoeff method");

            int ntot = GetTotPoints();
            Array<OneD, NekDouble> qmode(ntot);
            Array<OneD, NekDouble> emode(m_ncoeffs);

            returnval =
                MemoryManager<DNekMat>::AllocateSharedPtr(m_ncoeffs, m_ncoeffs);
            for (i = 0; i < m_ncoeffs; ++i)
            {
                // Get mode at quadrature points
                FillMode(i, qmode);

                // interpolate to equi spaced
                PhysInterpToSimplexEquiSpaced(qmode, emode, nummodes);

                // fill matrix
                Vmath::Vcopy(m_ncoeffs, &emode[0], 1,
                             returnval->GetRawPtr() + i * m_ncoeffs, 1);
            }
            // invert matrix
            returnval->Invert();
        }
        break;
        case eMass:
        case eHelmholtz:
        case eLaplacian:
        case eLaplacian00:
        case eLaplacian01:
        case eLaplacian02:
        case eLaplacian10:
        case eLaplacian11:
        case eLaplacian12:
        case eLaplacian20:
        case eLaplacian21:
        case eLaplacian22:
        case eWeakDeriv0:
        case eWeakDeriv1:
        case eWeakDeriv2:
        case eWeakDirectionalDeriv:
        case eMassLevelCurvature:
        case eLinearAdvection:
        case eLinearAdvectionReaction:
        case eLinearAdvectionDiffusionReaction:
        {
            Array<OneD, NekDouble> tmp(m_ncoeffs);
            returnval =
                MemoryManager<DNekMat>::AllocateSharedPtr(m_ncoeffs, m_ncoeffs);
            DNekMat &Mat = *returnval;

            for (i = 0; i < m_ncoeffs; ++i)
            {
                Vmath::Zero(m_ncoeffs, tmp, 1);
                tmp[i] = 1.0;

                GeneralMatrixOp_MatFree(tmp, tmp, mkey);

                Vmath::Vcopy(m_ncoeffs, &tmp[0], 1,
                             &(Mat.GetPtr())[0] + i * m_ncoeffs, 1);
            }
        }
        break;
        default:
        {
            NEKERROR(ErrorUtil::efatal,
                     "This type of matrix, " +
                         static_cast<std::string>(
                             MatrixTypeMap[mkey.GetMatrixType()]) +
                         ", can not be created using a general approach");
        }
        break;
    }

    return returnval;
}

void StdExpansion::GeneralMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                   Array<OneD, NekDouble> &outarray,
                                   const StdMatrixKey &mkey)
{
    switch (mkey.GetMatrixType())
    {
        case eMass:
            MassMatrixOp(inarray, outarray, mkey);
            break;
        case eWeakDeriv0:
            WeakDerivMatrixOp(0, inarray, outarray, mkey);
            break;
        case eWeakDeriv1:
            WeakDerivMatrixOp(1, inarray, outarray, mkey);
            break;
        case eWeakDeriv2:
            WeakDerivMatrixOp(2, inarray, outarray, mkey);
            break;
        case eWeakDirectionalDeriv:
            WeakDirectionalDerivMatrixOp(inarray, outarray, mkey);
            break;
        case eMassLevelCurvature:
            MassLevelCurvatureMatrixOp(inarray, outarray, mkey);
            break;
        case eLinearAdvection:
            LinearAdvectionMatrixOp(inarray, outarray, mkey);
            break;
        case eLinearAdvectionReaction:
            LinearAdvectionDiffusionReactionMatrixOp(inarray, outarray, mkey,
                                                     false);
            break;
        case eLinearAdvectionDiffusionReaction:
            LinearAdvectionDiffusionReactionMatrixOp(inarray, outarray, mkey);
            break;
        case eLaplacian:
            LaplacianMatrixOp(inarray, outarray, mkey);
            break;
        case eLaplacian00:
            LaplacianMatrixOp(0, 0, inarray, outarray, mkey);
            break;
        case eLaplacian01:
            LaplacianMatrixOp(0, 1, inarray, outarray, mkey);
            break;
        case eLaplacian02:
            LaplacianMatrixOp(0, 2, inarray, outarray, mkey);
            break;
        case eLaplacian10:
            LaplacianMatrixOp(1, 0, inarray, outarray, mkey);
            break;
        case eLaplacian11:
            LaplacianMatrixOp(1, 1, inarray, outarray, mkey);
            break;
        case eLaplacian12:
            LaplacianMatrixOp(1, 2, inarray, outarray, mkey);
            break;
        case eLaplacian20:
            LaplacianMatrixOp(2, 0, inarray, outarray, mkey);
            break;
        case eLaplacian21:
            LaplacianMatrixOp(2, 1, inarray, outarray, mkey);
            break;
        case eLaplacian22:
            LaplacianMatrixOp(2, 2, inarray, outarray, mkey);
            break;
        case eHelmholtz:
            HelmholtzMatrixOp(inarray, outarray, mkey);
            break;
        default:
            NEKERROR(ErrorUtil::efatal,
                     "This matrix does not have an operator");
            break;
    }
}

void StdExpansion::GeneralMatrixOp_MatFree(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const StdMatrixKey &mkey)
{
    switch (mkey.GetMatrixType())
    {
        case eMass:
            MassMatrixOp_MatFree(inarray, outarray, mkey);
            break;
        case eWeakDeriv0:
            WeakDerivMatrixOp_MatFree(0, inarray, outarray, mkey);
            break;
        case eWeakDeriv1:
            WeakDerivMatrixOp_MatFree(1, inarray, outarray, mkey);
            break;
        case eWeakDeriv2:
            WeakDerivMatrixOp_MatFree(2, inarray, outarray, mkey);
            break;
        case eWeakDirectionalDeriv:
            WeakDirectionalDerivMatrixOp_MatFree(inarray, outarray, mkey);
            break;
        case eMassLevelCurvature:
            MassLevelCurvatureMatrixOp_MatFree(inarray, outarray, mkey);
            break;
        case eLinearAdvection:
            LinearAdvectionMatrixOp_MatFree(inarray, outarray, mkey);
            break;
        case eLinearAdvectionReaction:
            LinearAdvectionDiffusionReactionMatrixOp_MatFree(inarray, outarray,
                                                             mkey, false);
            break;
        case eLinearAdvectionDiffusionReaction:
            LinearAdvectionDiffusionReactionMatrixOp_MatFree(inarray, outarray,
                                                             mkey);
            break;
        case eLaplacian:
            LaplacianMatrixOp_MatFree(inarray, outarray, mkey);
            break;
        case eLaplacian00:
            LaplacianMatrixOp_MatFree(0, 0, inarray, outarray, mkey);
            break;
        case eLaplacian01:
            LaplacianMatrixOp_MatFree(0, 1, inarray, outarray, mkey);
            break;
        case eLaplacian02:
            LaplacianMatrixOp_MatFree(0, 2, inarray, outarray, mkey);
            break;
        case eLaplacian10:
            LaplacianMatrixOp_MatFree(1, 0, inarray, outarray, mkey);
            break;
        case eLaplacian11:
            LaplacianMatrixOp_MatFree(1, 1, inarray, outarray, mkey);
            break;
        case eLaplacian12:
            LaplacianMatrixOp_MatFree(1, 2, inarray, outarray, mkey);
            break;
        case eLaplacian20:
            LaplacianMatrixOp_MatFree(2, 0, inarray, outarray, mkey);
            break;
        case eLaplacian21:
            LaplacianMatrixOp_MatFree(2, 1, inarray, outarray, mkey);
            break;
        case eLaplacian22:
            LaplacianMatrixOp_MatFree(2, 2, inarray, outarray, mkey);
            break;
        case eHelmholtz:
            HelmholtzMatrixOp_MatFree(inarray, outarray, mkey);
            break;
        default:
            NEKERROR(ErrorUtil::efatal,
                     "This matrix does not have an operator");
            break;
    }
}

void StdExpansion::MassMatrixOp_MatFree(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const StdMatrixKey &mkey)
{
    int nq = GetTotPoints();
    Array<OneD, NekDouble> tmp(nq);

    v_BwdTrans(inarray, tmp);

    if (mkey.HasVarCoeff(eVarCoeffMass))
    {
        Vmath::Vmul(nq, mkey.GetVarCoeff(eVarCoeffMass), 1, tmp, 1, tmp, 1);
    }

    v_IProductWRTBase(tmp, outarray);
}

void StdExpansion::LaplacianMatrixOp_MatFree(
    const int k1, const int k2, const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const StdMatrixKey &mkey)
{
    ASSERTL1(k1 >= 0 && k1 < GetCoordim(), "invalid first  argument");
    ASSERTL1(k2 >= 0 && k2 < GetCoordim(), "invalid second argument");

    int nq = GetTotPoints();
    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> dtmp(nq);
    VarCoeffType varcoefftypes[3][3] = {
        {eVarCoeffD00, eVarCoeffD01, eVarCoeffD02},
        {eVarCoeffD10, eVarCoeffD11, eVarCoeffD12},
        {eVarCoeffD20, eVarCoeffD21, eVarCoeffD22}};

    ConstFactorType constcoefftypes[3][3] = {
        {eFactorCoeffD00, eFactorCoeffD01, eFactorCoeffD02},
        {eFactorCoeffD01, eFactorCoeffD11, eFactorCoeffD12},
        {eFactorCoeffD02, eFactorCoeffD12, eFactorCoeffD22}};

    v_BwdTrans(inarray, tmp);
    v_PhysDeriv(k2, tmp, dtmp);
    if (mkey.GetNVarCoeff() && (!mkey.ConstFactorExists(eFactorSVVDiffCoeff)))
    {
        if (k1 == k2)
        {
            // By default, k1 == k2 has \sigma = 1 (diagonal entries)
            if (mkey.HasVarCoeff(varcoefftypes[k1][k1]))
            {
                Vmath::Vmul(nq, mkey.GetVarCoeff(varcoefftypes[k1][k1]), 1,
                            dtmp, 1, dtmp, 1);
            }
            v_IProductWRTDerivBase_SumFac(k1, dtmp, outarray);
        }
        else
        {
            // By default, k1 != k2 has \sigma = 0 (off-diagonal entries)
            if (mkey.HasVarCoeff(varcoefftypes[k1][k2]))
            {
                Vmath::Vmul(nq, mkey.GetVarCoeff(varcoefftypes[k1][k2]), 1,
                            dtmp, 1, dtmp, 1);
                v_IProductWRTDerivBase_SumFac(k1, dtmp, outarray);
            }
            else if (mkey.HasVarCoeff(
                         varcoefftypes[k2][k1])) // Check symmetric varcoeff
            {
                Vmath::Vmul(nq, mkey.GetVarCoeff(varcoefftypes[k2][k1]), 1,
                            dtmp, 1, dtmp, 1);
                v_IProductWRTDerivBase_SumFac(k1, dtmp, outarray);
            }
            else
            {
                Vmath::Zero(GetNcoeffs(), outarray, 1);
            }
        }
    }
    else if (mkey.ConstFactorExists(eFactorCoeffD00) &&
             (!mkey.ConstFactorExists(eFactorSVVDiffCoeff)))
    {
        if (k1 == k2)
        {
            // By default, k1 == k2 has \sigma = 1 (diagonal entries)
            if (mkey.ConstFactorExists(constcoefftypes[k1][k1]))
            {
                Vmath::Smul(nq, mkey.GetConstFactor(constcoefftypes[k1][k1]),
                            dtmp, 1, dtmp, 1);
            }
            v_IProductWRTDerivBase(k1, dtmp, outarray);
        }
        else
        {
            // By default, k1 != k2 has \sigma = 0 (off-diagonal entries)
            if (mkey.ConstFactorExists(constcoefftypes[k1][k2]))
            {
                Vmath::Smul(nq, mkey.GetConstFactor(constcoefftypes[k1][k2]),
                            dtmp, 1, dtmp, 1);
                v_IProductWRTDerivBase(k1, dtmp, outarray);
            }
            else
            {
                Vmath::Zero(GetNcoeffs(), outarray, 1);
            }
        }
    }
    else
    {
        // Multiply by svv tensor
        if (mkey.ConstFactorExists(eFactorSVVDiffCoeff))
        {
            Vmath::Vcopy(nq, dtmp, 1, tmp, 1);
            SVVLaplacianFilter(dtmp, mkey);
            Vmath::Vadd(nq, tmp, 1, dtmp, 1, dtmp, 1);
        }
        v_IProductWRTDerivBase(k1, dtmp, outarray);
    }
}

void StdExpansion::LaplacianMatrixOp_MatFree_GenericImpl(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const StdMatrixKey &mkey)
{
    const int dim = GetCoordim();

    int i, j;

    Array<OneD, NekDouble> store(m_ncoeffs);
    Array<OneD, NekDouble> store2(m_ncoeffs, 0.0);

    if ((mkey.GetNVarCoeff() == 0 &&
         !mkey.ConstFactorExists(eFactorCoeffD00)) ||
        mkey.ConstFactorExists(eFactorSVVDiffCoeff))
    {
        // just call diagonal matrix form of laplcian operator
        for (i = 0; i < dim; ++i)
        {
            LaplacianMatrixOp(i, i, inarray, store, mkey);
            Vmath::Vadd(m_ncoeffs, store, 1, store2, 1, store2, 1);
        }
    }
    else
    {
        const MatrixType mtype[3][3] = {
            {eLaplacian00, eLaplacian01, eLaplacian02},
            {eLaplacian10, eLaplacian11, eLaplacian12},
            {eLaplacian20, eLaplacian21, eLaplacian22}};
        StdMatrixKeySharedPtr mkeyij;

        for (i = 0; i < dim; i++)
        {
            for (j = 0; j < dim; j++)
            {
                mkeyij = MemoryManager<StdMatrixKey>::AllocateSharedPtr(
                    mkey, mtype[i][j]);
                LaplacianMatrixOp(i, j, inarray, store, *mkeyij);
                Vmath::Vadd(m_ncoeffs, store, 1, store2, 1, store2, 1);
            }
        }
    }

    Vmath::Vcopy(m_ncoeffs, store2.data(), 1, outarray.data(), 1);
}

void StdExpansion::WeakDerivMatrixOp_MatFree(
    const int k1, const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const StdMatrixKey &mkey)
{
    Array<OneD, NekDouble> tmp(GetTotPoints());
    int nq = GetTotPoints();

    v_BwdTrans(inarray, tmp);
    v_PhysDeriv(k1, tmp, tmp);

    VarCoeffType keys[] = {eVarCoeffD00, eVarCoeffD11, eVarCoeffD22};
    if (mkey.HasVarCoeff(keys[k1]))
    {
        Vmath::Vmul(nq, &(mkey.GetVarCoeff(keys[k1]))[0], 1, &tmp[0], 1,
                    &tmp[0], 1);
    }

    v_IProductWRTBase(tmp, outarray);
}

void StdExpansion::WeakDirectionalDerivMatrixOp_MatFree(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const StdMatrixKey &mkey)
{
    int nq = GetTotPoints();

    Array<OneD, NekDouble> tmp(nq), Dtmp(nq);
    Array<OneD, NekDouble> Mtmp(nq), Mout(m_ncoeffs);

    v_BwdTrans(inarray, tmp);
    v_PhysDirectionalDeriv(tmp, mkey.GetVarCoeff(eVarCoeffMF), Dtmp);

    v_IProductWRTBase(Dtmp, outarray);

    // Compte M_{div tv}
    Vmath::Vmul(nq, &(mkey.GetVarCoeff(eVarCoeffMFDiv))[0], 1, &tmp[0], 1,
                &Mtmp[0], 1);

    v_IProductWRTBase(Mtmp, Mout);

    // Add D_tv + M_{div tv}
    Vmath::Vadd(m_ncoeffs, &Mout[0], 1, &outarray[0], 1, &outarray[0], 1);
}

void StdExpansion::MassLevelCurvatureMatrixOp_MatFree(
    [[maybe_unused]] const Array<OneD, const NekDouble> &inarray,
    [[maybe_unused]] Array<OneD, NekDouble> &outarray,
    [[maybe_unused]] const StdMatrixKey &mkey)
{
}

void StdExpansion::LinearAdvectionMatrixOp_MatFree(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const StdMatrixKey &mkey)
{
    int i, ndir = 0;

    VarCoeffType varcoefftypes[] = {eVarCoeffVelX, eVarCoeffVelY,
                                    eVarCoeffVelZ};
    // Count advection velocities
    for (auto &x : varcoefftypes)
    {
        if (mkey.HasVarCoeff(x))
        {
            ndir++;
        }
    }

    ASSERTL0(ndir, "Must define at least one advection velocity");
    ASSERTL1(ndir <= GetCoordim(),
             "Number of constants is larger than coordinate dimensions");

    int totpts = GetTotPoints();
    Array<OneD, NekDouble> tmp(3 * totpts);
    Array<OneD, NekDouble> tmp_deriv = tmp + totpts;
    Array<OneD, NekDouble> tmp_adv   = tmp_deriv + totpts;

    v_BwdTrans(inarray, tmp); // transform to PhysSpace

    // Evaluate advection (u dx + v dy + w dz)
    Vmath::Zero(totpts, tmp_adv, 1);
    for (i = 0; i < ndir; ++i)
    {
        v_PhysDeriv(i, tmp, tmp_deriv);
        Vmath::Vvtvp(totpts, mkey.GetVarCoeff(varcoefftypes[i]), 1, tmp_deriv,
                     1, tmp_adv, 1, tmp_adv, 1);
    }

    v_IProductWRTBase(tmp_adv, outarray);
}

void StdExpansion::LinearAdvectionDiffusionReactionMatrixOp_MatFree(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const StdMatrixKey &mkey,
    bool addDiffusionTerm)
{
    int i, ndir = 0;

    VarCoeffType varcoefftypes[] = {eVarCoeffVelX, eVarCoeffVelY,
                                    eVarCoeffVelZ};
    // Count advection velocities
    for (auto &x : varcoefftypes)
    {
        if (mkey.HasVarCoeff(x))
        {
            ndir++;
        }
    }

    ASSERTL0(ndir, "Must define at least one advection velocity");
    ASSERTL1(ndir <= GetCoordim(),
             "Number of constants is larger than coordinate dimensions");

    NekDouble lambda = mkey.GetConstFactor(eFactorLambda);
    int totpts       = GetTotPoints();
    Array<OneD, NekDouble> tmp(3 * totpts);
    Array<OneD, NekDouble> tmp_deriv = tmp + totpts;
    Array<OneD, NekDouble> tmp_adv   = tmp_deriv + totpts;

    v_BwdTrans(inarray, tmp); // transform this mode \phi_i into PhysSpace

    // calculate advection u dx + v dy + ..
    Vmath::Zero(totpts, tmp_adv, 1);
    for (i = 0; i < ndir; ++i)
    {
        v_PhysDeriv(i, tmp, tmp_deriv);
        Vmath::Vvtvp(totpts, mkey.GetVarCoeff(varcoefftypes[i]), 1, tmp_deriv,
                     1, tmp_adv, 1, tmp_adv, 1);
    }

    // Add reaction term if lambda != 0.0
    if (lambda)
    {
        // Add mass varcoeff
        if (mkey.HasVarCoeff(eVarCoeffMass))
        {
            Vmath::Vmul(totpts, mkey.GetVarCoeff(eVarCoeffMass), 1, tmp, 1, tmp,
                        1);
        }

        Vmath::Svtvp(totpts, -lambda, tmp, 1, tmp_adv, 1, tmp_adv, 1);
    }

    // Create mass matrix = Advection - Reaction
    v_IProductWRTBase(tmp_adv,
                      outarray); // Create mass matrix of Advection - Reaction

    // Add Laplacian matrix
    if (addDiffusionTerm)
    {
        Array<OneD, NekDouble> lap(m_ncoeffs);
        StdMatrixKey mkeylap(eLaplacian, DetShapeType(), *this,
                             mkey.GetConstFactors(), mkey.GetVarCoeffs(),
                             mkey.GetNodalPointsType());
        LaplacianMatrixOp(inarray, lap, mkeylap);

        Vmath::Vadd(m_ncoeffs, lap, 1, outarray, 1, outarray,
                    1); // += Laplacian
    }
}

void StdExpansion::HelmholtzMatrixOp_MatFree_GenericImpl(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const StdMatrixKey &mkey)
{
    NekDouble lambda = mkey.GetConstFactor(eFactorLambda);
    Array<OneD, NekDouble> tmp(m_ncoeffs);
    StdMatrixKey mkeymass(eMass, DetShapeType(), *this);
    StdMatrixKey mkeylap(eLaplacian, DetShapeType(), *this,
                         mkey.GetConstFactors(), mkey.GetVarCoeffs(),
                         mkey.GetNodalPointsType());

    MassMatrixOp(inarray, tmp, mkeymass);
    LaplacianMatrixOp(inarray, outarray, mkeylap);

    Blas::Daxpy(m_ncoeffs, lambda, tmp, 1, outarray, 1);
}

// VIRTUAL INLINE FUNCTIONS FROM HEADER FILE
NekDouble StdExpansion::StdPhysEvaluate(
    const Array<OneD, const NekDouble> &Lcoord,
    const Array<OneD, const NekDouble> &physvals)
{
    return v_StdPhysEvaluate(Lcoord, physvals);
}

int StdExpansion::v_CalcNumberOfCoefficients(
    [[maybe_unused]] const std::vector<unsigned int> &nummodes,
    [[maybe_unused]] int &modes_offset)
{
    NEKERROR(ErrorUtil::efatal, "This function is not defined for this class");
    return 0;
}

void StdExpansion::v_NormVectorIProductWRTBase(
    [[maybe_unused]] const Array<OneD, const NekDouble> &Fx,
    [[maybe_unused]] Array<OneD, NekDouble> &outarray)
{
    NEKERROR(ErrorUtil::efatal, "This function is not valid for this class");
}

void StdExpansion::v_NormVectorIProductWRTBase(
    [[maybe_unused]] const Array<OneD, const NekDouble> &Fx,
    [[maybe_unused]] const Array<OneD, const NekDouble> &Fy,
    [[maybe_unused]] Array<OneD, NekDouble> &outarray)
{
    NEKERROR(ErrorUtil::efatal, "This function is not valid for this class");
}

void StdExpansion::v_NormVectorIProductWRTBase(
    [[maybe_unused]] const Array<OneD, const NekDouble> &Fx,
    [[maybe_unused]] const Array<OneD, const NekDouble> &Fy,
    [[maybe_unused]] const Array<OneD, const NekDouble> &Fz,
    [[maybe_unused]] Array<OneD, NekDouble> &outarray)
{
    NEKERROR(ErrorUtil::efatal, "This function is not valid for this class");
}

void StdExpansion::v_NormVectorIProductWRTBase(
    [[maybe_unused]] const Array<OneD, const Array<OneD, NekDouble>> &Fvec,
    [[maybe_unused]] Array<OneD, NekDouble> &outarray)
{
    NEKERROR(ErrorUtil::efatal, "This function is not valid for this class");
}

DNekScalBlkMatSharedPtr StdExpansion::v_GetLocStaticCondMatrix(
    [[maybe_unused]] const LocalRegions::MatrixKey &mkey)
{
    NEKERROR(ErrorUtil::efatal, "This function is only valid for LocalRegions");
    return NullDNekScalBlkMatSharedPtr;
}

void StdExpansion::v_DropLocStaticCondMatrix(
    [[maybe_unused]] const LocalRegions::MatrixKey &mkey)
{
    NEKERROR(ErrorUtil::efatal, "This function is only valid for LocalRegions");
}

void StdExpansion::v_SetCoeffsToOrientation(
    [[maybe_unused]] StdRegions::Orientation dir,
    [[maybe_unused]] Array<OneD, const NekDouble> &inarray,
    [[maybe_unused]] Array<OneD, NekDouble> &outarray)
{
    NEKERROR(ErrorUtil::efatal, "This function is not defined for this shape");
}

NekDouble StdExpansion::v_StdPhysEvaluate(
    [[maybe_unused]] const Array<OneD, const NekDouble> &Lcoord,
    [[maybe_unused]] const Array<OneD, const NekDouble> &physvals)

{
    NEKERROR(ErrorUtil::efatal, "This function is not defined for this shape");
    return 0;
}

void StdExpansion::v_LocCoordToLocCollapsed(
    [[maybe_unused]] const Array<OneD, const NekDouble> &xi,
    [[maybe_unused]] Array<OneD, NekDouble> &eta)
{
    NEKERROR(ErrorUtil::efatal, "This function is not defined for this shape");
}

void StdExpansion::v_LocCollapsedToLocCoord(
    [[maybe_unused]] const Array<OneD, const NekDouble> &eta,
    [[maybe_unused]] Array<OneD, NekDouble> &xi)
{
    NEKERROR(ErrorUtil::efatal, "This function is not defined for this shape");
}

void StdExpansion::v_PhysInterp(
    [[maybe_unused]] std::shared_ptr<StdExpansion> FromExp,
    [[maybe_unused]] const Array<OneD, const NekDouble> &fromData,
    [[maybe_unused]] Array<OneD, NekDouble> &toData)
{
    ASSERTL0(false, "This function is not valid or not defined");
}

const LibUtilities::BasisKey StdExpansion::v_GetTraceBasisKey(
    [[maybe_unused]] const int i, [[maybe_unused]] const int k,
    [[maybe_unused]] bool UseGLL) const
{
    ASSERTL0(false, "This function is not valid or not defined");
    return LibUtilities::NullBasisKey;
}

LibUtilities::PointsKey StdExpansion::v_GetTracePointsKey(
    [[maybe_unused]] const int i, [[maybe_unused]] const int j) const
{
    ASSERTL0(false, "This function is not valid or not defined");
    return LibUtilities::NullPointsKey;
}

const LibUtilities::PointsKey StdExpansion::v_GetNodalPointsKey() const
{
    ASSERTL0(false, "This function is not valid or not defined");

    return LibUtilities::NullPointsKey;
}

std::shared_ptr<StdExpansion> StdExpansion::v_GetStdExp(void) const
{
    ASSERTL0(false, "This method is not defined for this expansion");
    StdExpansionSharedPtr returnval;
    return returnval;
}

std::shared_ptr<StdExpansion> StdExpansion::v_GetLinStdExp(void) const
{
    ASSERTL0(false, "This method is not defined for this expansion");
    StdExpansionSharedPtr returnval;
    return returnval;
}

bool StdExpansion::v_IsBoundaryInteriorExpansion() const
{
    ASSERTL0(false, "This function has not been defined for this expansion");
    return false;
}

bool StdExpansion::v_IsNodalNonTensorialExp()
{
    return false;
}

void StdExpansion::v_IProductWRTDerivBase(
    [[maybe_unused]] const int dir,
    [[maybe_unused]] const Array<OneD, const NekDouble> &inarray,
    [[maybe_unused]] Array<OneD, NekDouble> &outarray)
{
    NEKERROR(ErrorUtil::efatal, "This method has not been defined");
}

/**
 *
 */
void StdExpansion::v_IProductWRTDirectionalDerivBase(
    [[maybe_unused]] const Array<OneD, const NekDouble> &direction,
    [[maybe_unused]] const Array<OneD, const NekDouble> &inarray,
    [[maybe_unused]] Array<OneD, NekDouble> &outarray)
{
    NEKERROR(ErrorUtil::efatal, "This method has not been defined");
}

/**
 *
 */
void StdExpansion::v_FwdTransBndConstrained(
    [[maybe_unused]] const Array<OneD, const NekDouble> &inarray,
    [[maybe_unused]] Array<OneD, NekDouble> &outarray)
{
    NEKERROR(ErrorUtil::efatal, "This method has not been defined");
}

/**
 * @brief Integrates the specified function over the domain.
 * @see StdRegions#StdExpansion#Integral.
 */
NekDouble StdExpansion::v_Integral(
    [[maybe_unused]] const Array<OneD, const NekDouble> &inarray)
{
    NEKERROR(ErrorUtil::efatal, "This function is only valid for "
                                "local expansions");
    return 0;
}

/**
 * @brief Calculate the derivative of the physical points
 * @see StdRegions#StdExpansion#PhysDeriv
 */
void StdExpansion::v_PhysDeriv(
    [[maybe_unused]] const Array<OneD, const NekDouble> &inarray,
    [[maybe_unused]] Array<OneD, NekDouble> &out_d1,
    [[maybe_unused]] Array<OneD, NekDouble> &out_d2,
    [[maybe_unused]] Array<OneD, NekDouble> &out_d3)
{
    NEKERROR(ErrorUtil::efatal, "This function is only valid for "
                                "local expansions");
}

void StdExpansion::v_PhysDeriv_s(
    [[maybe_unused]] const Array<OneD, const NekDouble> &inarray,
    [[maybe_unused]] Array<OneD, NekDouble> &out_ds)
{
    NEKERROR(ErrorUtil::efatal, "This function is only valid for "
                                "local expansions");
}
void StdExpansion::v_PhysDeriv_n(
    [[maybe_unused]] const Array<OneD, const NekDouble> &inarray,
    [[maybe_unused]] Array<OneD, NekDouble> &out_dn)
{
    NEKERROR(ErrorUtil::efatal, "This function is only valid for "
                                "local expansions");
}

/**
 * @brief Calculate the derivative of the physical points in a
 * given direction
 * @see StdRegions#StdExpansion#PhysDeriv
 */
void StdExpansion::v_PhysDeriv(
    [[maybe_unused]] const int dir,
    [[maybe_unused]] const Array<OneD, const NekDouble> &inarray,
    [[maybe_unused]] Array<OneD, NekDouble> &out_d0)

{
    NEKERROR(ErrorUtil::efatal, "This function is only valid for "
                                "specific element types");
}

/**
 * @brief Physical derivative along a direction vector.
 * @see StdRegions#StdExpansion#PhysDirectionalDeriv
 */
void StdExpansion::v_PhysDirectionalDeriv(
    [[maybe_unused]] const Array<OneD, const NekDouble> &inarray,
    [[maybe_unused]] const Array<OneD, const NekDouble> &direction,
    [[maybe_unused]] Array<OneD, NekDouble> &outarray)
{
    NEKERROR(ErrorUtil::efatal, "This function is only valid for "
                                "specific element types");
}

void StdExpansion::v_StdPhysDeriv(
    [[maybe_unused]] const Array<OneD, const NekDouble> &inarray,
    [[maybe_unused]] Array<OneD, NekDouble> &out_d1,
    [[maybe_unused]] Array<OneD, NekDouble> &out_d2,
    [[maybe_unused]] Array<OneD, NekDouble> &out_d3)
{
    NEKERROR(ErrorUtil::efatal, "Method does not exist for this shape");
}

void StdExpansion::v_StdPhysDeriv(
    [[maybe_unused]] const int dir,
    [[maybe_unused]] const Array<OneD, const NekDouble> &inarray,
    [[maybe_unused]] Array<OneD, NekDouble> &outarray)
{
    NEKERROR(ErrorUtil::efatal, "Method does not exist for this shape");
}

NekDouble StdExpansion::v_PhysEvaluateBasis(
    [[maybe_unused]] const Array<OneD, const NekDouble> &coords,
    [[maybe_unused]] int mode)
{
    NEKERROR(ErrorUtil::efatal, "Method does not exist for this shape");
    return 0;
}

NekDouble StdExpansion::v_PhysEvaluate(
    [[maybe_unused]] const Array<OneD, const NekDouble> &coords,
    [[maybe_unused]] const Array<OneD, const NekDouble> &physvals)
{
    NEKERROR(ErrorUtil::efatal, "Method does not exist for this shape");
    return 0;
}

NekDouble StdExpansion::v_PhysEvaluateInterp(
    [[maybe_unused]] const Array<OneD, DNekMatSharedPtr> &I,
    [[maybe_unused]] const Array<OneD, const NekDouble> &physvals)
{
    NEKERROR(ErrorUtil::efatal, "Method does not exist for this shape");
    return 0;
}

NekDouble StdExpansion::v_PhysEvalFirstDeriv(
    [[maybe_unused]] const Array<OneD, NekDouble> &coord,
    [[maybe_unused]] const Array<OneD, const NekDouble> &inarray,
    [[maybe_unused]] std::array<NekDouble, 3> &firstOrderDerivs)
{
    NEKERROR(ErrorUtil::efatal,
             "PhysEvaluate first order derivative method does not exist"
             " for this shape type: " +
                 static_cast<std::string>(
                     LibUtilities::ShapeTypeMap[DetShapeType()]));
    return 0;
}

NekDouble StdExpansion::v_PhysEvalFirstSecondDeriv(
    [[maybe_unused]] const Array<OneD, NekDouble> &coord,
    [[maybe_unused]] const Array<OneD, const NekDouble> &inarray,
    [[maybe_unused]] std::array<NekDouble, 3> &firstOrderDerivs,
    [[maybe_unused]] std::array<NekDouble, 6> &secondOrderDerivs)
{
    NEKERROR(ErrorUtil::efatal,
             "PhysEvaluate second order derivative method does not exist"
             " for this shape type: " +
                 static_cast<std::string>(
                     LibUtilities::ShapeTypeMap[DetShapeType()]));
    return 0;
}

void StdExpansion::v_FillMode([[maybe_unused]] const int mode,
                              [[maybe_unused]] Array<OneD, NekDouble> &outarray)
{
    NEKERROR(ErrorUtil::efatal, "This function has not "
                                "been defined for this shape");
}

DNekMatSharedPtr StdExpansion::v_GenMatrix(
    [[maybe_unused]] const StdMatrixKey &mkey)
{
    NEKERROR(ErrorUtil::efatal, "This function has not "
                                "been defined for this element");
    DNekMatSharedPtr returnval;
    return returnval;
}

DNekMatSharedPtr StdExpansion::v_CreateStdMatrix(
    [[maybe_unused]] const StdMatrixKey &mkey)
{
    NEKERROR(ErrorUtil::efatal, "This function has not "
                                "been defined for this element");
    DNekMatSharedPtr returnval;
    return returnval;
}

void StdExpansion::v_GetCoords(
    [[maybe_unused]] Array<OneD, NekDouble> &coords_0,
    [[maybe_unused]] Array<OneD, NekDouble> &coords_1,
    [[maybe_unused]] Array<OneD, NekDouble> &coords_2)
{
    NEKERROR(ErrorUtil::efatal, "Write coordinate definition method");
}

void StdExpansion::v_GetCoord(
    [[maybe_unused]] const Array<OneD, const NekDouble> &Lcoord,
    [[maybe_unused]] Array<OneD, NekDouble> &coord)
{
    NEKERROR(ErrorUtil::efatal, "Write coordinate definition method");
}

int StdExpansion::v_GetCoordim() const
{
    return GetShapeDimension();
}

void StdExpansion::v_GetBoundaryMap(
    [[maybe_unused]] Array<OneD, unsigned int> &outarray)
{
    NEKERROR(ErrorUtil::efatal, "Method does not exist for this shape");
}

void StdExpansion::v_GetInteriorMap(
    [[maybe_unused]] Array<OneD, unsigned int> &outarray)
{
    NEKERROR(ErrorUtil::efatal, "Method does not exist for this shape");
}

int StdExpansion::v_GetVertexMap([[maybe_unused]] const int localVertexId,
                                 [[maybe_unused]] bool useCoeffPacking)
{
    NEKERROR(ErrorUtil::efatal, "Method does not exist for this shape");
    return 0;
}

void StdExpansion::v_GetTraceToElementMap(
    [[maybe_unused]] const int tid,
    [[maybe_unused]] Array<OneD, unsigned int> &maparray,
    [[maybe_unused]] Array<OneD, int> &signarray,
    [[maybe_unused]] Orientation traceOrient, [[maybe_unused]] int P,
    [[maybe_unused]] int Q)
{
    NEKERROR(ErrorUtil::efatal, "Method does not exist for this shape");
}

void StdExpansion::v_GetTraceCoeffMap(
    [[maybe_unused]] const unsigned int traceid,
    [[maybe_unused]] Array<OneD, unsigned int> &maparray)
{
    NEKERROR(ErrorUtil::efatal, "Method does not exist for this shape");
}

void StdExpansion::v_GetElmtTraceToTraceMap(
    [[maybe_unused]] const unsigned int tid,
    [[maybe_unused]] Array<OneD, unsigned int> &maparray,
    [[maybe_unused]] Array<OneD, int> &signarray,
    [[maybe_unused]] Orientation traceOrient, [[maybe_unused]] int P,
    [[maybe_unused]] int Q)
{
    NEKERROR(ErrorUtil::efatal, "Method does not exist for this shape");
}

void StdExpansion::v_GetTraceInteriorToElementMap(
    [[maybe_unused]] const int tid,
    [[maybe_unused]] Array<OneD, unsigned int> &maparray,
    [[maybe_unused]] Array<OneD, int> &signarray,
    [[maybe_unused]] const Orientation traceOrient)
{
    NEKERROR(ErrorUtil::efatal, "Method does not exist for this shape");
}

void StdExpansion::v_GetTraceNumModes([[maybe_unused]] const int tid,
                                      [[maybe_unused]] int &numModes0,
                                      [[maybe_unused]] int &numModes1,
                                      [[maybe_unused]] Orientation traceOrient)
{
    NEKERROR(ErrorUtil::efatal, "Method does not exist for this shape");
}

void StdExpansion::v_GetVertexPhysVals(
    [[maybe_unused]] const int vertex,
    [[maybe_unused]] const Array<OneD, const NekDouble> &inarray,
    [[maybe_unused]] NekDouble &outarray)
{
    NEKERROR(ErrorUtil::efatal, "Method does not exist for "
                                "this shape or library");
}

void StdExpansion::v_MultiplyByQuadratureMetric(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray)
{
    v_MultiplyByStdQuadratureMetric(inarray, outarray);
}

void StdExpansion::v_MultiplyByStdQuadratureMetric(
    [[maybe_unused]] const Array<OneD, const NekDouble> &inarray,
    [[maybe_unused]] Array<OneD, NekDouble> &outarray)
{
    NEKERROR(ErrorUtil::efatal,
             "Method does not exist for this shape or library");
}

void StdExpansion::v_BwdTrans_SumFac(
    [[maybe_unused]] const Array<OneD, const NekDouble> &inarray,
    [[maybe_unused]] Array<OneD, NekDouble> &outarray)
{
    NEKERROR(ErrorUtil::efatal, "Method does not exist for this shape");
}

void StdExpansion::v_IProductWRTBase_SumFac(
    [[maybe_unused]] const Array<OneD, const NekDouble> &inarray,
    [[maybe_unused]] Array<OneD, NekDouble> &outarray,
    [[maybe_unused]] bool multiplybyweights)
{
    NEKERROR(ErrorUtil::efatal, "Method does not exist for this shape");
}

/**
 *
 */
void StdExpansion::v_IProductWRTDirectionalDerivBase_SumFac(
    [[maybe_unused]] const Array<OneD, const NekDouble> &direction,
    [[maybe_unused]] const Array<OneD, const NekDouble> &inarray,
    [[maybe_unused]] Array<OneD, NekDouble> &outarray)
{
    NEKERROR(ErrorUtil::efatal, "Method does not exist for this shape");
}

void StdExpansion::v_IProductWRTDerivBase_SumFac(
    [[maybe_unused]] const int dir,
    [[maybe_unused]] const Array<OneD, const NekDouble> &inarray,
    [[maybe_unused]] Array<OneD, NekDouble> &outarray)
{
    NEKERROR(ErrorUtil::efatal, "Method does not exist for this shape");
}

void StdExpansion::v_MassMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                  Array<OneD, NekDouble> &outarray,
                                  const StdMatrixKey &mkey)
{
    // If this function is not reimplemented on shape level, the function
    // below will be called
    MassMatrixOp_MatFree(inarray, outarray, mkey);
}

void StdExpansion::v_LaplacianMatrixOp(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const StdMatrixKey &mkey)
{
    // If this function is not reimplemented on shape level, the function
    // below will be called
    LaplacianMatrixOp_MatFree(inarray, outarray, mkey);
}

void StdExpansion::v_SVVLaplacianFilter(
    [[maybe_unused]] Array<OneD, NekDouble> &array,
    [[maybe_unused]] const StdMatrixKey &mkey)
{
    ASSERTL0(false, "This function is not defined in StdExpansion.");
}

void StdExpansion::v_ExponentialFilter(
    [[maybe_unused]] Array<OneD, NekDouble> &array,
    [[maybe_unused]] const NekDouble alpha,
    [[maybe_unused]] const NekDouble exponent,
    [[maybe_unused]] const NekDouble cutoff)
{
    ASSERTL0(false, "This function is not defined in StdExpansion.");
}

void StdExpansion::v_ReduceOrderCoeffs(
    [[maybe_unused]] int numMin,
    [[maybe_unused]] const Array<OneD, const NekDouble> &inarray,
    [[maybe_unused]] Array<OneD, NekDouble> &outarray)
{
    ASSERTL0(false, "This function is not defined in StdExpansion.");
}

void StdExpansion::v_LaplacianMatrixOp(
    const int k1, const int k2, const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const StdMatrixKey &mkey)
{
    // If this function is not reimplemented on shape level, the function
    // below will be called
    LaplacianMatrixOp_MatFree(k1, k2, inarray, outarray, mkey);
}

void StdExpansion::v_WeakDerivMatrixOp(
    const int i, const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const StdMatrixKey &mkey)
{
    // If this function is not reimplemented on shape level, the function
    // below will be called
    WeakDerivMatrixOp_MatFree(i, inarray, outarray, mkey);
}

void StdExpansion::v_WeakDirectionalDerivMatrixOp(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const StdMatrixKey &mkey)
{
    // If this function is not reimplemented on shape level, the function
    // below will be called
    WeakDirectionalDerivMatrixOp_MatFree(inarray, outarray, mkey);
}

void StdExpansion::v_MassLevelCurvatureMatrixOp(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const StdMatrixKey &mkey)
{
    // If this function is not reimplemented on shape level, the function
    // below will be called
    MassLevelCurvatureMatrixOp_MatFree(inarray, outarray, mkey);
}

void StdExpansion::v_LinearAdvectionMatrixOp(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const StdMatrixKey &mkey)
{
    // If this function is not reimplemented on shape level, the function
    // below will be called
    LinearAdvectionMatrixOp_MatFree(inarray, outarray, mkey);
}

void StdExpansion::v_LinearAdvectionDiffusionReactionMatrixOp(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const StdMatrixKey &mkey,
    bool addDiffusionTerm)
{
    // If this function is not reimplemented on shape level, the function
    // below will be called
    LinearAdvectionDiffusionReactionMatrixOp_MatFree(inarray, outarray, mkey,
                                                     addDiffusionTerm);
}

void StdExpansion::v_HelmholtzMatrixOp(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const StdMatrixKey &mkey)
{
    // If this function is not reimplemented on shape level, the function
    // below will be called
    HelmholtzMatrixOp_MatFree(inarray, outarray, mkey);
}

void StdExpansion::v_LaplacianMatrixOp_MatFree(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const StdMatrixKey &mkey)
{
    // If this function is not reimplemented on shape level, the function
    // below will be called
    LaplacianMatrixOp_MatFree_GenericImpl(inarray, outarray, mkey);
}

void StdExpansion::v_LaplacianMatrixOp_MatFree_Kernel(
    [[maybe_unused]] const Array<OneD, const NekDouble> &inarray,
    [[maybe_unused]] Array<OneD, NekDouble> &outarray,
    [[maybe_unused]] Array<OneD, NekDouble> &wsp)
{
    ASSERTL0(false, "Not implemented.");
}

void StdExpansion::v_HelmholtzMatrixOp_MatFree(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const StdMatrixKey &mkey)
{
    // If this function is not reimplemented on shape level, the function
    // below will be called
    HelmholtzMatrixOp_MatFree_GenericImpl(inarray, outarray, mkey);
}

DNekMatSharedPtr StdExpansion::v_BuildInverseTransformationMatrix(
    [[maybe_unused]] const DNekScalMatSharedPtr &m_transformationmatrix)
{
    NEKERROR(ErrorUtil::efatal, "This function is only valid for LocalRegions");
    return NullDNekMatSharedPtr;
}

void StdExpansion::PhysInterpToSimplexEquiSpaced(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, int npset)
{
    LibUtilities::ShapeType shape = DetShapeType();
    DNekMatSharedPtr intmat;

    int nqtot = GetTotPoints();
    int np    = 0;
    if (npset == -1) // use values from basis num points()
    {
        int nqbase;
        for (int i = 0; i < m_base.size(); ++i)
        {
            nqbase = m_base[i]->GetNumPoints();
            np     = std::max(np, nqbase);
        }

        StdMatrixKey Ikey(ePhysInterpToEquiSpaced, shape, *this);
        intmat = GetStdMatrix(Ikey);
    }
    else
    {
        np = npset;

        ConstFactorMap cmap;
        cmap[eFactorConst] = np;
        StdMatrixKey Ikey(ePhysInterpToEquiSpaced, shape, *this, cmap);
        intmat = GetStdMatrix(Ikey);
    }

    NekVector<NekDouble> in(nqtot, inarray, eWrapper);
    NekVector<NekDouble> out(
        LibUtilities::GetNumberOfCoefficients(shape, np, np, np), outarray,
        eWrapper);
    out = (*intmat) * in;
}

void StdExpansion::v_GetSimplexEquiSpacedConnectivity(
    [[maybe_unused]] Array<OneD, int> &conn, [[maybe_unused]] bool standard)
{
    NEKERROR(ErrorUtil::efatal,
             "GetSimplexEquiSpacedConnectivity not"
             " implemented for " +
                 static_cast<std::string>(
                     LibUtilities::ShapeTypeMap[DetShapeType()]));
}

void StdExpansion::EquiSpacedToCoeffs(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray)
{
    LibUtilities::ShapeType shape = DetShapeType();

    // inarray has to be consistent with NumModes definition
    // There is also a check in GetStdMatrix to see if all
    // modes are of the same size
    ConstFactorMap cmap;

    cmap[eFactorConst] = m_base[0]->GetNumModes();
    StdMatrixKey Ikey(eEquiSpacedToCoeffs, shape, *this, cmap);
    DNekMatSharedPtr intmat = GetStdMatrix(Ikey);

    NekVector<NekDouble> in(m_ncoeffs, inarray, eWrapper);
    NekVector<NekDouble> out(m_ncoeffs, outarray, eWrapper);
    out = (*intmat) * in;
}

} // namespace Nektar::StdRegions
