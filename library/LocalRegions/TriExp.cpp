///////////////////////////////////////////////////////////////////////////////
//
// File: TriExp.cpp
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
// Description: Expasion for triangular elements.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Foundations/Interp.h>
#include <LibUtilities/Foundations/InterpCoeff.h>
#include <LocalRegions/Expansion3D.h>
#include <LocalRegions/SegExp.h>
#include <LocalRegions/TriExp.h>
#include <StdRegions/StdNodalTriExp.h>

using namespace std;

namespace Nektar::LocalRegions
{
TriExp::TriExp(const LibUtilities::BasisKey &Ba,
               const LibUtilities::BasisKey &Bb,
               const SpatialDomains::Geometry2DSharedPtr &geom)
    : StdExpansion(LibUtilities::StdTriData::getNumberOfCoefficients(
                       Ba.GetNumModes(), (Bb.GetNumModes())),
                   2, Ba, Bb),
      StdExpansion2D(LibUtilities::StdTriData::getNumberOfCoefficients(
                         Ba.GetNumModes(), (Bb.GetNumModes())),
                     Ba, Bb),
      StdTriExp(Ba, Bb), Expansion(geom), Expansion2D(geom),
      m_matrixManager(
          std::bind(&Expansion2D::CreateMatrix, this, std::placeholders::_1),
          std::string("TriExpMatrix")),
      m_staticCondMatrixManager(std::bind(&Expansion::CreateStaticCondMatrix,
                                          this, std::placeholders::_1),
                                std::string("TriExpStaticCondMatrix"))
{
}

TriExp::TriExp(const TriExp &T)
    : StdExpansion(T), StdExpansion2D(T), StdTriExp(T), Expansion(T),
      Expansion2D(T), m_matrixManager(T.m_matrixManager),
      m_staticCondMatrixManager(T.m_staticCondMatrixManager)
{
}

NekDouble TriExp::v_Integral(const Array<OneD, const NekDouble> &inarray)
{
    int nquad0                       = m_base[0]->GetNumPoints();
    int nquad1                       = m_base[1]->GetNumPoints();
    Array<OneD, const NekDouble> jac = m_metricinfo->GetJac(GetPointsKeys());
    NekDouble ival;
    Array<OneD, NekDouble> tmp(nquad0 * nquad1);

    // multiply inarray with Jacobian
    if (m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
    {
        Vmath::Vmul(nquad0 * nquad1, jac, 1, inarray, 1, tmp, 1);
    }
    else
    {
        Vmath::Smul(nquad0 * nquad1, jac[0], inarray, 1, tmp, 1);
    }

    // call StdQuadExp version;
    ival = StdTriExp::v_Integral(tmp);
    return ival;
}

void TriExp::v_PhysDeriv(const Array<OneD, const NekDouble> &inarray,
                         Array<OneD, NekDouble> &out_d0,
                         Array<OneD, NekDouble> &out_d1,
                         Array<OneD, NekDouble> &out_d2)
{
    int nquad0 = m_base[0]->GetNumPoints();
    int nquad1 = m_base[1]->GetNumPoints();
    int nqtot  = nquad0 * nquad1;
    const Array<TwoD, const NekDouble> &df =
        m_metricinfo->GetDerivFactors(GetPointsKeys());

    Array<OneD, NekDouble> diff0(2 * nqtot);
    Array<OneD, NekDouble> diff1(diff0 + nqtot);

    StdTriExp::v_PhysDeriv(inarray, diff0, diff1);

    if (m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
    {
        if (out_d0.size())
        {
            Vmath::Vmul(nqtot, df[0], 1, diff0, 1, out_d0, 1);
            Vmath::Vvtvp(nqtot, df[1], 1, diff1, 1, out_d0, 1, out_d0, 1);
        }

        if (out_d1.size())
        {
            Vmath::Vmul(nqtot, df[2], 1, diff0, 1, out_d1, 1);
            Vmath::Vvtvp(nqtot, df[3], 1, diff1, 1, out_d1, 1, out_d1, 1);
        }

        if (out_d2.size())
        {
            Vmath::Vmul(nqtot, df[4], 1, diff0, 1, out_d2, 1);
            Vmath::Vvtvp(nqtot, df[5], 1, diff1, 1, out_d2, 1, out_d2, 1);
        }
    }
    else // regular geometry
    {
        if (out_d0.size())
        {
            Vmath::Smul(nqtot, df[0][0], diff0, 1, out_d0, 1);
            Blas::Daxpy(nqtot, df[1][0], diff1, 1, out_d0, 1);
        }

        if (out_d1.size())
        {
            Vmath::Smul(nqtot, df[2][0], diff0, 1, out_d1, 1);
            Blas::Daxpy(nqtot, df[3][0], diff1, 1, out_d1, 1);
        }

        if (out_d2.size())
        {
            Vmath::Smul(nqtot, df[4][0], diff0, 1, out_d2, 1);
            Blas::Daxpy(nqtot, df[5][0], diff1, 1, out_d2, 1);
        }
    }
}

void TriExp::v_PhysDeriv(const int dir,
                         const Array<OneD, const NekDouble> &inarray,
                         Array<OneD, NekDouble> &outarray)
{
    switch (dir)
    {
        case 0:
        {
            PhysDeriv(inarray, outarray, NullNekDouble1DArray,
                      NullNekDouble1DArray);
        }
        break;
        case 1:
        {
            PhysDeriv(inarray, NullNekDouble1DArray, outarray,
                      NullNekDouble1DArray);
        }
        break;
        case 2:
        {
            PhysDeriv(inarray, NullNekDouble1DArray, NullNekDouble1DArray,
                      outarray);
        }
        break;
        default:
        {
            ASSERTL1(false, "input dir is out of range");
        }
        break;
    }
}

void TriExp::v_PhysDirectionalDeriv(
    const Array<OneD, const NekDouble> &inarray,
    const Array<OneD, const NekDouble> &direction, Array<OneD, NekDouble> &out)
{
    if (!out.size())
    {
        return;
    }

    int nquad0 = m_base[0]->GetNumPoints();
    int nquad1 = m_base[1]->GetNumPoints();
    int nqtot  = nquad0 * nquad1;

    const Array<TwoD, const NekDouble> &df =
        m_metricinfo->GetDerivFactors(GetPointsKeys());

    Array<OneD, NekDouble> diff0(2 * nqtot);
    Array<OneD, NekDouble> diff1(diff0 + nqtot);

    // diff0 = du/d_xi, diff1 = du/d_eta
    StdTriExp::v_PhysDeriv(inarray, diff0, diff1);

    if (m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
    {
        Array<OneD, Array<OneD, NekDouble>> tangmat(2);

        // D^v_xi = v_x*d_xi/dx + v_y*d_xi/dy + v_z*d_xi/dz
        // D^v_eta = v_x*d_eta/dx + v_y*d_eta/dy + v_z*d_eta/dz
        for (int i = 0; i < 2; ++i)
        {
            tangmat[i] = Array<OneD, NekDouble>(nqtot, 0.0);
            for (int k = 0; k < (m_geom->GetCoordim()); ++k)
            {
                Vmath::Vvtvp(nqtot, &df[2 * k + i][0], 1, &direction[k * nqtot],
                             1, &tangmat[i][0], 1, &tangmat[i][0], 1);
            }
        }

        /// D_v = D^v_xi * du/d_xi + D^v_eta * du/d_eta
        Vmath::Vmul(nqtot, &tangmat[0][0], 1, &diff0[0], 1, &out[0], 1);
        Vmath::Vvtvp(nqtot, &tangmat[1][0], 1, &diff1[0], 1, &out[0], 1,
                     &out[0], 1);
    }
    else
    {
        Array<OneD, Array<OneD, NekDouble>> tangmat(2);

        for (int i = 0; i < 2; ++i)
        {
            tangmat[i] = Array<OneD, NekDouble>(nqtot, 0.0);
            for (int k = 0; k < (m_geom->GetCoordim()); ++k)
            {
                Vmath::Svtvp(nqtot, df[2 * k + i][0], &direction[k * nqtot], 1,
                             &tangmat[i][0], 1, &tangmat[i][0], 1);
            }
        }

        /// D_v = D^v_xi * du/d_xi + D^v_eta * du/d_eta
        Vmath::Vmul(nqtot, &tangmat[0][0], 1, &diff0[0], 1, &out[0], 1);

        Vmath::Vvtvp(nqtot, &tangmat[1][0], 1, &diff1[0], 1, &out[0], 1,
                     &out[0], 1);
    }
}

void TriExp::v_FwdTrans(const Array<OneD, const NekDouble> &inarray,
                        Array<OneD, NekDouble> &outarray)
{
    IProductWRTBase(inarray, outarray);

    // get Mass matrix inverse
    MatrixKey masskey(StdRegions::eInvMass, DetShapeType(), *this);
    DNekScalMatSharedPtr matsys = m_matrixManager[masskey];

    // copy inarray in case inarray == outarray
    NekVector<NekDouble> in(m_ncoeffs, outarray, eCopy);
    NekVector<NekDouble> out(m_ncoeffs, outarray, eWrapper);

    out = (*matsys) * in;
}

void TriExp::v_FwdTransBndConstrained(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray)
{
    int i, j;
    int npoints[2] = {m_base[0]->GetNumPoints(), m_base[1]->GetNumPoints()};
    int nmodes[2]  = {m_base[0]->GetNumModes(), m_base[1]->GetNumModes()};

    fill(outarray.data(), outarray.data() + m_ncoeffs, 0.0);

    if (nmodes[0] == 1 && nmodes[1] == 1)
    {
        outarray[0] = inarray[0];
        return;
    }

    Array<OneD, NekDouble> physEdge[3];
    Array<OneD, NekDouble> coeffEdge[3];
    for (i = 0; i < 3; i++)
    {
        // define physEdge and add 1 so can interpolate grl10 points if
        // necessary
        physEdge[i]  = Array<OneD, NekDouble>(max(npoints[i != 0], npoints[0]));
        coeffEdge[i] = Array<OneD, NekDouble>(nmodes[i != 0]);
    }

    for (i = 0; i < npoints[0]; i++)
    {
        physEdge[0][i] = inarray[i];
    }

    // extract data in cartesian directions
    for (i = 0; i < npoints[1]; i++)
    {
        physEdge[1][i] = inarray[npoints[0] - 1 + i * npoints[0]];
        physEdge[2][i] = inarray[i * npoints[0]];
    }

    SegExpSharedPtr segexp[3];
    segexp[0] = MemoryManager<LocalRegions::SegExp>::AllocateSharedPtr(
        m_base[0]->GetBasisKey(), GetGeom2D()->GetEdge(0));

    if (m_base[1]->GetPointsType() == LibUtilities::eGaussLobattoLegendre)
    {
        for (i = 1; i < 3; i++)
        {
            segexp[i] = MemoryManager<LocalRegions::SegExp>::AllocateSharedPtr(
                m_base[i != 0]->GetBasisKey(), GetGeom2D()->GetEdge(i));
        }
    }
    else // interploate using edge 0 GLL distribution
    {
        for (i = 1; i < 3; i++)
        {
            segexp[i] = MemoryManager<LocalRegions::SegExp>::AllocateSharedPtr(
                m_base[0]->GetBasisKey(), GetGeom2D()->GetEdge(i));

            LibUtilities::Interp1D(m_base[1]->GetPointsKey(), physEdge[i],
                                   m_base[0]->GetPointsKey(), physEdge[i]);
        }
        npoints[1] = npoints[0];
    }

    Array<OneD, unsigned int> mapArray;
    Array<OneD, int> signArray;
    NekDouble sign;
    // define an orientation to get EdgeToElmtMapping from Cartesian data
    StdRegions::Orientation orient[3] = {
        StdRegions::eForwards, StdRegions::eForwards, StdRegions::eForwards};

    for (i = 0; i < 3; i++)
    {
        segexp[i]->FwdTransBndConstrained(physEdge[i], coeffEdge[i]);

        // this orient goes with the one above and so could
        // probably set both to eForwards
        GetTraceToElementMap(i, mapArray, signArray, orient[i]);
        for (j = 0; j < nmodes[i != 0]; j++)
        {
            sign                  = (NekDouble)signArray[j];
            outarray[mapArray[j]] = sign * coeffEdge[i][j];
        }
    }

    int nBoundaryDofs = NumBndryCoeffs();
    int nInteriorDofs = m_ncoeffs - nBoundaryDofs;

    if (nInteriorDofs > 0)
    {
        Array<OneD, NekDouble> tmp0(m_ncoeffs);
        Array<OneD, NekDouble> tmp1(m_ncoeffs);

        StdRegions::StdMatrixKey stdmasskey(StdRegions::eMass, DetShapeType(),
                                            *this);
        MassMatrixOp(outarray, tmp0, stdmasskey);
        IProductWRTBase(inarray, tmp1);

        Vmath::Vsub(m_ncoeffs, tmp1, 1, tmp0, 1, tmp1, 1);

        // get Mass matrix inverse (only of interior DOF)
        // use block (1,1) of the static condensed system
        // note: this block alreay contains the inverse matrix
        MatrixKey masskey(StdRegions::eMass, DetShapeType(), *this);
        DNekScalMatSharedPtr matsys =
            (m_staticCondMatrixManager[masskey])->GetBlock(1, 1);

        Array<OneD, NekDouble> rhs(nInteriorDofs);
        Array<OneD, NekDouble> result(nInteriorDofs);

        GetInteriorMap(mapArray);

        for (i = 0; i < nInteriorDofs; i++)
        {
            rhs[i] = tmp1[mapArray[i]];
        }

        Blas::Dgemv('N', nInteriorDofs, nInteriorDofs, matsys->Scale(),
                    &((matsys->GetOwnedMatrix())->GetPtr())[0], nInteriorDofs,
                    rhs.data(), 1, 0.0, result.data(), 1);

        for (i = 0; i < nInteriorDofs; i++)
        {
            outarray[mapArray[i]] = result[i];
        }
    }
}

void TriExp::v_IProductWRTBase(const Array<OneD, const NekDouble> &inarray,
                               Array<OneD, NekDouble> &outarray)
{
    IProductWRTBase_SumFac(inarray, outarray);
}

void TriExp::v_IProductWRTDerivBase(const int dir,
                                    const Array<OneD, const NekDouble> &inarray,
                                    Array<OneD, NekDouble> &outarray)
{
    IProductWRTDerivBase_SumFac(dir, inarray, outarray);
}

void TriExp::v_IProductWRTBase_SumFac(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, bool multiplybyweights)
{
    int nquad0 = m_base[0]->GetNumPoints();
    int nquad1 = m_base[1]->GetNumPoints();
    int order0 = m_base[0]->GetNumModes();

    if (multiplybyweights)
    {
        Array<OneD, NekDouble> tmp(nquad0 * nquad1 + nquad1 * order0);
        Array<OneD, NekDouble> wsp(tmp + nquad0 * nquad1);

        MultiplyByQuadratureMetric(inarray, tmp);
        IProductWRTBase_SumFacKernel(m_base[0]->GetBdata(),
                                     m_base[1]->GetBdata(), tmp, outarray, wsp);
    }
    else
    {
        Array<OneD, NekDouble> wsp(+nquad1 * order0);

        IProductWRTBase_SumFacKernel(m_base[0]->GetBdata(),
                                     m_base[1]->GetBdata(), inarray, outarray,
                                     wsp);
    }
}

void TriExp::v_IProductWRTDerivBase_SumFac(
    const int dir, const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray)
{
    int nquad0  = m_base[0]->GetNumPoints();
    int nquad1  = m_base[1]->GetNumPoints();
    int nqtot   = nquad0 * nquad1;
    int nmodes0 = m_base[0]->GetNumModes();
    int wspsize = max(max(nqtot, m_ncoeffs), nquad1 * nmodes0);

    Array<OneD, NekDouble> tmp0(4 * wspsize);
    Array<OneD, NekDouble> tmp1(tmp0 + wspsize);
    Array<OneD, NekDouble> tmp2(tmp0 + 2 * wspsize);
    Array<OneD, NekDouble> tmp3(tmp0 + 3 * wspsize);

    Array<OneD, Array<OneD, NekDouble>> tmp2D{2};
    tmp2D[0] = tmp1;
    tmp2D[1] = tmp2;

    TriExp::v_AlignVectorToCollapsedDir(dir, inarray, tmp2D);

    MultiplyByQuadratureMetric(tmp1, tmp1);
    MultiplyByQuadratureMetric(tmp2, tmp2);

    IProductWRTBase_SumFacKernel(m_base[0]->GetDbdata(), m_base[1]->GetBdata(),
                                 tmp1, tmp3, tmp0);
    IProductWRTBase_SumFacKernel(m_base[0]->GetBdata(), m_base[1]->GetDbdata(),
                                 tmp2, outarray, tmp0);
    Vmath::Vadd(m_ncoeffs, tmp3, 1, outarray, 1, outarray, 1);
}

void TriExp::v_AlignVectorToCollapsedDir(
    const int dir, const Array<OneD, const NekDouble> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    ASSERTL1((dir == 0) || (dir == 1) || (dir == 2), "Invalid direction.");
    ASSERTL1((dir == 2) ? (m_geom->GetCoordim() == 3) : true,
             "Invalid direction.");

    int nquad0  = m_base[0]->GetNumPoints();
    int nquad1  = m_base[1]->GetNumPoints();
    int nqtot   = nquad0 * nquad1;
    int nmodes0 = m_base[0]->GetNumModes();
    int wspsize = max(max(nqtot, m_ncoeffs), nquad1 * nmodes0);

    const Array<TwoD, const NekDouble> &df =
        m_metricinfo->GetDerivFactors(GetPointsKeys());

    Array<OneD, NekDouble> tmp0(wspsize);
    Array<OneD, NekDouble> tmp3(wspsize);
    Array<OneD, NekDouble> gfac0(wspsize);
    Array<OneD, NekDouble> gfac1(wspsize);

    Array<OneD, NekDouble> tmp1 = outarray[0];
    Array<OneD, NekDouble> tmp2 = outarray[1];

    const Array<OneD, const NekDouble> &z0 = m_base[0]->GetZ();
    const Array<OneD, const NekDouble> &z1 = m_base[1]->GetZ();

    // set up geometric factor: 2/(1-z1)
    for (int i = 0; i < nquad1; ++i)
    {
        gfac0[i] = 2.0 / (1 - z1[i]);
    }
    for (int i = 0; i < nquad0; ++i)
    {
        gfac1[i] = 0.5 * (1 + z0[i]);
    }

    for (int i = 0; i < nquad1; ++i)
    {
        Vmath::Smul(nquad0, gfac0[i], &inarray[0] + i * nquad0, 1,
                    &tmp0[0] + i * nquad0, 1);
    }

    for (int i = 0; i < nquad1; ++i)
    {
        Vmath::Vmul(nquad0, &gfac1[0], 1, &tmp0[0] + i * nquad0, 1,
                    &tmp1[0] + i * nquad0, 1);
    }

    if (m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
    {
        Vmath::Vmul(nqtot, &df[2 * dir][0], 1, &tmp0[0], 1, &tmp0[0], 1);
        Vmath::Vmul(nqtot, &df[2 * dir + 1][0], 1, &tmp1[0], 1, &tmp1[0], 1);
        Vmath::Vmul(nqtot, &df[2 * dir + 1][0], 1, &inarray[0], 1, &tmp2[0], 1);
    }
    else
    {
        Vmath::Smul(nqtot, df[2 * dir][0], tmp0, 1, tmp0, 1);
        Vmath::Smul(nqtot, df[2 * dir + 1][0], tmp1, 1, tmp1, 1);
        Vmath::Smul(nqtot, df[2 * dir + 1][0], inarray, 1, tmp2, 1);
    }
    Vmath::Vadd(nqtot, tmp0, 1, tmp1, 1, tmp1, 1);
}

void TriExp::v_IProductWRTDirectionalDerivBase(
    const Array<OneD, const NekDouble> &direction,
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray)
{
    IProductWRTDirectionalDerivBase_SumFac(direction, inarray, outarray);
}

/**
 * @brief Directinoal Derivative in the modal space in the dir
 * direction of varcoeffs.
 */
void TriExp::v_IProductWRTDirectionalDerivBase_SumFac(
    const Array<OneD, const NekDouble> &direction,
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray)
{
    int i;
    int shapedim = 2;
    int nquad0   = m_base[0]->GetNumPoints();
    int nquad1   = m_base[1]->GetNumPoints();
    int nqtot    = nquad0 * nquad1;
    int nmodes0  = m_base[0]->GetNumModes();
    int wspsize  = max(max(nqtot, m_ncoeffs), nquad1 * nmodes0);

    const Array<TwoD, const NekDouble> &df =
        m_metricinfo->GetDerivFactors(GetPointsKeys());

    Array<OneD, NekDouble> tmp0(6 * wspsize);
    Array<OneD, NekDouble> tmp1(tmp0 + wspsize);
    Array<OneD, NekDouble> tmp2(tmp0 + 2 * wspsize);
    Array<OneD, NekDouble> tmp3(tmp0 + 3 * wspsize);
    Array<OneD, NekDouble> gfac0(tmp0 + 4 * wspsize);
    Array<OneD, NekDouble> gfac1(tmp0 + 5 * wspsize);

    const Array<OneD, const NekDouble> &z0 = m_base[0]->GetZ();
    const Array<OneD, const NekDouble> &z1 = m_base[1]->GetZ();

    // set up geometric factor: 2/(1-z1)
    for (i = 0; i < nquad1; ++i)
    {
        gfac0[i] = 2.0 / (1 - z1[i]);
    }
    for (i = 0; i < nquad0; ++i)
    {
        gfac1[i] = 0.5 * (1 + z0[i]);
    }
    for (i = 0; i < nquad1; ++i)
    {
        Vmath::Smul(nquad0, gfac0[i], &inarray[0] + i * nquad0, 1,
                    &tmp0[0] + i * nquad0, 1);
    }
    for (i = 0; i < nquad1; ++i)
    {
        Vmath::Vmul(nquad0, &gfac1[0], 1, &tmp0[0] + i * nquad0, 1,
                    &tmp1[0] + i * nquad0, 1);
    }

    // Compute gmat \cdot e^j
    Array<OneD, Array<OneD, NekDouble>> dfdir(shapedim);
    Expansion::ComputeGmatcdotMF(df, direction, dfdir);

    Vmath::Vmul(nqtot, &dfdir[0][0], 1, &tmp0[0], 1, &tmp0[0], 1);
    Vmath::Vmul(nqtot, &dfdir[1][0], 1, &tmp1[0], 1, &tmp1[0], 1);
    Vmath::Vmul(nqtot, &dfdir[1][0], 1, &inarray[0], 1, &tmp2[0], 1);

    Vmath::Vadd(nqtot, &tmp0[0], 1, &tmp1[0], 1, &tmp1[0], 1);

    MultiplyByQuadratureMetric(tmp1, tmp1);
    MultiplyByQuadratureMetric(tmp2, tmp2);

    IProductWRTBase_SumFacKernel(m_base[0]->GetDbdata(), m_base[1]->GetBdata(),
                                 tmp1, tmp3, tmp0);
    IProductWRTBase_SumFacKernel(m_base[0]->GetBdata(), m_base[1]->GetDbdata(),
                                 tmp2, outarray, tmp0);
    Vmath::Vadd(m_ncoeffs, tmp3, 1, outarray, 1, outarray, 1);
}

void TriExp::v_NormVectorIProductWRTBase(const Array<OneD, const NekDouble> &Fx,
                                         const Array<OneD, const NekDouble> &Fy,
                                         const Array<OneD, const NekDouble> &Fz,
                                         Array<OneD, NekDouble> &outarray)
{
    int nq = m_base[0]->GetNumPoints() * m_base[1]->GetNumPoints();
    Array<OneD, NekDouble> Fn(nq);

    const Array<OneD, const Array<OneD, NekDouble>> &normals =
        GetLeftAdjacentElementExp()->GetTraceNormal(
            GetLeftAdjacentElementTrace());

    if (m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
    {
        Vmath::Vvtvvtp(nq, &normals[0][0], 1, &Fx[0], 1, &normals[1][0], 1,
                       &Fy[0], 1, &Fn[0], 1);
        Vmath::Vvtvp(nq, &normals[2][0], 1, &Fz[0], 1, &Fn[0], 1, &Fn[0], 1);
    }
    else
    {
        Vmath::Svtsvtp(nq, normals[0][0], &Fx[0], 1, normals[1][0], &Fy[0], 1,
                       &Fn[0], 1);
        Vmath::Svtvp(nq, normals[2][0], &Fz[0], 1, &Fn[0], 1, &Fn[0], 1);
    }

    IProductWRTBase(Fn, outarray);
}

void TriExp::v_NormVectorIProductWRTBase(
    const Array<OneD, const Array<OneD, NekDouble>> &Fvec,
    Array<OneD, NekDouble> &outarray)
{
    NormVectorIProductWRTBase(Fvec[0], Fvec[1], Fvec[2], outarray);
}

StdRegions::StdExpansionSharedPtr TriExp::v_GetStdExp(void) const
{

    return MemoryManager<StdRegions::StdTriExp>::AllocateSharedPtr(
        m_base[0]->GetBasisKey(), m_base[1]->GetBasisKey());
}

StdRegions::StdExpansionSharedPtr TriExp::v_GetLinStdExp(void) const
{
    LibUtilities::BasisKey bkey0(m_base[0]->GetBasisType(), 2,
                                 m_base[0]->GetPointsKey());
    LibUtilities::BasisKey bkey1(m_base[1]->GetBasisType(), 2,
                                 m_base[1]->GetPointsKey());

    return MemoryManager<StdRegions::StdTriExp>::AllocateSharedPtr(bkey0,
                                                                   bkey1);
}

void TriExp::v_GetCoord(const Array<OneD, const NekDouble> &Lcoords,
                        Array<OneD, NekDouble> &coords)
{
    int i;

    ASSERTL1(Lcoords[0] >= -1.0 && Lcoords[1] <= 1.0 && Lcoords[1] >= -1.0 &&
                 Lcoords[1] <= 1.0,
             "Local coordinates are not in region [-1,1]");

    m_geom->FillGeom();

    for (i = 0; i < m_geom->GetCoordim(); ++i)
    {
        coords[i] = m_geom->GetCoord(i, Lcoords);
    }
}

void TriExp::v_GetCoords(Array<OneD, NekDouble> &coords_0,
                         Array<OneD, NekDouble> &coords_1,
                         Array<OneD, NekDouble> &coords_2)
{
    Expansion::v_GetCoords(coords_0, coords_1, coords_2);
}

/**
 * Given the local cartesian coordinate \a Lcoord evaluate the
 * value of physvals at this point by calling through to the
 * StdExpansion method
 */
NekDouble TriExp::v_StdPhysEvaluate(
    const Array<OneD, const NekDouble> &Lcoord,
    const Array<OneD, const NekDouble> &physvals)
{
    // Evaluate point in local (eta) coordinates.
    return StdExpansion2D::v_PhysEvaluate(Lcoord, physvals);
}

NekDouble TriExp::v_PhysEvaluate(const Array<OneD, const NekDouble> &coord,
                                 const Array<OneD, const NekDouble> &physvals)
{
    Array<OneD, NekDouble> Lcoord = Array<OneD, NekDouble>(2);

    ASSERTL0(m_geom, "m_geom not defined");
    m_geom->GetLocCoords(coord, Lcoord);

    return StdExpansion2D::v_PhysEvaluate(Lcoord, physvals);
}

NekDouble TriExp::v_PhysEvalFirstDeriv(
    const Array<OneD, NekDouble> &coord,
    const Array<OneD, const NekDouble> &inarray,
    std::array<NekDouble, 3> &firstOrderDerivs)
{
    Array<OneD, NekDouble> Lcoord(2);
    ASSERTL0(m_geom, "m_geom not defined");
    m_geom->GetLocCoords(coord, Lcoord);
    return StdTriExp::v_PhysEvalFirstDeriv(Lcoord, inarray, firstOrderDerivs);
}

void TriExp::v_GetTracePhysVals(
    const int edge, const StdRegions::StdExpansionSharedPtr &EdgeExp,
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, StdRegions::Orientation orient)
{
    int nquad0 = m_base[0]->GetNumPoints();
    int nquad1 = m_base[1]->GetNumPoints();
    int nt     = 0;
    // Extract in Cartesian direction because we have to deal with
    // e.g. Gauss-Radau points.
    switch (edge)
    {
        case 0:
            Vmath::Vcopy(nquad0, &(inarray[0]), 1, &(outarray[0]), 1);
            nt = nquad0;
            break;
        case 1:
            Vmath::Vcopy(nquad1, &(inarray[0]) + (nquad0 - 1), nquad0,
                         &(outarray[0]), 1);
            nt = nquad1;
            break;
        case 2:
            Vmath::Vcopy(nquad1, &(inarray[0]), nquad0, &(outarray[0]), 1);
            nt = nquad1;
            break;
        default:
            ASSERTL0(false, "edge value (< 3) is out of range");
            break;
    }

    ASSERTL1(EdgeExp->GetBasis(0)->GetPointsType() ==
                 LibUtilities::eGaussLobattoLegendre,
             "Edge expansion should be GLL");

    // Interpolate if required
    if (m_base[edge ? 1 : 0]->GetPointsKey() !=
        EdgeExp->GetBasis(0)->GetPointsKey())
    {
        Array<OneD, NekDouble> outtmp(max(nquad0, nquad1));

        Vmath::Vcopy(nt, outarray, 1, outtmp, 1);

        LibUtilities::Interp1D(m_base[edge ? 1 : 0]->GetPointsKey(), outtmp,
                               EdgeExp->GetBasis(0)->GetPointsKey(), outarray);
    }

    if (orient == StdRegions::eNoOrientation)
    {
        orient = GetTraceOrient(edge);
    }

    // Reverse data if necessary
    if (orient == StdRegions::eBackwards)
    {
        Vmath::Reverse(EdgeExp->GetNumPoints(0), &outarray[0], 1, &outarray[0],
                       1);
    }
}

void TriExp::v_GetTraceQFactors(
    [[maybe_unused]] const int edge,
    [[maybe_unused]] Array<OneD, NekDouble> &outarray)
{
    ASSERTL0(false, "Routine not implemented for triangular elements");
}

void TriExp::v_GetTracePhysMap(const int edge, Array<OneD, int> &outarray)
{
    int nquad0 = m_base[0]->GetNumPoints();
    int nquad1 = m_base[1]->GetNumPoints();

    // Get points in Cartesian orientation
    switch (edge)
    {
        case 0:
            outarray = Array<OneD, int>(nquad0);
            for (int i = 0; i < nquad0; ++i)
            {
                outarray[i] = i;
            }
            break;
        case 1:
            outarray = Array<OneD, int>(nquad1);
            for (int i = 0; i < nquad1; ++i)
            {
                outarray[i] = (nquad0 - 1) + i * nquad0;
            }
            break;
        case 2:
            outarray = Array<OneD, int>(nquad1);
            for (int i = 0; i < nquad1; ++i)
            {
                outarray[i] = i * nquad0;
            }
            break;
        default:
            ASSERTL0(false, "edge value (< 3) is out of range");
            break;
    }
}

void TriExp::v_ComputeTraceNormal(const int edge)
{
    int i;
    const SpatialDomains::GeomFactorsSharedPtr &geomFactors =
        GetGeom()->GetMetricInfo();

    LibUtilities::PointsKeyVector ptsKeys = GetPointsKeys();
    for (i = 0; i < ptsKeys.size(); ++i)
    {
        // Need at least 2 points for computing normals
        if (ptsKeys[i].GetNumPoints() == 1)
        {
            LibUtilities::PointsKey pKey(2, ptsKeys[i].GetPointsType());
            ptsKeys[i] = pKey;
        }
    }

    const SpatialDomains::GeomType type = geomFactors->GetGtype();
    const Array<TwoD, const NekDouble> &df =
        geomFactors->GetDerivFactors(ptsKeys);
    const Array<OneD, const NekDouble> &jac = geomFactors->GetJac(ptsKeys);

    // The points of normals should follow trace basis, not local basis.
    LibUtilities::BasisKey tobasis = GetTraceBasisKey(edge);

    int nqe = tobasis.GetNumPoints();
    int dim = GetCoordim();

    m_traceNormals[edge] = Array<OneD, Array<OneD, NekDouble>>(dim);
    Array<OneD, Array<OneD, NekDouble>> &normal = m_traceNormals[edge];
    for (i = 0; i < dim; ++i)
    {
        normal[i] = Array<OneD, NekDouble>(nqe);
    }

    size_t nqb                     = nqe;
    size_t nbnd                    = edge;
    m_elmtBndNormDirElmtLen[nbnd]  = Array<OneD, NekDouble>{nqb, 0.0};
    Array<OneD, NekDouble> &length = m_elmtBndNormDirElmtLen[nbnd];

    // Regular geometry case
    if ((type == SpatialDomains::eRegular) ||
        (type == SpatialDomains::eMovingRegular))
    {
        NekDouble fac;
        // Set up normals
        switch (edge)
        {
            case 0:
                for (i = 0; i < GetCoordim(); ++i)
                {
                    Vmath::Fill(nqe, -df[2 * i + 1][0], normal[i], 1);
                }
                break;
            case 1:
                for (i = 0; i < GetCoordim(); ++i)
                {
                    Vmath::Fill(nqe, df[2 * i + 1][0] + df[2 * i][0], normal[i],
                                1);
                }
                break;
            case 2:
                for (i = 0; i < GetCoordim(); ++i)
                {
                    Vmath::Fill(nqe, -df[2 * i][0], normal[i], 1);
                }
                break;
            default:
                ASSERTL0(false, "Edge is out of range (edge < 3)");
        }

        // normalise
        fac = 0.0;
        for (i = 0; i < GetCoordim(); ++i)
        {
            fac += normal[i][0] * normal[i][0];
        }
        fac = 1.0 / sqrt(fac);

        Vmath::Fill(nqb, fac, length, 1);

        for (i = 0; i < GetCoordim(); ++i)
        {
            Vmath::Smul(nqe, fac, normal[i], 1, normal[i], 1);
        }
    }
    else // Set up deformed normals
    {
        int j;

        int nquad0 = ptsKeys[0].GetNumPoints();
        int nquad1 = ptsKeys[1].GetNumPoints();

        LibUtilities::PointsKey from_key;

        Array<OneD, NekDouble> normals(GetCoordim() * max(nquad0, nquad1), 0.0);
        Array<OneD, NekDouble> edgejac(GetCoordim() * max(nquad0, nquad1), 0.0);

        // Extract Jacobian along edges and recover local
        // derivates (dx/dr) for polynomial interpolation by
        // multiplying m_gmat by jacobian
        switch (edge)
        {
            case 0:
                for (j = 0; j < nquad0; ++j)
                {
                    edgejac[j] = jac[j];
                    for (i = 0; i < GetCoordim(); ++i)
                    {
                        normals[i * nquad0 + j] =
                            -df[2 * i + 1][j] * edgejac[j];
                    }
                }
                from_key = ptsKeys[0];
                break;
            case 1:
                for (j = 0; j < nquad1; ++j)
                {
                    edgejac[j] = jac[nquad0 * j + nquad0 - 1];
                    for (i = 0; i < GetCoordim(); ++i)
                    {
                        normals[i * nquad1 + j] =
                            (df[2 * i][nquad0 * j + nquad0 - 1] +
                             df[2 * i + 1][nquad0 * j + nquad0 - 1]) *
                            edgejac[j];
                    }
                }
                from_key = ptsKeys[1];
                break;
            case 2:
                for (j = 0; j < nquad1; ++j)
                {
                    edgejac[j] = jac[nquad0 * j];
                    for (i = 0; i < GetCoordim(); ++i)
                    {
                        normals[i * nquad1 + j] =
                            -df[2 * i][nquad0 * j] * edgejac[j];
                    }
                }
                from_key = ptsKeys[1];
                break;
            default:
                ASSERTL0(false, "edge is out of range (edge < 3)");
        }

        int nq = from_key.GetNumPoints();
        Array<OneD, NekDouble> work(nqe, 0.0);

        // interpolate Jacobian and invert
        LibUtilities::Interp1D(from_key, jac, tobasis.GetPointsKey(), work);
        Vmath::Sdiv(nqe, 1.0, &work[0], 1, &work[0], 1);

        // interpolate
        for (i = 0; i < GetCoordim(); ++i)
        {
            LibUtilities::Interp1D(from_key, &normals[i * nq],
                                   tobasis.GetPointsKey(), &normal[i][0]);
            Vmath::Vmul(nqe, work, 1, normal[i], 1, normal[i], 1);
        }

        // normalise normal vectors
        Vmath::Zero(nqe, work, 1);
        for (i = 0; i < GetCoordim(); ++i)
        {
            Vmath::Vvtvp(nqe, normal[i], 1, normal[i], 1, work, 1, work, 1);
        }

        Vmath::Vsqrt(nqe, work, 1, work, 1);
        Vmath::Sdiv(nqe, 1.0, work, 1, work, 1);

        Vmath::Vcopy(nqb, work, 1, length, 1);

        for (i = 0; i < GetCoordim(); ++i)
        {
            Vmath::Vmul(nqe, normal[i], 1, work, 1, normal[i], 1);
        }
    }

    if (GetGeom()->GetEorient(edge) == StdRegions::eBackwards)
    {
        for (i = 0; i < GetCoordim(); ++i)
        {
            if (geomFactors->GetGtype() == SpatialDomains::eDeformed)
            {
                Vmath::Reverse(nqe, normal[i], 1, normal[i], 1);
            }
        }
    }
}

void TriExp::v_ExtractDataToCoeffs(
    const NekDouble *data, const std::vector<unsigned int> &nummodes,
    const int mode_offset, NekDouble *coeffs,
    [[maybe_unused]] std::vector<LibUtilities::BasisType> &fromType)
{
    int data_order0 = nummodes[mode_offset];
    int fillorder0  = min(m_base[0]->GetNumModes(), data_order0);
    int data_order1 = nummodes[mode_offset + 1];
    int order1      = m_base[1]->GetNumModes();
    int fillorder1  = min(order1, data_order1);

    switch (m_base[0]->GetBasisType())
    {
        case LibUtilities::eModified_A:
        case LibUtilities::eOrtho_A:
        {
            int i;
            int cnt  = 0;
            int cnt1 = 0;

            ASSERTL1(m_base[1]->GetBasisType() == LibUtilities::eModified_B ||
                         m_base[1]->GetBasisType() == LibUtilities::eOrtho_B,
                     "Extraction routine not set up for this basis");

            Vmath::Zero(m_ncoeffs, coeffs, 1);
            for (i = 0; i < fillorder0; ++i)
            {
                Vmath::Vcopy(fillorder1 - i, &data[cnt], 1, &coeffs[cnt1], 1);
                cnt += data_order1 - i;
                cnt1 += order1 - i;
            }
        }
        break;
        default:
            ASSERTL0(false, "basis is either not set up or not hierarchicial");
    }
}

StdRegions::Orientation TriExp::v_GetTraceOrient(int edge)
{
    return GetGeom2D()->GetEorient(edge);
}

DNekMatSharedPtr TriExp::v_GenMatrix(const StdRegions::StdMatrixKey &mkey)
{
    DNekMatSharedPtr returnval;
    switch (mkey.GetMatrixType())
    {
        case StdRegions::eHybridDGHelmholtz:
        case StdRegions::eHybridDGLamToU:
        case StdRegions::eHybridDGLamToQ0:
        case StdRegions::eHybridDGLamToQ1:
        case StdRegions::eHybridDGLamToQ2:
        case StdRegions::eHybridDGHelmBndLam:
        case StdRegions::eInvLaplacianWithUnityMean:
            returnval = Expansion2D::v_GenMatrix(mkey);
            break;
        default:
            returnval = StdTriExp::v_GenMatrix(mkey);
            break;
    }

    return returnval;
}

DNekMatSharedPtr TriExp::v_CreateStdMatrix(const StdRegions::StdMatrixKey &mkey)
{
    LibUtilities::BasisKey bkey0 = m_base[0]->GetBasisKey();
    LibUtilities::BasisKey bkey1 = m_base[1]->GetBasisKey();
    StdRegions::StdTriExpSharedPtr tmp =
        MemoryManager<StdTriExp>::AllocateSharedPtr(bkey0, bkey1);

    return tmp->GetStdMatrix(mkey);
}

DNekScalMatSharedPtr TriExp::v_GetLocMatrix(const MatrixKey &mkey)
{
    return m_matrixManager[mkey];
}

void TriExp::v_DropLocMatrix(const MatrixKey &mkey)
{
    m_matrixManager.DeleteObject(mkey);
}

DNekScalBlkMatSharedPtr TriExp::v_GetLocStaticCondMatrix(const MatrixKey &mkey)
{
    return m_staticCondMatrixManager[mkey];
}

void TriExp::v_DropLocStaticCondMatrix(const MatrixKey &mkey)
{
    m_staticCondMatrixManager.DeleteObject(mkey);
}

void TriExp::v_MassMatrixOp(const Array<OneD, const NekDouble> &inarray,
                            Array<OneD, NekDouble> &outarray,
                            const StdRegions::StdMatrixKey &mkey)
{
    StdExpansion::MassMatrixOp_MatFree(inarray, outarray, mkey);
}

void TriExp::v_LaplacianMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                 Array<OneD, NekDouble> &outarray,
                                 const StdRegions::StdMatrixKey &mkey)
{
    TriExp::LaplacianMatrixOp_MatFree(inarray, outarray, mkey);
}

void TriExp::v_LaplacianMatrixOp(const int k1, const int k2,
                                 const Array<OneD, const NekDouble> &inarray,
                                 Array<OneD, NekDouble> &outarray,
                                 const StdRegions::StdMatrixKey &mkey)
{
    StdExpansion::LaplacianMatrixOp_MatFree(k1, k2, inarray, outarray, mkey);
}

void TriExp::v_WeakDerivMatrixOp(const int i,
                                 const Array<OneD, const NekDouble> &inarray,
                                 Array<OneD, NekDouble> &outarray,
                                 const StdRegions::StdMatrixKey &mkey)
{
    StdExpansion::WeakDerivMatrixOp_MatFree(i, inarray, outarray, mkey);
}

void TriExp::v_WeakDirectionalDerivMatrixOp(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const StdRegions::StdMatrixKey &mkey)
{
    StdExpansion::WeakDirectionalDerivMatrixOp_MatFree(inarray, outarray, mkey);
}

void TriExp::v_MassLevelCurvatureMatrixOp(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const StdRegions::StdMatrixKey &mkey)
{
    StdExpansion::MassLevelCurvatureMatrixOp_MatFree(inarray, outarray, mkey);
}

void TriExp::v_HelmholtzMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                 Array<OneD, NekDouble> &outarray,
                                 const StdRegions::StdMatrixKey &mkey)
{
    TriExp::HelmholtzMatrixOp_MatFree(inarray, outarray, mkey);
}

void TriExp::v_LaplacianMatrixOp_MatFree_Kernel(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, Array<OneD, NekDouble> &wsp)
{
    if (m_metrics.count(eMetricLaplacian00) == 0)
    {
        ComputeLaplacianMetric();
    }

    int nquad0  = m_base[0]->GetNumPoints();
    int nquad1  = m_base[1]->GetNumPoints();
    int nqtot   = nquad0 * nquad1;
    int nmodes0 = m_base[0]->GetNumModes();
    int nmodes1 = m_base[1]->GetNumModes();
    int wspsize =
        max(max(max(nqtot, m_ncoeffs), nquad1 * nmodes0), nquad0 * nmodes1);

    ASSERTL1(wsp.size() >= 3 * wspsize, "Workspace is of insufficient size.");

    const Array<OneD, const NekDouble> &base0  = m_base[0]->GetBdata();
    const Array<OneD, const NekDouble> &base1  = m_base[1]->GetBdata();
    const Array<OneD, const NekDouble> &dbase0 = m_base[0]->GetDbdata();
    const Array<OneD, const NekDouble> &dbase1 = m_base[1]->GetDbdata();
    const Array<OneD, const NekDouble> &metric00 =
        m_metrics[eMetricLaplacian00];
    const Array<OneD, const NekDouble> &metric01 =
        m_metrics[eMetricLaplacian01];
    const Array<OneD, const NekDouble> &metric11 =
        m_metrics[eMetricLaplacian11];

    // Allocate temporary storage
    Array<OneD, NekDouble> wsp0(wsp);
    Array<OneD, NekDouble> wsp1(wsp + wspsize);
    Array<OneD, NekDouble> wsp2(wsp + 2 * wspsize);

    StdExpansion2D::PhysTensorDeriv(inarray, wsp1, wsp2);

    // wsp0 = k = g0 * wsp1 + g1 * wsp2 = g0 * du_dxi1 + g1 * du_dxi2
    // wsp2 = l = g1 * wsp1 + g2 * wsp2 = g0 * du_dxi1 + g1 * du_dxi2
    // where g0, g1 and g2 are the metric terms set up in the GeomFactors class
    // especially for this purpose
    Vmath::Vvtvvtp(nqtot, &metric00[0], 1, &wsp1[0], 1, &metric01[0], 1,
                   &wsp2[0], 1, &wsp0[0], 1);
    Vmath::Vvtvvtp(nqtot, &metric01[0], 1, &wsp1[0], 1, &metric11[0], 1,
                   &wsp2[0], 1, &wsp2[0], 1);

    // outarray = m = (D_xi1 * B)^T * k
    // wsp1     = n = (D_xi2 * B)^T * l
    IProductWRTBase_SumFacKernel(dbase0, base1, wsp0, outarray, wsp1);
    IProductWRTBase_SumFacKernel(base0, dbase1, wsp2, wsp1, wsp0);

    // outarray = outarray + wsp1
    //          = L * u_hat
    Vmath::Vadd(m_ncoeffs, wsp1.data(), 1, outarray.data(), 1, outarray.data(),
                1);
}

void TriExp::v_ComputeLaplacianMetric()
{
    if (m_metrics.count(eMetricQuadrature) == 0)
    {
        ComputeQuadratureMetric();
    }

    unsigned int i, j;
    const SpatialDomains::GeomType type = m_metricinfo->GetGtype();
    const unsigned int nqtot            = GetTotPoints();
    const unsigned int dim              = 2;
    const MetricType m[3][3]            = {
        {eMetricLaplacian00, eMetricLaplacian01, eMetricLaplacian02},
        {eMetricLaplacian01, eMetricLaplacian11, eMetricLaplacian12},
        {eMetricLaplacian02, eMetricLaplacian12, eMetricLaplacian22}};

    Array<OneD, NekDouble> dEta_dXi[2] = {Array<OneD, NekDouble>(nqtot, 1.0),
                                          Array<OneD, NekDouble>(nqtot, 1.0)};

    for (i = 0; i < dim; ++i)
    {
        for (j = i; j < dim; ++j)
        {
            m_metrics[m[i][j]] = Array<OneD, NekDouble>(nqtot);
        }
    }

    const Array<OneD, const NekDouble> &z0 = m_base[0]->GetZ();
    const Array<OneD, const NekDouble> &z1 = m_base[1]->GetZ();
    const unsigned int nquad0              = m_base[0]->GetNumPoints();
    const unsigned int nquad1              = m_base[1]->GetNumPoints();
    const Array<TwoD, const NekDouble> &df =
        m_metricinfo->GetDerivFactors(GetPointsKeys());

    for (i = 0; i < nquad1; i++)
    {
        Blas::Dscal(nquad0, 2.0 / (1 - z1[i]), &dEta_dXi[0][0] + i * nquad0, 1);
        Blas::Dscal(nquad0, 2.0 / (1 - z1[i]), &dEta_dXi[1][0] + i * nquad0, 1);
    }
    for (i = 0; i < nquad0; i++)
    {
        Blas::Dscal(nquad1, 0.5 * (1 + z0[i]), &dEta_dXi[1][0] + i, nquad0);
    }

    Array<OneD, NekDouble> tmp(nqtot);
    if ((type == SpatialDomains::eRegular ||
         type == SpatialDomains::eMovingRegular))
    {
        Vmath::Smul(nqtot, df[0][0], &dEta_dXi[0][0], 1, &tmp[0], 1);
        Vmath::Svtvp(nqtot, df[1][0], &dEta_dXi[1][0], 1, &tmp[0], 1, &tmp[0],
                     1);

        Vmath::Vmul(nqtot, &tmp[0], 1, &tmp[0], 1,
                    &m_metrics[eMetricLaplacian00][0], 1);
        Vmath::Smul(nqtot, df[1][0], &tmp[0], 1,
                    &m_metrics[eMetricLaplacian01][0], 1);

        Vmath::Smul(nqtot, df[2][0], &dEta_dXi[0][0], 1, &tmp[0], 1);
        Vmath::Svtvp(nqtot, df[3][0], &dEta_dXi[1][0], 1, &tmp[0], 1, &tmp[0],
                     1);

        Vmath::Vvtvp(nqtot, &tmp[0], 1, &tmp[0], 1,
                     &m_metrics[eMetricLaplacian00][0], 1,
                     &m_metrics[eMetricLaplacian00][0], 1);
        Vmath::Svtvp(nqtot, df[3][0], &tmp[0], 1,
                     &m_metrics[eMetricLaplacian01][0], 1,
                     &m_metrics[eMetricLaplacian01][0], 1);

        if (GetCoordim() == 3)
        {
            Vmath::Smul(nqtot, df[4][0], &dEta_dXi[0][0], 1, &tmp[0], 1);
            Vmath::Svtvp(nqtot, df[5][0], &dEta_dXi[1][0], 1, &tmp[0], 1,
                         &tmp[0], 1);

            Vmath::Vvtvp(nqtot, &tmp[0], 1, &tmp[0], 1,
                         &m_metrics[eMetricLaplacian00][0], 1,
                         &m_metrics[eMetricLaplacian00][0], 1);
            Vmath::Svtvp(nqtot, df[5][0], &tmp[0], 1,
                         &m_metrics[eMetricLaplacian01][0], 1,
                         &m_metrics[eMetricLaplacian01][0], 1);
        }

        NekDouble g2 = df[1][0] * df[1][0] + df[3][0] * df[3][0];
        if (GetCoordim() == 3)
        {
            g2 += df[5][0] * df[5][0];
        }
        Vmath::Fill(nqtot, g2, &m_metrics[eMetricLaplacian11][0], 1);
    }
    else
    {

        Vmath::Vmul(nqtot, &df[0][0], 1, &dEta_dXi[0][0], 1, &tmp[0], 1);
        Vmath::Vvtvp(nqtot, &df[1][0], 1, &dEta_dXi[1][0], 1, &tmp[0], 1,
                     &tmp[0], 1);

        Vmath::Vmul(nqtot, &tmp[0], 1, &tmp[0], 1,
                    &m_metrics[eMetricLaplacian00][0], 1);
        Vmath::Vmul(nqtot, &df[1][0], 1, &tmp[0], 1,
                    &m_metrics[eMetricLaplacian01][0], 1);
        Vmath::Vmul(nqtot, &df[1][0], 1, &df[1][0], 1,
                    &m_metrics[eMetricLaplacian11][0], 1);

        Vmath::Vmul(nqtot, &df[2][0], 1, &dEta_dXi[0][0], 1, &tmp[0], 1);
        Vmath::Vvtvp(nqtot, &df[3][0], 1, &dEta_dXi[1][0], 1, &tmp[0], 1,
                     &tmp[0], 1);

        Vmath::Vvtvp(nqtot, &tmp[0], 1, &tmp[0], 1,
                     &m_metrics[eMetricLaplacian00][0], 1,
                     &m_metrics[eMetricLaplacian00][0], 1);
        Vmath::Vvtvp(nqtot, &df[3][0], 1, &tmp[0], 1,
                     &m_metrics[eMetricLaplacian01][0], 1,
                     &m_metrics[eMetricLaplacian01][0], 1);
        Vmath::Vvtvp(nqtot, &df[3][0], 1, &df[3][0], 1,
                     &m_metrics[eMetricLaplacian11][0], 1,
                     &m_metrics[eMetricLaplacian11][0], 1);

        if (GetCoordim() == 3)
        {
            Vmath::Vmul(nqtot, &df[4][0], 1, &dEta_dXi[0][0], 1, &tmp[0], 1);
            Vmath::Vvtvp(nqtot, &df[5][0], 1, &dEta_dXi[1][0], 1, &tmp[0], 1,
                         &tmp[0], 1);

            Vmath::Vvtvp(nqtot, &tmp[0], 1, &tmp[0], 1,
                         &m_metrics[eMetricLaplacian00][0], 1,
                         &m_metrics[eMetricLaplacian00][0], 1);
            Vmath::Vvtvp(nqtot, &df[5][0], 1, &tmp[0], 1,
                         &m_metrics[eMetricLaplacian01][0], 1,
                         &m_metrics[eMetricLaplacian01][0], 1);
            Vmath::Vvtvp(nqtot, &df[5][0], 1, &df[5][0], 1,
                         &m_metrics[eMetricLaplacian11][0], 1,
                         &m_metrics[eMetricLaplacian11][0], 1);
        }
    }

    for (unsigned int i = 0; i < dim; ++i)
    {
        for (unsigned int j = i; j < dim; ++j)
        {
            MultiplyByQuadratureMetric(m_metrics[m[i][j]], m_metrics[m[i][j]]);
        }
    }
}

/**
 * Function is used to compute exactly the advective numerical flux on
 * theinterface of two elements with different expansions, hence an
 * appropriate number of Gauss points has to be used. The number of
 * Gauss points has to be equal to the number used by the highest
 * polynomial degree of the two adjacent elements. Furthermore, this
 * function is used to compute the sensor value in each element.
 *
 * @param   numMin     Is the reduced polynomial order
 * @param   inarray    Input array of coefficients
 * @param   dumpVar    Output array of reduced coefficients.
 */
void TriExp::v_ReduceOrderCoeffs(int numMin,
                                 const Array<OneD, const NekDouble> &inarray,
                                 Array<OneD, NekDouble> &outarray)
{
    int n_coeffs = inarray.size();
    int nquad0   = m_base[0]->GetNumPoints();
    int nquad1   = m_base[1]->GetNumPoints();
    int nqtot    = nquad0 * nquad1;
    int nmodes0  = m_base[0]->GetNumModes();
    int nmodes1  = m_base[1]->GetNumModes();
    int numMin2  = nmodes0, i;

    Array<OneD, NekDouble> coeff(n_coeffs, 0.0);
    Array<OneD, NekDouble> phys_tmp(nqtot, 0.0);
    Array<OneD, NekDouble> tmp, tmp2;

    const LibUtilities::PointsKey Pkey0 = m_base[0]->GetPointsKey();
    const LibUtilities::PointsKey Pkey1 = m_base[1]->GetPointsKey();

    LibUtilities::BasisKey b0(m_base[0]->GetBasisType(), nmodes0, Pkey0);
    LibUtilities::BasisKey b1(m_base[1]->GetBasisType(), nmodes1, Pkey1);
    LibUtilities::BasisKey bortho0(LibUtilities::eOrtho_A, nmodes0, Pkey0);
    LibUtilities::BasisKey bortho1(LibUtilities::eOrtho_B, nmodes1, Pkey1);

    // Check if it is also possible to use the same InterCoeff routine
    // which is also used for Quadrilateral and Hexagonal shaped
    // elements

    // For now, set up the used basis on the standard element to
    // calculate the phys values, set up the orthogonal basis to do a
    // forward transform, to obtain the coefficients in orthogonal
    // coefficient space
    StdRegions::StdTriExpSharedPtr m_OrthoTriExp;
    StdRegions::StdTriExpSharedPtr m_TriExp;

    m_TriExp = MemoryManager<StdRegions::StdTriExp>::AllocateSharedPtr(b0, b1);
    m_OrthoTriExp = MemoryManager<StdRegions::StdTriExp>::AllocateSharedPtr(
        bortho0, bortho1);

    m_TriExp->BwdTrans(inarray, phys_tmp);
    m_OrthoTriExp->FwdTrans(phys_tmp, coeff);

    for (i = 0; i < n_coeffs; i++)
    {
        if (i == numMin)
        {
            coeff[i] = 0.0;
            numMin += numMin2 - 1;
            numMin2 -= 1.0;
        }
    }

    m_OrthoTriExp->BwdTrans(coeff, phys_tmp);
    m_TriExp->FwdTrans(phys_tmp, outarray);
}

void TriExp::v_SVVLaplacianFilter(Array<OneD, NekDouble> &array,
                                  const StdRegions::StdMatrixKey &mkey)
{
    int nq = GetTotPoints();

    // Calculate sqrt of the Jacobian
    Array<OneD, const NekDouble> jac = m_metricinfo->GetJac(GetPointsKeys());
    Array<OneD, NekDouble> sqrt_jac(nq);
    if (m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
    {
        Vmath::Vsqrt(nq, jac, 1, sqrt_jac, 1);
    }
    else
    {
        Vmath::Fill(nq, sqrt(jac[0]), sqrt_jac, 1);
    }

    // Multiply array by sqrt(Jac)
    Vmath::Vmul(nq, sqrt_jac, 1, array, 1, array, 1);

    // Apply std region filter
    StdTriExp::v_SVVLaplacianFilter(array, mkey);

    // Divide by sqrt(Jac)
    Vmath::Vdiv(nq, array, 1, sqrt_jac, 1, array, 1);
}

/** @brief: This method gets all of the factors which are
    required as part of the Gradient Jump Penalty
    stabilisation and involves the product of the normal and
    geometric factors along the element trace.
*/
void TriExp::v_NormalTraceDerivFactors(
    Array<OneD, Array<OneD, NekDouble>> &d0factors,
    Array<OneD, Array<OneD, NekDouble>> &d1factors,
    [[maybe_unused]] Array<OneD, Array<OneD, NekDouble>> &d2factors)
{
    int nquad0 = GetNumPoints(0);
    int nquad1 = GetNumPoints(1);

    const Array<TwoD, const NekDouble> &df =
        m_metricinfo->GetDerivFactors(GetPointsKeys());

    if (d0factors.size() != 3)
    {
        d0factors = Array<OneD, Array<OneD, NekDouble>>(3);
        d1factors = Array<OneD, Array<OneD, NekDouble>>(3);
    }

    if (d0factors[0].size() != nquad0)
    {
        d0factors[0] = Array<OneD, NekDouble>(nquad0);
        d1factors[0] = Array<OneD, NekDouble>(nquad0);
    }

    if (d0factors[1].size() != nquad1)
    {
        d0factors[1] = Array<OneD, NekDouble>(nquad1);
        d0factors[2] = Array<OneD, NekDouble>(nquad1);
        d1factors[1] = Array<OneD, NekDouble>(nquad1);
        d1factors[2] = Array<OneD, NekDouble>(nquad1);
    }

    // Outwards normals
    const Array<OneD, const Array<OneD, NekDouble>> &normal_0 =
        GetTraceNormal(0);
    const Array<OneD, const Array<OneD, NekDouble>> &normal_1 =
        GetTraceNormal(1);
    const Array<OneD, const Array<OneD, NekDouble>> &normal_2 =
        GetTraceNormal(2);

    int ncoords = normal_0.size();

    if (m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
    {

        // d xi_2/dx n_x
        for (int i = 0; i < nquad0; ++i)
        {
            d1factors[0][i] = df[1][i] * normal_0[0][i];
        }

        // d xi_1/dx n_x
        for (int i = 0; i < nquad1; ++i)
        {
            d0factors[1][i] = df[0][(i + 1) * nquad0 - 1] * normal_1[0][i];
            d0factors[2][i] = df[0][i * nquad0] * normal_2[0][i];
        }

        for (int n = 1; n < ncoords; ++n)
        {
            // d xi_2/dy n_y
            // needs checking for 3D coords
            for (int i = 0; i < nquad0; ++i)
            {
                d1factors[0][i] += df[2 * n + 1][i] * normal_0[n][i];
            }

            // d xi_1/dy n_y
            // needs checking for 3D coords
            for (int i = 0; i < nquad1; ++i)
            {
                d0factors[1][i] +=
                    df[2 * n][(i + 1) * nquad0 - 1] * normal_1[n][i];
                d0factors[2][i] += df[2 * n][i * nquad0] * normal_2[n][i];
            }
        }

        // d0 factors
        // d xi_1/dx n_x
        for (int i = 0; i < nquad0; ++i)
        {
            d0factors[0][i] = df[0][i] * normal_0[0][i];
        }

        // d xi_2/dx n_x
        for (int i = 0; i < nquad1; ++i)
        {
            d1factors[1][i] = df[1][(i + 1) * nquad0 - 1] * normal_1[0][i];
            d1factors[2][i] = df[1][i * nquad0] * normal_2[0][i];
        }

        for (int n = 1; n < ncoords; ++n)
        {
            // d xi_1/dy n_y
            // needs checking for 3D coords
            for (int i = 0; i < nquad0; ++i)
            {
                d0factors[0][i] += df[2 * n][i] * normal_0[n][i];
            }

            // d xi_2/dy n_y
            // needs checking for 3D coords
            for (int i = 0; i < nquad1; ++i)
            {
                d1factors[1][i] +=
                    df[2 * n + 1][(i + 1) * nquad0 - 1] * normal_1[n][i];
                d1factors[2][i] += df[2 * n + 1][i * nquad0] * normal_2[n][i];
            }
        }
    }
    else
    {
        // d xi_2/dx n_x
        for (int i = 0; i < nquad0; ++i)
        {
            d1factors[0][i] = df[1][0] * normal_0[0][i];
        }

        // d xi_1/dx n_x
        for (int i = 0; i < nquad1; ++i)
        {
            d0factors[1][i] = df[0][0] * normal_1[0][i];
            d0factors[2][i] = df[0][0] * normal_2[0][i];
        }

        for (int n = 1; n < ncoords; ++n)
        {
            // d xi_2/dy n_y
            // needs checking for 3D coords
            for (int i = 0; i < nquad0; ++i)
            {
                d1factors[0][i] += df[2 * n + 1][0] * normal_0[n][i];
            }

            // d xi_1/dy n_y
            // needs checking for 3D coords
            for (int i = 0; i < nquad1; ++i)
            {
                d0factors[1][i] += df[2 * n][0] * normal_1[n][i];
                d0factors[2][i] += df[2 * n][0] * normal_2[n][i];
            }
        }

        // d1factors
        // d xi_1/dx n_x
        for (int i = 0; i < nquad0; ++i)
        {
            d0factors[0][i] = df[0][0] * normal_0[0][i];
        }

        // d xi_2/dx n_x
        for (int i = 0; i < nquad1; ++i)
        {
            d1factors[1][i] = df[1][0] * normal_1[0][i];
            d1factors[2][i] = df[1][0] * normal_2[0][i];
        }

        for (int n = 1; n < ncoords; ++n)
        {
            // d xi_1/dy n_y
            // needs checking for 3D coords
            for (int i = 0; i < nquad0; ++i)
            {
                d0factors[0][i] += df[2 * n][0] * normal_0[n][i];
            }

            // d xi_2/dy n_y
            // needs checking for 3D coords
            for (int i = 0; i < nquad1; ++i)
            {
                d1factors[1][i] += df[2 * n + 1][0] * normal_1[n][i];
                d1factors[2][i] += df[2 * n + 1][0] * normal_2[n][i];
            }
        }
    }
}
} // namespace Nektar::LocalRegions
