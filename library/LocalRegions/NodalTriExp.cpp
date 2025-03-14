///////////////////////////////////////////////////////////////////////////////
//
// File: NodalTriExp.cpp
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
// Description: NodalTriExp routines
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Foundations/Interp.h>
#include <LocalRegions/NodalTriExp.h>

using namespace std;

namespace Nektar::LocalRegions
{
NodalTriExp::NodalTriExp(const LibUtilities::BasisKey &Ba,
                         const LibUtilities::BasisKey &Bb,
                         const LibUtilities::PointsType Ntype,
                         const SpatialDomains::TriGeomSharedPtr &geom)
    : StdExpansion(LibUtilities::StdTriData::getNumberOfCoefficients(
                       Ba.GetNumModes(), (Bb.GetNumModes())),
                   2, Ba, Bb),
      StdExpansion2D(LibUtilities::StdTriData::getNumberOfCoefficients(
                         Ba.GetNumModes(), (Bb.GetNumModes())),
                     Ba, Bb),
      StdNodalTriExp(Ba, Bb, Ntype), Expansion(geom), Expansion2D(geom),
      m_matrixManager(
          std::bind(&Expansion2D::CreateMatrix, this, std::placeholders::_1),
          std::string("NodalTriExpMatrix")),
      m_staticCondMatrixManager(std::bind(&Expansion::CreateStaticCondMatrix,
                                          this, std::placeholders::_1),
                                std::string("NodalTriExpStaticCondMatrix"))
{
}

NodalTriExp::NodalTriExp(const NodalTriExp &T)
    : StdExpansion(T), StdExpansion2D(T), StdRegions::StdTriExp(T),
      StdRegions::StdNodalTriExp(T), Expansion(T), Expansion2D(T),
      m_matrixManager(T.m_matrixManager),
      m_staticCondMatrixManager(T.m_staticCondMatrixManager)
{
}

NekDouble NodalTriExp::v_Integral(const Array<OneD, const NekDouble> &inarray)
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
    ival = StdNodalTriExp::v_Integral(tmp);
    return ival;
}

void NodalTriExp::v_PhysDeriv(const Array<OneD, const NekDouble> &inarray,
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

    StdNodalTriExp::v_PhysDeriv(inarray, diff0, diff1);

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

void NodalTriExp::v_PhysDeriv(const int dir,
                              const Array<OneD, const NekDouble> &inarray,
                              Array<OneD, NekDouble> &outarray)
{
    Array<OneD, NekDouble> tmp;
    switch (dir)
    {
        case 0:
        {
            PhysDeriv(inarray, outarray, tmp);
        }
        break;
        case 1:
        {
            PhysDeriv(inarray, tmp, outarray);
        }
        break;
        default:
        {
            ASSERTL1(dir >= 0 && dir < 2, "input dir is out of range");
        }
        break;
    }
}

void NodalTriExp::v_FwdTrans(const Array<OneD, const NekDouble> &inarray,
                             Array<OneD, NekDouble> &outarray)
{
    IProductWRTBase(inarray, outarray);

    // get Mass matrix inverse
    MatrixKey masskey(StdRegions::eInvMass, DetShapeType(), *this,
                      StdRegions::NullConstFactorMap,
                      StdRegions::NullVarCoeffMap,
                      m_nodalPointsKey.GetPointsType());
    DNekScalMatSharedPtr matsys = m_matrixManager[masskey];

    // copy inarray in case inarray == outarray
    NekVector<NekDouble> in(m_ncoeffs, outarray, eCopy);
    NekVector<NekDouble> out(m_ncoeffs, outarray, eWrapper);

    out = (*matsys) * in;
}

void NodalTriExp::v_IProductWRTBase(const Array<OneD, const NekDouble> &inarray,
                                    Array<OneD, NekDouble> &outarray)
{
    v_IProductWRTBase_SumFac(inarray, outarray);
}

void NodalTriExp::v_IProductWRTDerivBase(
    const int dir, const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray)
{
    v_IProductWRTDerivBase_SumFac(dir, inarray, outarray);
}

void NodalTriExp::v_IProductWRTBase_SumFac(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, bool multiplybyweights)
{
    int nquad0 = m_base[0]->GetNumPoints();
    int nquad1 = m_base[1]->GetNumPoints();
    int order1 = m_base[1]->GetNumModes();

    if (multiplybyweights)
    {
        Array<OneD, NekDouble> tmp(nquad0 * nquad1 + nquad0 * order1);
        Array<OneD, NekDouble> wsp(tmp + nquad0 * nquad1);

        MultiplyByQuadratureMetric(inarray, tmp);
        StdTriExp::IProductWRTBase_SumFacKernel(
            m_base[0]->GetBdata(), m_base[1]->GetBdata(), tmp, outarray, wsp);
        NodalToModalTranspose(outarray, outarray);
    }
    else
    {
        Array<OneD, NekDouble> wsp(nquad0 * order1);

        StdTriExp::IProductWRTBase_SumFacKernel(m_base[0]->GetBdata(),
                                                m_base[1]->GetBdata(), inarray,
                                                outarray, wsp);
        NodalToModalTranspose(outarray, outarray);
    }
}

void NodalTriExp::v_IProductWRTDerivBase_SumFac(
    const int dir, const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray)
{
    int nquad0  = m_base[0]->GetNumPoints();
    int nquad1  = m_base[1]->GetNumPoints();
    int nqtot   = nquad0 * nquad1;
    int wspsize = max(nqtot, m_ncoeffs);

    Array<OneD, NekDouble> tmp0(4 * wspsize);
    Array<OneD, NekDouble> tmp1(tmp0 + wspsize);
    Array<OneD, NekDouble> tmp2(tmp0 + 2 * wspsize);
    Array<OneD, NekDouble> tmp3(tmp0 + 3 * wspsize);

    Array<OneD, Array<OneD, NekDouble>> tmp2D{2};
    tmp2D[0] = tmp1;
    tmp2D[1] = tmp2;

    AlignVectorToCollapsedDir(dir, inarray, tmp2D);

    MultiplyByQuadratureMetric(tmp1, tmp1);
    MultiplyByQuadratureMetric(tmp2, tmp2);

    IProductWRTBase_SumFacKernel(m_base[0]->GetDbdata(), m_base[1]->GetBdata(),
                                 tmp1, tmp3, tmp0);
    IProductWRTBase_SumFacKernel(m_base[0]->GetBdata(), m_base[1]->GetDbdata(),
                                 tmp2, outarray, tmp0);
    Vmath::Vadd(m_ncoeffs, tmp3, 1, outarray, 1, outarray, 1);

    NodalToModalTranspose(outarray, outarray);
}

void NodalTriExp::v_AlignVectorToCollapsedDir(
    const int dir, const Array<OneD, const NekDouble> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    ASSERTL1((dir == 0) || (dir == 1) || (dir == 2), "Invalid direction.");
    ASSERTL1((dir == 2) ? (m_geom->GetCoordim() == 3) : true,
             "Invalid direction.");

    int nquad0  = m_base[0]->GetNumPoints();
    int nquad1  = m_base[1]->GetNumPoints();
    int nqtot   = nquad0 * nquad1;
    int wspsize = max(nqtot, m_ncoeffs);

    const Array<TwoD, const NekDouble> &df =
        m_metricinfo->GetDerivFactors(GetPointsKeys());

    Array<OneD, NekDouble> tmp0(4 * wspsize);
    Array<OneD, NekDouble> tmp3(tmp0 + wspsize);
    Array<OneD, NekDouble> gfac0(tmp0 + 2 * wspsize);
    Array<OneD, NekDouble> gfac1(tmp0 + 3 * wspsize);

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

StdRegions::StdExpansionSharedPtr NodalTriExp::v_GetStdExp(void) const
{

    return MemoryManager<StdRegions::StdNodalTriExp>::AllocateSharedPtr(
        m_base[0]->GetBasisKey(), m_base[1]->GetBasisKey(),
        m_nodalPointsKey.GetPointsType());
}

StdRegions::StdExpansionSharedPtr NodalTriExp::v_GetLinStdExp(void) const
{
    LibUtilities::BasisKey bkey0(m_base[0]->GetBasisType(), 2,
                                 m_base[0]->GetPointsKey());
    LibUtilities::BasisKey bkey1(m_base[1]->GetBasisType(), 2,
                                 m_base[1]->GetPointsKey());

    return MemoryManager<StdRegions::StdNodalTriExp>::AllocateSharedPtr(
        bkey0, bkey1, m_nodalPointsKey.GetPointsType());
}

void NodalTriExp::v_GetCoords(Array<OneD, NekDouble> &coords_0,
                              Array<OneD, NekDouble> &coords_1,
                              Array<OneD, NekDouble> &coords_2)
{
    Expansion::v_GetCoords(coords_0, coords_1, coords_2);
}

void NodalTriExp::v_GetCoord(const Array<OneD, const NekDouble> &Lcoords,
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

NekDouble NodalTriExp::v_PhysEvaluate(
    const Array<OneD, const NekDouble> &coord,
    const Array<OneD, const NekDouble> &physvals)

{
    Array<OneD, NekDouble> Lcoord = Array<OneD, NekDouble>(2);

    ASSERTL0(m_geom, "m_geom not defined");
    m_geom->GetLocCoords(coord, Lcoord);

    return StdExpansion2D::v_PhysEvaluate(Lcoord, physvals);
}

void NodalTriExp::v_ComputeTraceNormal(const int edge)
{
    int i;
    const SpatialDomains::GeomFactorsSharedPtr &geomFactors =
        GetGeom()->GetMetricInfo();
    const SpatialDomains::GeomType type = geomFactors->GetGtype();

    LibUtilities::PointsKeyVector ptsKeys = GetPointsKeys();
    const Array<TwoD, const NekDouble> &df =
        geomFactors->GetDerivFactors(ptsKeys);
    const Array<OneD, const NekDouble> &jac = geomFactors->GetJac(ptsKeys);
    int nqe                                 = m_base[0]->GetNumPoints();
    int dim                                 = GetCoordim();

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
        LibUtilities::Interp1D(from_key, jac, m_base[0]->GetPointsKey(), work);
        Vmath::Sdiv(nq, 1.0, &work[0], 1, &work[0], 1);

        // interpolate
        for (i = 0; i < GetCoordim(); ++i)
        {
            LibUtilities::Interp1D(from_key, &normals[i * nq],
                                   m_base[0]->GetPointsKey(), &normal[i][0]);
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

        // Reverse direction so that points are in
        // anticlockwise direction if edge >=2
        if (edge >= 2)
        {
            for (i = 0; i < GetCoordim(); ++i)
            {
                Vmath::Reverse(nqe, normal[i], 1, normal[i], 1);
            }
        }
    }
}

void NodalTriExp::v_ExtractDataToCoeffs(
    const NekDouble *data, const std::vector<unsigned int> &nummodes,
    const int mode_offset, NekDouble *coeffs,
    [[maybe_unused]] std::vector<LibUtilities::BasisType> &fromType)
{
    Array<OneD, NekDouble> modes(m_ncoeffs);
    Expansion::ExtractDataToCoeffs(data, nummodes, mode_offset, &modes[0],
                                   fromType);

    Array<OneD, NekDouble> nodes(m_ncoeffs, coeffs, eArrayWrapper);
    ModalToNodal(modes, nodes);
}

void NodalTriExp::v_GetTracePhysVals(
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

DNekMatSharedPtr NodalTriExp::v_GenMatrix(const StdRegions::StdMatrixKey &mkey)
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
            returnval = Expansion2D::v_GenMatrix(mkey);
            break;
        default:
            returnval = StdNodalTriExp::v_GenMatrix(mkey);
            break;
    }
    return returnval;
}

DNekMatSharedPtr NodalTriExp::v_CreateStdMatrix(
    const StdRegions::StdMatrixKey &mkey)
{
    LibUtilities::BasisKey bkey0   = m_base[0]->GetBasisKey();
    LibUtilities::BasisKey bkey1   = m_base[1]->GetBasisKey();
    LibUtilities::PointsType ntype = m_nodalPointsKey.GetPointsType();
    StdRegions::StdNodalTriExpSharedPtr tmp =
        MemoryManager<StdNodalTriExp>::AllocateSharedPtr(bkey0, bkey1, ntype);

    return tmp->GetStdMatrix(mkey);
}

DNekScalMatSharedPtr NodalTriExp::v_GetLocMatrix(const MatrixKey &mkey)
{
    return m_matrixManager[mkey];
}
DNekScalBlkMatSharedPtr NodalTriExp::v_GetLocStaticCondMatrix(
    const MatrixKey &mkey)
{
    return m_staticCondMatrixManager[mkey];
}
void NodalTriExp::v_DropLocMatrix(const MatrixKey &mkey)
{
    m_matrixManager.DeleteObject(mkey);
}

void NodalTriExp::v_MassMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                 Array<OneD, NekDouble> &outarray,
                                 const StdRegions::StdMatrixKey &mkey)
{
    StdExpansion::MassMatrixOp_MatFree(inarray, outarray, mkey);
}

void NodalTriExp::v_LaplacianMatrixOp(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const StdRegions::StdMatrixKey &mkey)
{
    StdExpansion::LaplacianMatrixOp_MatFree_GenericImpl(inarray, outarray,
                                                        mkey);
}

void NodalTriExp::v_LaplacianMatrixOp(
    const int k1, const int k2, const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const StdRegions::StdMatrixKey &mkey)
{
    StdExpansion::LaplacianMatrixOp_MatFree(k1, k2, inarray, outarray, mkey);
}

void NodalTriExp::v_WeakDerivMatrixOp(
    const int i, const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const StdRegions::StdMatrixKey &mkey)
{
    StdExpansion::WeakDerivMatrixOp_MatFree(i, inarray, outarray, mkey);
}

void NodalTriExp::v_HelmholtzMatrixOp(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const StdRegions::StdMatrixKey &mkey)
{
    StdExpansion::HelmholtzMatrixOp_MatFree_GenericImpl(inarray, outarray,
                                                        mkey);
}

} // namespace Nektar::LocalRegions
