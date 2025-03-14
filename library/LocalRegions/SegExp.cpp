///////////////////////////////////////////////////////////////////////////////
//
// File: SegExp.cpp
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
// Description: SegExp routines
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Foundations/Interp.h>
#include <LocalRegions/Expansion2D.h>
#include <LocalRegions/SegExp.h>

using namespace std;

namespace Nektar::LocalRegions
{

/**
 * @class SegExp
 * Defines a Segment local expansion.
 */

/// Constructor using BasisKey class for quadrature points and
/// order definition.
/**
 * @param   Ba          Basis key of segment expansion.
 * @param   geom        Description of geometry.
 */
SegExp::SegExp(const LibUtilities::BasisKey &Ba,
               const SpatialDomains::Geometry1DSharedPtr &geom)
    : StdExpansion(Ba.GetNumModes(), 1, Ba),
      StdExpansion1D(Ba.GetNumModes(), Ba), StdRegions::StdSegExp(Ba),
      Expansion(geom), Expansion1D(geom),
      m_matrixManager(
          std::bind(&SegExp::CreateMatrix, this, std::placeholders::_1),
          std::string("SegExpMatrix")),
      m_staticCondMatrixManager(std::bind(&Expansion::CreateStaticCondMatrix,
                                          this, std::placeholders::_1),
                                std::string("SegExpStaticCondMatrix"))
{
}

/// Copy Constructor
/**
 * @param   S           Existing segment to duplicate.
 */
SegExp::SegExp(const SegExp &S)
    : StdExpansion(S), StdExpansion1D(S), StdRegions::StdSegExp(S),
      Expansion(S), Expansion1D(S), m_matrixManager(S.m_matrixManager),
      m_staticCondMatrixManager(S.m_staticCondMatrixManager)
{
}

//----------------------------
// Integration Methods
//----------------------------

/** \brief Integrate the physical point list \a inarray over region
    and return the value

    Inputs:\n

    - \a inarray: definition of function to be returned at
    quadrature point of expansion.

    Outputs:\n

    - returns \f$\int^1_{-1} u(\xi_1)d \xi_1 \f$ where \f$inarray[i]
    = u(\xi_{1i}) \f$
*/

NekDouble SegExp::v_Integral(const Array<OneD, const NekDouble> &inarray)
{
    int nquad0                       = m_base[0]->GetNumPoints();
    Array<OneD, const NekDouble> jac = m_metricinfo->GetJac(GetPointsKeys());
    NekDouble ival;
    Array<OneD, NekDouble> tmp(nquad0);

    // multiply inarray with Jacobian
    if (m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
    {
        Vmath::Vmul(nquad0, jac, 1, inarray, 1, tmp, 1);
    }
    else
    {
        Vmath::Smul(nquad0, jac[0], inarray, 1, tmp, 1);
    }

    // call StdSegExp version;
    ival = StdSegExp::v_Integral(tmp);
    // ival = StdSegExp::Integral(tmp);
    return ival;
}

//-----------------------------
// Differentiation Methods
//-----------------------------

/** \brief Evaluate the derivative \f$ d/d{\xi_1} \f$ at the
    physical quadrature points given by \a inarray and return in \a
    outarray.

    This is a wrapper around StdExpansion1D::Tensor_Deriv

    Input:\n

    - \a n: number of derivatives to be evaluated where \f$ n \leq  dim\f$

    - \a inarray: array of function evaluated at the quadrature points

    Output: \n

    - \a outarray: array of the derivatives \f$
    du/d_{\xi_1}|_{\xi_{1i}} d\xi_1/dx,
    du/d_{\xi_1}|_{\xi_{1i}} d\xi_1/dy,
    du/d_{\xi_1}|_{\xi_{1i}} d\xi_1/dz,
    \f$ depending on value of \a dim
*/
void SegExp::v_PhysDeriv(const Array<OneD, const NekDouble> &inarray,
                         Array<OneD, NekDouble> &out_d0,
                         Array<OneD, NekDouble> &out_d1,
                         Array<OneD, NekDouble> &out_d2)
{
    int nquad0 = m_base[0]->GetNumPoints();
    Array<TwoD, const NekDouble> gmat =
        m_metricinfo->GetDerivFactors(GetPointsKeys());
    Array<OneD, NekDouble> diff(nquad0);

    // StdExpansion1D::PhysTensorDeriv(inarray,diff);
    PhysTensorDeriv(inarray, diff);
    if (m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
    {
        if (out_d0.size())
        {
            Vmath::Vmul(nquad0, &gmat[0][0], 1, &diff[0], 1, &out_d0[0], 1);
        }

        if (out_d1.size())
        {
            Vmath::Vmul(nquad0, &gmat[1][0], 1, &diff[0], 1, &out_d1[0], 1);
        }

        if (out_d2.size())
        {
            Vmath::Vmul(nquad0, &gmat[2][0], 1, &diff[0], 1, &out_d2[0], 1);
        }
    }
    else
    {
        if (out_d0.size())
        {
            Vmath::Smul(nquad0, gmat[0][0], diff, 1, out_d0, 1);
        }

        if (out_d1.size())
        {
            Vmath::Smul(nquad0, gmat[1][0], diff, 1, out_d1, 1);
        }

        if (out_d2.size())
        {
            Vmath::Smul(nquad0, gmat[2][0], diff, 1, out_d2, 1);
        }
    }
}

/**
 *\brief Evaluate the derivative along a line:
 * \f$ d/ds=\frac{spacedim}{||tangent||}d/d{\xi}  \f$.
 * The derivative is calculated performing
 *the product \f$ du/d{s}=\nabla u \cdot tangent \f$.
 *\param inarray function to derive
 *\param out_ds result of the derivative operation
 **/
void SegExp::v_PhysDeriv_s(const Array<OneD, const NekDouble> &inarray,
                           Array<OneD, NekDouble> &out_ds)
{
    int nquad0  = m_base[0]->GetNumPoints();
    int coordim = m_geom->GetCoordim();
    Array<OneD, NekDouble> diff(nquad0);
    // this operation is needed if you put out_ds==inarray
    Vmath::Zero(nquad0, out_ds, 1);
    switch (coordim)
    {
        case 2:
            // diff= dU/de
            Array<OneD, NekDouble> diff(nquad0);

            PhysTensorDeriv(inarray, diff);

            // get dS/de= (Jac)^-1
            Array<OneD, NekDouble> Jac = m_metricinfo->GetJac(GetPointsKeys());
            if (m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                // calculate the derivative as (dU/de)*(Jac)^-1
                Vmath::Vdiv(nquad0, diff, 1, Jac, 1, out_ds, 1);
            }
            else
            {
                NekDouble invJac = 1 / Jac[0];
                Vmath::Smul(nquad0, invJac, diff, 1, out_ds, 1);
            }
    }
}

/**
 *\brief Evaluate the derivative normal to a line:
 * \f$ d/dn=\frac{spacedim}{||normal||}d/d{\xi}  \f$.
 * The derivative is calculated performing
 *the product \f$ du/d{s}=\nabla u \cdot normal \f$.
 *\param inarray function to derive
 *\param out_dn result of the derivative operation
 **/
void SegExp::v_PhysDeriv_n(const Array<OneD, const NekDouble> &inarray,
                           Array<OneD, NekDouble> &out_dn)
{
    int nquad0 = m_base[0]->GetNumPoints();
    Array<TwoD, const NekDouble> gmat =
        m_metricinfo->GetDerivFactors(GetPointsKeys());
    int coordim = m_geom->GetCoordim();
    Array<OneD, NekDouble> out_dn_tmp(nquad0, 0.0);
    switch (coordim)
    {
        case 2:

            Array<OneD, NekDouble> inarray_d0(nquad0);
            Array<OneD, NekDouble> inarray_d1(nquad0);

            v_PhysDeriv(inarray, inarray_d0, inarray_d1);
            Array<OneD, Array<OneD, NekDouble>> normals;
            normals = Array<OneD, Array<OneD, NekDouble>>(coordim);
            cout << "der_n" << endl;
            for (int k = 0; k < coordim; ++k)
            {
                normals[k] = Array<OneD, NekDouble>(nquad0);
            }
            // @TODO: this routine no longer makes sense, since normals are not
            // unique on
            //        an edge
            //        normals = GetMetricInfo()->GetNormal();
            for (int i = 0; i < nquad0; i++)
            {
                cout << "nx= " << normals[0][i] << "  ny=" << normals[1][i]
                     << endl;
            }
            ASSERTL0(normals != NullNekDoubleArrayOfArray,
                     "normal vectors do not exist: check if a"
                     "boundary region is defined as I ");
            // \nabla u \cdot normal
            Vmath::Vmul(nquad0, normals[0], 1, inarray_d0, 1, out_dn_tmp, 1);
            Vmath::Vadd(nquad0, out_dn_tmp, 1, out_dn, 1, out_dn, 1);
            Vmath::Zero(nquad0, out_dn_tmp, 1);
            Vmath::Vmul(nquad0, normals[1], 1, inarray_d1, 1, out_dn_tmp, 1);
            Vmath::Vadd(nquad0, out_dn_tmp, 1, out_dn, 1, out_dn, 1);

            for (int i = 0; i < nquad0; i++)
            {
                cout << "deps/dx =" << inarray_d0[i]
                     << "  deps/dy=" << inarray_d1[i] << endl;
            }
    }
}
void SegExp::v_PhysDeriv(const int dir,
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

//-----------------------------
// Transforms
//-----------------------------

/** \brief Forward transform from physical quadrature space
    stored in \a inarray and evaluate the expansion coefficients and
    store in \a outarray

    Perform a forward transform using a Galerkin projection by
    taking the inner product of the physical points and multiplying
    by the inverse of the mass matrix using the Solve method of the
    standard matrix container holding the local mass matrix, i.e.
    \f$ {\bf \hat{u}} = {\bf M}^{-1} {\bf I} \f$ where \f$ {\bf I}[p] =
    \int^1_{-1} \phi_p(\xi_1) u(\xi_1) d\xi_1 \f$

    Inputs:\n

    - \a inarray: array of physical quadrature points to be transformed

    Outputs:\n

    - \a outarray: updated array of expansion coefficients.

*/
// need to sort out family of matrices
void SegExp::v_FwdTrans(const Array<OneD, const NekDouble> &inarray,
                        Array<OneD, NekDouble> &outarray)
{
    if (m_base[0]->Collocation())
    {
        Vmath::Vcopy(m_ncoeffs, inarray, 1, outarray, 1);
    }
    else
    {
        v_IProductWRTBase(inarray, outarray);

        // get Mass matrix inverse
        MatrixKey masskey(StdRegions::eInvMass, DetShapeType(), *this);
        DNekScalMatSharedPtr matsys = m_matrixManager[masskey];

        // copy inarray in case inarray == outarray
        NekVector<NekDouble> in(m_ncoeffs, outarray, eCopy);
        NekVector<NekDouble> out(m_ncoeffs, outarray, eWrapper);

        out = (*matsys) * in;
    }
}

void SegExp::v_FwdTransBndConstrained(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray)
{
    if (m_base[0]->Collocation())
    {
        Vmath::Vcopy(m_ncoeffs, inarray, 1, outarray, 1);
    }
    else
    {
        int nInteriorDofs = m_ncoeffs - 2;
        int offset        = 0;

        switch (m_base[0]->GetBasisType())
        {
            case LibUtilities::eGLL_Lagrange:
            {
                offset = 1;
            }
            break;
            case LibUtilities::eGauss_Lagrange:
            {
                nInteriorDofs = m_ncoeffs;
                offset        = 0;
            }
            break;
            case LibUtilities::eModified_A:
            case LibUtilities::eModified_B:
            {
                ASSERTL1(
                    m_base[0]->GetPointsType() ==
                            LibUtilities::eGaussLobattoLegendre ||
                        m_base[0]->GetPointsType() ==
                            LibUtilities::ePolyEvenlySpaced,
                    "Cannot use FwdTrans_BndConstrained with these points.");
                offset = 2;
            }
            break;
            default:
                ASSERTL0(false, "This type of FwdTrans is not defined"
                                "for this expansion type");
        }

        fill(outarray.data(), outarray.data() + m_ncoeffs, 0.0);

        if (m_base[0]->GetBasisType() != LibUtilities::eGauss_Lagrange)
        {

            outarray[GetVertexMap(0)] = inarray[0];
            outarray[GetVertexMap(1)] = inarray[m_base[0]->GetNumPoints() - 1];

            if (m_ncoeffs > 2)
            {
                // ideally, we would like to have tmp0 to be replaced
                // by outarray (currently MassMatrixOp does not allow
                // aliasing)
                Array<OneD, NekDouble> tmp0(m_ncoeffs);
                Array<OneD, NekDouble> tmp1(m_ncoeffs);

                StdRegions::StdMatrixKey stdmasskey(StdRegions::eMass,
                                                    DetShapeType(), *this);
                MassMatrixOp(outarray, tmp0, stdmasskey);
                v_IProductWRTBase(inarray, tmp1);

                Vmath::Vsub(m_ncoeffs, tmp1, 1, tmp0, 1, tmp1, 1);

                // get Mass matrix inverse (only of interior DOF)
                MatrixKey masskey(StdRegions::eMass, DetShapeType(), *this);
                DNekScalMatSharedPtr matsys =
                    (m_staticCondMatrixManager[masskey])->GetBlock(1, 1);

                Blas::Dgemv('N', nInteriorDofs, nInteriorDofs, matsys->Scale(),
                            &((matsys->GetOwnedMatrix())->GetPtr())[0],
                            nInteriorDofs, tmp1.data() + offset, 1, 0.0,
                            outarray.data() + offset, 1);
            }
        }
        else
        {
            SegExp::v_FwdTrans(inarray, outarray);
        }
    }
}

//-----------------------------
// Inner product functions
//-----------------------------

/** \brief  Inner product of \a inarray over region with respect to
    the expansion basis (this)->_Base[0] and return in \a outarray

    Wrapper call to SegExp::IProduct_WRT_B

    Input:\n

    - \a inarray: array of function evaluated at the physical
    collocation points

    Output:\n

    - \a outarray: array of inner product with respect to each
    basis over region
*/
void SegExp::v_IProductWRTBase(const Array<OneD, const NekDouble> &inarray,
                               Array<OneD, NekDouble> &outarray)
{
    v_IProductWRTBase(m_base[0]->GetBdata(), inarray, outarray, 1);
}

/**
   \brief  Inner product of \a inarray over region with respect to
   expansion basis \a base and return in \a outarray

   Calculate \f$ I[p] = \int^{1}_{-1} \phi_p(\xi_1) u(\xi_1) d\xi_1
   = \sum_{i=0}^{nq-1} \phi_p(\xi_{1i}) u(\xi_{1i}) w_i \f$ where
   \f$ outarray[p] = I[p], inarray[i] = u(\xi_{1i}), base[p*nq+i] =
   \phi_p(\xi_{1i}) \f$.

   Inputs: \n

   - \a base: an array definiing the local basis for the inner
   product usually passed from Basis->get_bdata() or
   Basis->get_Dbdata()
   - \a inarray: physical point array of function to be integrated
   \f$ u(\xi_1) \f$
   - \a coll_check: Flag to identify when a Basis->collocation()
   call should be performed to see if this is a GLL_Lagrange basis
   with a collocation property. (should be set to 0 if taking the
   inner product with respect to the derivative of basis)

   Output: \n

   - \a outarray: array of coefficients representing the inner
   product of function with ever  mode in the exapnsion

**/
void SegExp::v_IProductWRTBase(const Array<OneD, const NekDouble> &base,
                               const Array<OneD, const NekDouble> &inarray,
                               Array<OneD, NekDouble> &outarray, int coll_check)
{
    int nquad0                       = m_base[0]->GetNumPoints();
    Array<OneD, const NekDouble> jac = m_metricinfo->GetJac(GetPointsKeys());
    Array<OneD, NekDouble> tmp(nquad0);

    // multiply inarray with Jacobian
    if (m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
    {
        Vmath::Vmul(nquad0, jac, 1, inarray, 1, tmp, 1);
    }
    else
    {
        Vmath::Smul(nquad0, jac[0], inarray, 1, tmp, 1);
    }
    StdSegExp::v_IProductWRTBase(base, tmp, outarray, coll_check);
}

void SegExp::v_IProductWRTDerivBase(const int dir,
                                    const Array<OneD, const NekDouble> &inarray,
                                    Array<OneD, NekDouble> &outarray)
{
    ASSERTL1(dir < 3, "input dir is out of range");
    ASSERTL1((dir == 2) ? m_geom->GetCoordim() == 3 : true,
             "input dir is out of range");

    int nquad = m_base[0]->GetNumPoints();
    const Array<TwoD, const NekDouble> &gmat =
        m_metricinfo->GetDerivFactors(GetPointsKeys());

    Array<OneD, NekDouble> tmp1(nquad);

    if (m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
    {
        Vmath::Vmul(nquad, gmat[dir], 1, inarray, 1, tmp1, 1);
    }
    else
    {
        Vmath::Smul(nquad, gmat[dir][0], inarray, 1, tmp1, 1);
    }

    v_IProductWRTBase(m_base[0]->GetDbdata(), tmp1, outarray, 1);
}

void SegExp::v_NormVectorIProductWRTBase(const Array<OneD, const NekDouble> &Fx,
                                         const Array<OneD, const NekDouble> &Fy,
                                         Array<OneD, NekDouble> &outarray)
{
    int nq = m_base[0]->GetNumPoints();
    Array<OneD, NekDouble> Fn(nq);

    // @TODO: This routine no longer makes sense as a normal is not unique to an
    // edge
    const Array<OneD, const Array<OneD, NekDouble>> &normals =
        GetLeftAdjacentElementExp()->GetTraceNormal(
            GetLeftAdjacentElementTrace());
    Vmath::Vmul(nq, &Fx[0], 1, &normals[0][0], 1, &Fn[0], 1);
    Vmath::Vvtvp(nq, &Fy[0], 1, &normals[1][0], 1, &Fn[0], 1, &Fn[0], 1);

    v_IProductWRTBase(Fn, outarray);
}

void SegExp::v_NormVectorIProductWRTBase(
    const Array<OneD, const Array<OneD, NekDouble>> &Fvec,
    Array<OneD, NekDouble> &outarray)
{
    NormVectorIProductWRTBase(Fvec[0], Fvec[1], outarray);
}

//-----------------------------
// Evaluation functions
//-----------------------------

/**
 * Given the local cartesian coordinate \a Lcoord evaluate the
 * value of physvals at this point by calling through to the
 * StdExpansion method
 */
NekDouble SegExp::v_StdPhysEvaluate(
    const Array<OneD, const NekDouble> &Lcoord,
    const Array<OneD, const NekDouble> &physvals)
{
    // Evaluate point in local (eta) coordinates.
    return StdExpansion1D::v_PhysEvaluate(Lcoord, physvals);
}

NekDouble SegExp::v_PhysEvaluate(const Array<OneD, const NekDouble> &coord,
                                 const Array<OneD, const NekDouble> &physvals)
{
    Array<OneD, NekDouble> Lcoord = Array<OneD, NekDouble>(1);

    ASSERTL0(m_geom, "m_geom not defined");
    m_geom->GetLocCoords(coord, Lcoord);

    return StdExpansion1D::v_PhysEvaluate(Lcoord, physvals);
}

NekDouble SegExp::v_PhysEvalFirstDeriv(
    const Array<OneD, NekDouble> &coord,
    const Array<OneD, const NekDouble> &inarray,
    std::array<NekDouble, 3> &firstOrderDerivs)
{
    Array<OneD, NekDouble> Lcoord(1);
    ASSERTL0(m_geom, "m_geom not defined");
    m_geom->GetLocCoords(coord, Lcoord);
    return StdSegExp::v_PhysEvalFirstDeriv(Lcoord, inarray, firstOrderDerivs);
}

NekDouble SegExp::v_PhysEvalFirstSecondDeriv(
    const Array<OneD, NekDouble> &coord,
    const Array<OneD, const NekDouble> &inarray,
    std::array<NekDouble, 3> &firstOrderDerivs,
    std::array<NekDouble, 6> &secondOrderDerivs)
{
    Array<OneD, NekDouble> Lcoord(1);
    ASSERTL0(m_geom, "m_geom not defined");
    m_geom->GetLocCoords(coord, Lcoord);
    return StdSegExp::v_PhysEvalFirstSecondDeriv(
        Lcoord, inarray, firstOrderDerivs, secondOrderDerivs);
}

void SegExp::v_GetCoord(const Array<OneD, const NekDouble> &Lcoords,
                        Array<OneD, NekDouble> &coords)
{
    int i;

    ASSERTL1(Lcoords[0] >= -1.0 && Lcoords[0] <= 1.0,
             "Local coordinates are not in region [-1,1]");

    m_geom->FillGeom();
    for (i = 0; i < m_geom->GetCoordim(); ++i)
    {
        coords[i] = m_geom->GetCoord(i, Lcoords);
    }
}

void SegExp::v_GetCoords(Array<OneD, NekDouble> &coords_0,
                         Array<OneD, NekDouble> &coords_1,
                         Array<OneD, NekDouble> &coords_2)
{
    Expansion::v_GetCoords(coords_0, coords_1, coords_2);
}

// Get vertex value from the 1D Phys space.
void SegExp::v_GetVertexPhysVals(const int vertex,
                                 const Array<OneD, const NekDouble> &inarray,
                                 NekDouble &outarray)
{
    int nquad = m_base[0]->GetNumPoints();

    if (m_base[0]->GetPointsType() != LibUtilities::eGaussGaussLegendre)
    {
        switch (vertex)
        {
            case 0:
                outarray = inarray[0];
                break;
            case 1:
                outarray = inarray[nquad - 1];
                break;
        }
    }
    else
    {
        StdRegions::ConstFactorMap factors;
        factors[StdRegions::eFactorGaussVertex] = vertex;

        StdRegions::StdMatrixKey key(StdRegions::eInterpGauss, DetShapeType(),
                                     *this, factors);

        DNekScalMatSharedPtr mat_gauss = m_matrixManager[key];

        outarray =
            Blas::Ddot(nquad, mat_gauss->GetOwnedMatrix()->GetPtr().data(), 1,
                       &inarray[0], 1);
    }
}

// Get vertex value from the 1D Phys space.
void SegExp::v_GetTracePhysVals(
    const int edge,
    [[maybe_unused]] const StdRegions::StdExpansionSharedPtr &EdgeExp,
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray,
    [[maybe_unused]] StdRegions::Orientation orient)
{
    NekDouble result;
    v_GetVertexPhysVals(edge, inarray, result);
    outarray[0] = result;
}

// Get vertex map from the 1D Phys space.
void SegExp::v_GetTracePhysMap(const int vertex, Array<OneD, int> &map)
{
    int nquad = m_base[0]->GetNumPoints();

    ASSERTL1(vertex == 0 || vertex == 1, "Vertex value should be 0 or 1");

    map = Array<OneD, int>(1);

    map[0] = vertex == 0 ? 0 : nquad - 1;
}

//-----------------------------
// Helper functions
//-----------------------------

void SegExp::v_SetCoeffsToOrientation(StdRegions::Orientation dir,
                                      Array<OneD, const NekDouble> &inarray,
                                      Array<OneD, NekDouble> &outarray)
{

    if (dir == StdRegions::eBackwards)
    {
        if (&inarray[0] != &outarray[0])
        {
            Array<OneD, NekDouble> intmp(inarray);
            ReverseCoeffsAndSign(intmp, outarray);
        }
        else
        {
            ReverseCoeffsAndSign(inarray, outarray);
        }
    }
}

StdRegions::StdExpansionSharedPtr SegExp::v_GetStdExp() const
{
    return MemoryManager<StdRegions::StdSegExp>::AllocateSharedPtr(
        m_base[0]->GetBasisKey());
}

StdRegions::StdExpansionSharedPtr SegExp::v_GetLinStdExp(void) const
{
    LibUtilities::BasisKey bkey0(m_base[0]->GetBasisType(), 2,
                                 m_base[0]->GetPointsKey());

    return MemoryManager<StdRegions::StdSegExp>::AllocateSharedPtr(bkey0);
}

const Array<OneD, const NekDouble> &SegExp::v_GetPhysNormals(void)
{
    NEKERROR(ErrorUtil::efatal, "Got to SegExp");
    return NullNekDouble1DArray;
}

int SegExp::v_NumBndryCoeffs() const
{
    return 2;
}

int SegExp::v_NumDGBndryCoeffs() const
{
    return 2;
}

/// Unpack data from input file assuming it comes from
// the same expansion type
void SegExp::v_ExtractDataToCoeffs(
    const NekDouble *data, const std::vector<unsigned int> &nummodes,
    const int mode_offset, NekDouble *coeffs,
    [[maybe_unused]] std::vector<LibUtilities::BasisType> &fromType)
{
    switch (m_base[0]->GetBasisType())
    {
        case LibUtilities::eModified_A:
        {
            int fillorder = min((int)nummodes[mode_offset], m_ncoeffs);

            Vmath::Zero(m_ncoeffs, coeffs, 1);
            Vmath::Vcopy(fillorder, &data[0], 1, &coeffs[0], 1);
        }
        break;
        case LibUtilities::eGLL_Lagrange:
        {
            // Assume that input is also Gll_Lagrange
            // but no way to check;
            LibUtilities::PointsKey f0(nummodes[mode_offset],
                                       LibUtilities::eGaussLobattoLegendre);
            LibUtilities::PointsKey t0(m_base[0]->GetNumModes(),
                                       LibUtilities::eGaussLobattoLegendre);
            LibUtilities::Interp1D(f0, data, t0, coeffs);
        }
        break;
        case LibUtilities::eGauss_Lagrange:
        {
            // Assume that input is also Gauss_Lagrange
            // but no way to check;
            LibUtilities::PointsKey f0(nummodes[mode_offset],
                                       LibUtilities::eGaussGaussLegendre);
            LibUtilities::PointsKey t0(m_base[0]->GetNumModes(),
                                       LibUtilities::eGaussGaussLegendre);
            LibUtilities::Interp1D(f0, data, t0, coeffs);
        }
        break;
        default:
            ASSERTL0(false, "basis is either not set up or not hierarchicial");
    }
}

void SegExp::v_ComputeTraceNormal(const int vertex)
{
    int i;
    const SpatialDomains::GeomFactorsSharedPtr &geomFactors =
        GetGeom()->GetMetricInfo();
    SpatialDomains::GeomType type = geomFactors->GetGtype();
    const Array<TwoD, const NekDouble> &gmat =
        geomFactors->GetDerivFactors(GetPointsKeys());
    int nqe       = 1;
    int vCoordDim = GetCoordim();

    m_traceNormals[vertex] = Array<OneD, Array<OneD, NekDouble>>(vCoordDim);
    Array<OneD, Array<OneD, NekDouble>> &normal = m_traceNormals[vertex];
    for (i = 0; i < vCoordDim; ++i)
    {
        normal[i] = Array<OneD, NekDouble>(nqe);
    }

    size_t nqb                     = nqe;
    size_t nbnd                    = vertex;
    m_elmtBndNormDirElmtLen[nbnd]  = Array<OneD, NekDouble>{nqb, 0.0};
    Array<OneD, NekDouble> &length = m_elmtBndNormDirElmtLen[nbnd];

    // Regular geometry case
    if ((type == SpatialDomains::eRegular) ||
        (type == SpatialDomains::eMovingRegular))
    {
        NekDouble vert;
        // Set up normals
        switch (vertex)
        {
            case 0:
                for (i = 0; i < vCoordDim; ++i)
                {
                    Vmath::Fill(nqe, -gmat[i][0], normal[i], 1);
                }
                break;
            case 1:
                for (i = 0; i < vCoordDim; ++i)
                {
                    Vmath::Fill(nqe, gmat[i][0], normal[i], 1);
                }
                break;
            default:
                ASSERTL0(false, "point is out of range (point < 2)");
        }

        // normalise
        vert = 0.0;
        for (i = 0; i < vCoordDim; ++i)
        {
            vert += normal[i][0] * normal[i][0];
        }
        vert = 1.0 / sqrt(vert);

        Vmath::Fill(nqb, vert, length, 1);

        for (i = 0; i < vCoordDim; ++i)
        {
            Vmath::Smul(nqe, vert, normal[i], 1, normal[i], 1);
        }
    }
}

//-----------------------------
// Operator creation functions
//-----------------------------

void SegExp::v_LaplacianMatrixOp(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray,
    [[maybe_unused]] const StdRegions::StdMatrixKey &mkey)
{
    int nquad = m_base[0]->GetNumPoints();
    const Array<TwoD, const NekDouble> &gmat =
        m_metricinfo->GetDerivFactors(GetPointsKeys());

    Array<OneD, NekDouble> physValues(nquad);
    Array<OneD, NekDouble> dPhysValuesdx(nquad);

    BwdTrans(inarray, physValues);

    // Laplacian matrix operation
    switch (m_geom->GetCoordim())
    {
        case 1:
        {
            PhysDeriv(physValues, dPhysValuesdx);

            // multiply with the proper geometric factors
            if (m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                Vmath::Vmul(nquad, &gmat[0][0], 1, dPhysValuesdx.data(), 1,
                            dPhysValuesdx.data(), 1);
            }
            else
            {
                Blas::Dscal(nquad, gmat[0][0], dPhysValuesdx.data(), 1);
            }
        }
        break;
        case 2:
        {
            Array<OneD, NekDouble> dPhysValuesdy(nquad);

            PhysDeriv(physValues, dPhysValuesdx, dPhysValuesdy);

            // multiply with the proper geometric factors
            if (m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                Vmath::Vmul(nquad, &gmat[0][0], 1, dPhysValuesdx.data(), 1,
                            dPhysValuesdx.data(), 1);
                Vmath::Vvtvp(nquad, &gmat[1][0], 1, dPhysValuesdy.data(), 1,
                             dPhysValuesdx.data(), 1, dPhysValuesdx.data(), 1);
            }
            else
            {
                Blas::Dscal(nquad, gmat[0][0], dPhysValuesdx.data(), 1);
                Blas::Daxpy(nquad, gmat[1][0], dPhysValuesdy.data(), 1,
                            dPhysValuesdx.data(), 1);
            }
        }
        break;
        case 3:
        {
            Array<OneD, NekDouble> dPhysValuesdy(nquad);
            Array<OneD, NekDouble> dPhysValuesdz(nquad);

            PhysDeriv(physValues, dPhysValuesdx, dPhysValuesdy, dPhysValuesdz);

            // multiply with the proper geometric factors
            if (m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                Vmath::Vmul(nquad, &gmat[0][0], 1, dPhysValuesdx.data(), 1,
                            dPhysValuesdx.data(), 1);
                Vmath::Vvtvp(nquad, &gmat[1][0], 1, dPhysValuesdy.data(), 1,
                             dPhysValuesdx.data(), 1, dPhysValuesdx.data(), 1);
                Vmath::Vvtvp(nquad, &gmat[2][0], 1, dPhysValuesdz.data(), 1,
                             dPhysValuesdx.data(), 1, dPhysValuesdx.data(), 1);
            }
            else
            {
                Blas::Dscal(nquad, gmat[0][0], dPhysValuesdx.data(), 1);
                Blas::Daxpy(nquad, gmat[1][0], dPhysValuesdy.data(), 1,
                            dPhysValuesdx.data(), 1);
                Blas::Daxpy(nquad, gmat[2][0], dPhysValuesdz.data(), 1,
                            dPhysValuesdx.data(), 1);
            }
        }
        break;
        default:
            ASSERTL0(false, "Wrong number of dimensions");
            break;
    }

    v_IProductWRTBase(m_base[0]->GetDbdata(), dPhysValuesdx, outarray, 1);
}

void SegExp::v_LaplacianMatrixOp(const int k1, const int k2,
                                 const Array<OneD, const NekDouble> &inarray,
                                 Array<OneD, NekDouble> &outarray,
                                 const StdRegions::StdMatrixKey &mkey)
{
    StdExpansion::LaplacianMatrixOp_MatFree(k1, k2, inarray, outarray, mkey);
}

void SegExp::v_HelmholtzMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                 Array<OneD, NekDouble> &outarray,
                                 const StdRegions::StdMatrixKey &mkey)
{
    int nquad = m_base[0]->GetNumPoints();
    const Array<TwoD, const NekDouble> &gmat =
        m_metricinfo->GetDerivFactors(GetPointsKeys());
    const NekDouble lambda = mkey.GetConstFactor(StdRegions::eFactorLambda);

    Array<OneD, NekDouble> physValues(nquad);
    Array<OneD, NekDouble> dPhysValuesdx(nquad);
    Array<OneD, NekDouble> wsp(m_ncoeffs);

    BwdTrans(inarray, physValues);

    // mass matrix operation
    v_IProductWRTBase((m_base[0]->GetBdata()), physValues, wsp, 1);

    // Laplacian matrix operation
    switch (m_geom->GetCoordim())
    {
        case 1:
        {
            PhysDeriv(physValues, dPhysValuesdx);

            // multiply with the proper geometric factors
            if (m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                Vmath::Vmul(nquad, &gmat[0][0], 1, dPhysValuesdx.data(), 1,
                            dPhysValuesdx.data(), 1);
            }
            else
            {
                Blas::Dscal(nquad, gmat[0][0], dPhysValuesdx.data(), 1);
            }
        }
        break;
        case 2:
        {
            Array<OneD, NekDouble> dPhysValuesdy(nquad);

            PhysDeriv(physValues, dPhysValuesdx, dPhysValuesdy);

            // multiply with the proper geometric factors
            if (m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                Vmath::Vmul(nquad, &gmat[0][0], 1, dPhysValuesdx.data(), 1,
                            dPhysValuesdx.data(), 1);
                Vmath::Vvtvp(nquad, &gmat[1][0], 1, dPhysValuesdy.data(), 1,
                             dPhysValuesdx.data(), 1, dPhysValuesdx.data(), 1);
            }
            else
            {
                Blas::Dscal(nquad, gmat[0][0], dPhysValuesdx.data(), 1);
                Blas::Daxpy(nquad, gmat[1][0], dPhysValuesdy.data(), 1,
                            dPhysValuesdx.data(), 1);
            }
        }
        break;
        case 3:
        {
            Array<OneD, NekDouble> dPhysValuesdy(nquad);
            Array<OneD, NekDouble> dPhysValuesdz(nquad);

            PhysDeriv(physValues, dPhysValuesdx, dPhysValuesdy, dPhysValuesdz);

            // multiply with the proper geometric factors
            if (m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                Vmath::Vmul(nquad, &gmat[0][0], 1, dPhysValuesdx.data(), 1,
                            dPhysValuesdx.data(), 1);
                Vmath::Vvtvp(nquad, &gmat[1][0], 1, dPhysValuesdy.data(), 1,
                             dPhysValuesdx.data(), 1, dPhysValuesdx.data(), 1);
                Vmath::Vvtvp(nquad, &gmat[2][0], 1, dPhysValuesdz.data(), 1,
                             dPhysValuesdx.data(), 1, dPhysValuesdx.data(), 1);
            }
            else
            {
                Blas::Dscal(nquad, gmat[0][0], dPhysValuesdx.data(), 1);
                Blas::Daxpy(nquad, gmat[1][0], dPhysValuesdy.data(), 1,
                            dPhysValuesdx.data(), 1);
                Blas::Daxpy(nquad, gmat[2][0], dPhysValuesdz.data(), 1,
                            dPhysValuesdx.data(), 1);
            }
        }
        break;
        default:
            ASSERTL0(false, "Wrong number of dimensions");
            break;
    }

    v_IProductWRTBase(m_base[0]->GetDbdata(), dPhysValuesdx, outarray, 1);
    Blas::Daxpy(m_ncoeffs, lambda, wsp.data(), 1, outarray.data(), 1);
}

//-----------------------------
// Matrix creation functions
//-----------------------------

DNekScalBlkMatSharedPtr SegExp::v_GetLocStaticCondMatrix(const MatrixKey &mkey)
{
    return m_staticCondMatrixManager[mkey];
}

void SegExp::v_DropLocStaticCondMatrix(const MatrixKey &mkey)
{
    m_staticCondMatrixManager.DeleteObject(mkey);
}

DNekScalMatSharedPtr SegExp::v_GetLocMatrix(const MatrixKey &mkey)
{
    return m_matrixManager[mkey];
}

void SegExp::v_DropLocMatrix(const MatrixKey &mkey)
{
    m_matrixManager.DeleteObject(mkey);
}

DNekMatSharedPtr SegExp::v_CreateStdMatrix(const StdRegions::StdMatrixKey &mkey)
{
    LibUtilities::BasisKey bkey = m_base[0]->GetBasisKey();
    StdRegions::StdSegExpSharedPtr tmp =
        MemoryManager<StdSegExp>::AllocateSharedPtr(bkey);

    return tmp->GetStdMatrix(mkey);
}

DNekScalMatSharedPtr SegExp::CreateMatrix(const MatrixKey &mkey)
{
    DNekScalMatSharedPtr returnval;
    NekDouble fac;
    LibUtilities::PointsKeyVector ptsKeys = GetPointsKeys();

    ASSERTL2(m_metricinfo->GetGtype() != SpatialDomains::eNoGeomType,
             "Geometric information is not set up");

    switch (mkey.GetMatrixType())
    {
        case StdRegions::eMass:
        {
            if ((m_metricinfo->GetGtype() == SpatialDomains::eDeformed) ||
                (mkey.GetNVarCoeff()))
            {
                fac = 1.0;
                goto UseLocRegionsMatrix;
            }
            else
            {
                fac = (m_metricinfo->GetJac(ptsKeys))[0];
                goto UseStdRegionsMatrix;
            }
        }
        break;
        case StdRegions::eInvMass:
        {
            if ((m_metricinfo->GetGtype() == SpatialDomains::eDeformed) ||
                (mkey.GetNVarCoeff()))
            {
                NekDouble one = 1.0;
                StdRegions::StdMatrixKey masskey(StdRegions::eMass,
                                                 DetShapeType(), *this);
                DNekMatSharedPtr mat = GenMatrix(masskey);
                mat->Invert();

                returnval =
                    MemoryManager<DNekScalMat>::AllocateSharedPtr(one, mat);
            }
            else
            {
                fac = 1.0 / (m_metricinfo->GetJac(ptsKeys))[0];
                goto UseStdRegionsMatrix;
            }
        }
        break;
        case StdRegions::eWeakDeriv0:
        case StdRegions::eWeakDeriv1:
        case StdRegions::eWeakDeriv2:
        {
            if (m_metricinfo->GetGtype() == SpatialDomains::eDeformed ||
                mkey.GetNVarCoeff())
            {
                fac = 1.0;
                goto UseLocRegionsMatrix;
            }
            else
            {
                int dir = 0;
                switch (mkey.GetMatrixType())
                {
                    case StdRegions::eWeakDeriv0:
                        dir = 0;
                        break;
                    case StdRegions::eWeakDeriv1:
                        ASSERTL1(m_geom->GetCoordim() >= 2,
                                 "Cannot call eWeakDeriv2 in a "
                                 "coordinate system which is not at "
                                 "least two-dimensional");
                        dir = 1;
                        break;
                    case StdRegions::eWeakDeriv2:
                        ASSERTL1(m_geom->GetCoordim() == 3,
                                 "Cannot call eWeakDeriv2 in a "
                                 "coordinate system which is not "
                                 "three-dimensional");
                        dir = 2;
                        break;
                    default:
                        break;
                }

                MatrixKey deriv0key(StdRegions::eWeakDeriv0,
                                    mkey.GetShapeType(), *this);

                DNekMatSharedPtr WeakDerivStd = GetStdMatrix(deriv0key);
                fac = m_metricinfo->GetDerivFactors(ptsKeys)[dir][0] *
                      m_metricinfo->GetJac(ptsKeys)[0];

                returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(
                    fac, WeakDerivStd);
            }
        }
        break;
        case StdRegions::eLaplacian:
        {
            if (m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                fac = 1.0;
                goto UseLocRegionsMatrix;
            }
            else
            {
                int coordim = m_geom->GetCoordim();
                fac         = 0.0;
                for (int i = 0; i < coordim; ++i)
                {
                    fac += m_metricinfo->GetDerivFactors(ptsKeys)[i][0] *
                           m_metricinfo->GetDerivFactors(ptsKeys)[i][0];
                }
                fac *= m_metricinfo->GetJac(ptsKeys)[0];
                goto UseStdRegionsMatrix;
            }
        }
        break;
        case StdRegions::eHelmholtz:
        {
            NekDouble factor = mkey.GetConstFactor(StdRegions::eFactorLambda);
            MatrixKey masskey(StdRegions::eMass, mkey.GetShapeType(), *this);
            DNekScalMat &MassMat = *(this->m_matrixManager[masskey]);
            MatrixKey lapkey(StdRegions::eLaplacian, mkey.GetShapeType(), *this,
                             mkey.GetConstFactors(), mkey.GetVarCoeffs());
            DNekScalMat &LapMat = *(this->m_matrixManager[lapkey]);

            int rows = LapMat.GetRows();
            int cols = LapMat.GetColumns();

            DNekMatSharedPtr helm =
                MemoryManager<DNekMat>::AllocateSharedPtr(rows, cols);

            NekDouble one = 1.0;
            (*helm)       = LapMat + factor * MassMat;

            returnval =
                MemoryManager<DNekScalMat>::AllocateSharedPtr(one, helm);
        }
        break;
        case StdRegions::eHybridDGHelmholtz:
        case StdRegions::eHybridDGLamToU:
        case StdRegions::eHybridDGLamToQ0:
        case StdRegions::eHybridDGHelmBndLam:
        {
            NekDouble one = 1.0;

            DNekMatSharedPtr mat = GenMatrix(mkey);
            returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one, mat);
        }
        break;
        case StdRegions::eInvHybridDGHelmholtz:
        {
            NekDouble one = 1.0;

            //                  StdRegions::StdMatrixKey
            //                  hkey(StdRegions::eHybridDGHelmholtz,
            //                                                DetShapeType(),*this,
            //                                                mkey.GetConstant(0),
            //                                                mkey.GetConstant(1));
            MatrixKey hkey(StdRegions::eHybridDGHelmholtz, DetShapeType(),
                           *this, mkey.GetConstFactors(), mkey.GetVarCoeffs());
            DNekMatSharedPtr mat = GenMatrix(hkey);

            mat->Invert();
            returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one, mat);
        }
        break;
        case StdRegions::eInterpGauss:
        {
            DNekMatSharedPtr m_Ix;
            Array<OneD, NekDouble> coords(1, 0.0);
            StdRegions::ConstFactorMap factors = mkey.GetConstFactors();
            int vertex = (int)factors[StdRegions::eFactorGaussVertex];

            coords[0] = (vertex == 0) ? -1.0 : 1.0;

            m_Ix = m_base[0]->GetI(coords);
            returnval =
                MemoryManager<DNekScalMat>::AllocateSharedPtr(1.0, m_Ix);
        }
        break;

        UseLocRegionsMatrix:
        {
            DNekMatSharedPtr mat = GenMatrix(mkey);
            returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(fac, mat);
        }
        break;
        UseStdRegionsMatrix:
        {
            DNekMatSharedPtr mat = GetStdMatrix(mkey);
            returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(fac, mat);
        }
        break;
        default:
            NEKERROR(ErrorUtil::efatal, "Matrix creation not defined");
            break;
    }

    return returnval;
}

DNekMatSharedPtr SegExp::v_GenMatrix(const StdRegions::StdMatrixKey &mkey)
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
            returnval = Expansion1D::v_GenMatrix(mkey);
            break;
        default:
            returnval = StdSegExp::v_GenMatrix(mkey);
            break;
    }

    return returnval;
}

//-----------------------------
// Private methods
//-----------------------------

///  Reverse the coefficients in a boundary interior expansion
///  this routine is of use when we need the segment
///  coefficients corresponding to a expansion in the reverse
///  coordinate direction
void SegExp::ReverseCoeffsAndSign(const Array<OneD, NekDouble> &inarray,
                                  Array<OneD, NekDouble> &outarray)
{

    int m;
    NekDouble sgn = 1;

    ASSERTL1(&inarray[0] != &outarray[0],
             "inarray and outarray can not be the same");
    switch (GetBasisType(0))
    {
        case LibUtilities::eModified_A:
            // Swap vertices
            outarray[0] = inarray[1];
            outarray[1] = inarray[0];
            // negate odd modes
            for (m = 2; m < m_ncoeffs; ++m)
            {
                outarray[m] = sgn * inarray[m];
                sgn         = -sgn;
            }
            break;
        case LibUtilities::eGLL_Lagrange:
        case LibUtilities::eGauss_Lagrange:
            for (m = 0; m < m_ncoeffs; ++m)
            {
                outarray[m_ncoeffs - 1 - m] = inarray[m];
            }
            break;
        default:
            ASSERTL0(false, "This basis is not allowed in this method");
            break;
    }
}

/* \brief Mass inversion product from \a inarray to \a outarray
 *
 *     Multiply by the inverse of the mass matrix
 *     \f$ {\bf \hat{u}} = {\bf M}^{-1} {\bf I} \f$
 *
 **/
void SegExp::MultiplyByElmtInvMass(const Array<OneD, const NekDouble> &inarray,
                                   Array<OneD, NekDouble> &outarray)
{
    // get Mass matrix inverse
    MatrixKey masskey(StdRegions::eInvMass, DetShapeType(), *this);
    DNekScalMatSharedPtr matsys = m_matrixManager[masskey];

    NekVector<NekDouble> in(m_ncoeffs, inarray, eCopy);
    NekVector<NekDouble> out(m_ncoeffs, outarray, eWrapper);

    out = (*matsys) * in;
}

} // namespace Nektar::LocalRegions
