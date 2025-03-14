///////////////////////////////////////////////////////////////////////////////
//
// File: IProductWRTBase.cpp
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
// Description: IProductWRTBase operator implementations
//
///////////////////////////////////////////////////////////////////////////////

#include <Collections/Collection.h>
#include <Collections/IProduct.h>
#include <Collections/MatrixFreeBase.h>
#include <Collections/Operator.h>
#include <MatrixFreeOps/Operator.hpp>

using namespace std;

namespace Nektar::Collections
{

using LibUtilities::eHexahedron;
using LibUtilities::ePrism;
using LibUtilities::ePyramid;
using LibUtilities::eQuadrilateral;
using LibUtilities::eSegment;
using LibUtilities::eTetrahedron;
using LibUtilities::eTriangle;

/**
 * @brief Inner product help class to calculate the size of the collection
 * that is given as an input and as an output to the IProductWRTBase Operator.
 * The size evaluation takes into account the conversion from the physical space
 * to the coefficient space.
 */
class IProductWRTBase_Helper : virtual public Operator
{
protected:
    IProductWRTBase_Helper()
    {
        // expect input to be number of elements by the number of quad points
        m_inputSize = m_numElmt * m_stdExp->GetTotPoints();
        // expect input to be number of elements by the number of coefficients
        m_outputSize = m_numElmt * m_stdExp->GetNcoeffs();
    }
};

/**
 * @brief Inner product operator using standard matrix approach
 */
class IProductWRTBase_StdMat final : virtual public Operator,
                                     virtual public IProductWRTBase_Helper
{
public:
    OPERATOR_CREATE(IProductWRTBase_StdMat)

    ~IProductWRTBase_StdMat() final = default;

    void operator()(const Array<OneD, const NekDouble> &input,
                    Array<OneD, NekDouble> &output,
                    [[maybe_unused]] Array<OneD, NekDouble> &output1,
                    [[maybe_unused]] Array<OneD, NekDouble> &output2,
                    Array<OneD, NekDouble> &wsp) final
    {
        ASSERTL1(wsp.size() == m_wspSize, "Incorrect workspace size");

        if (m_isDeformed)
        {
            Vmath::Vmul(m_jac.size(), m_jac, 1, input, 1, wsp, 1);
        }
        else
        {
            Array<OneD, NekDouble> tmp;
            for (int e = 0; e < m_numElmt; ++e)
            {
                Vmath::Smul(m_nqe, m_jac[e], input + e * m_nqe, 1,
                            tmp = wsp + e * m_nqe, 1);
            }
        }

        Blas::Dgemm('N', 'N', m_mat->GetRows(), m_numElmt, m_mat->GetColumns(),
                    1.0, m_mat->GetRawPtr(), m_mat->GetRows(), wsp.data(),
                    m_stdExp->GetTotPoints(), 0.0, output.data(),
                    m_stdExp->GetNcoeffs());
    }

    void operator()([[maybe_unused]] int dir,
                    [[maybe_unused]] const Array<OneD, const NekDouble> &input,
                    [[maybe_unused]] Array<OneD, NekDouble> &output,
                    [[maybe_unused]] Array<OneD, NekDouble> &wsp) final
    {
        NEKERROR(ErrorUtil::efatal, "Not valid for this operator.");
    }

protected:
    DNekMatSharedPtr m_mat;
    Array<OneD, const NekDouble> m_jac;

private:
    IProductWRTBase_StdMat(vector<StdRegions::StdExpansionSharedPtr> pCollExp,
                           CoalescedGeomDataSharedPtr pGeomData,
                           StdRegions::FactorMap factors)
        : Operator(pCollExp, pGeomData, factors), IProductWRTBase_Helper()
    {
        m_jac = pGeomData->GetJac(pCollExp);
        StdRegions::StdMatrixKey key(StdRegions::eIProductWRTBase,
                                     m_stdExp->DetShapeType(), *m_stdExp);
        m_mat     = m_stdExp->GetStdMatrix(key);
        m_nqe     = m_stdExp->GetTotPoints();
        m_wspSize = m_nqe * m_numElmt;
    }
};

/// Factory initialisation for the IProductWRTBase_StdMat operators
OperatorKey IProductWRTBase_StdMat::m_typeArr[] = {
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eSegment, eIProductWRTBase, eStdMat, false),
        IProductWRTBase_StdMat::create, "IProductWRTBase_StdMat_Seg"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTriangle, eIProductWRTBase, eStdMat, false),
        IProductWRTBase_StdMat::create, "IProductWRTBase_StdMat_Tri"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTriangle, eIProductWRTBase, eStdMat, true),
        IProductWRTBase_StdMat::create, "IProductWRTBase_StdMat_NodalTri"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eQuadrilateral, eIProductWRTBase, eStdMat, false),
        IProductWRTBase_StdMat::create, "IProductWRTBase_StdMat_Quad"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTetrahedron, eIProductWRTBase, eStdMat, false),
        IProductWRTBase_StdMat::create, "IProductWRTBase_StdMat_Tet"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTetrahedron, eIProductWRTBase, eStdMat, true),
        IProductWRTBase_StdMat::create, "IProductWRTBase_StdMat_NodalTet"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePyramid, eIProductWRTBase, eStdMat, false),
        IProductWRTBase_StdMat::create, "IProductWRTBase_StdMat_Pyr"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePrism, eIProductWRTBase, eStdMat, false),
        IProductWRTBase_StdMat::create, "IProductWRTBase_StdMat_Prism"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePrism, eIProductWRTBase, eStdMat, true),
        IProductWRTBase_StdMat::create, "IProductWRTBase_StdMat_NodalPrism"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eHexahedron, eIProductWRTBase, eStdMat, false),
        IProductWRTBase_StdMat::create, "IProductWRTBase_StdMat_Hex"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePyramid, eIProductWRTBase, eSumFac, false),
        IProductWRTBase_StdMat::create, "IProductWRTBase_SumFac_Pyr")};

/**
 * @brief Inner product operator using operator using matrix free operators.
 */
class IProductWRTBase_MatrixFree final : virtual public Operator,
                                         MatrixFreeBase,
                                         virtual public IProductWRTBase_Helper
{
public:
    OPERATOR_CREATE(IProductWRTBase_MatrixFree)

    ~IProductWRTBase_MatrixFree() final = default;

    void operator()(const Array<OneD, const NekDouble> &input,
                    Array<OneD, NekDouble> &output,
                    [[maybe_unused]] Array<OneD, NekDouble> &output1,
                    [[maybe_unused]] Array<OneD, NekDouble> &output2,
                    [[maybe_unused]] Array<OneD, NekDouble> &wsp) final
    {
        (*m_oper)(input, output);
    }

    void operator()([[maybe_unused]] int dir,
                    [[maybe_unused]] const Array<OneD, const NekDouble> &input,
                    [[maybe_unused]] Array<OneD, NekDouble> &output,
                    [[maybe_unused]] Array<OneD, NekDouble> &wsp) final
    {
        NEKERROR(ErrorUtil::efatal, "Not valid for this operator.");
    }

private:
    std::shared_ptr<MatrixFree::IProduct> m_oper;

    IProductWRTBase_MatrixFree(
        vector<StdRegions::StdExpansionSharedPtr> pCollExp,
        CoalescedGeomDataSharedPtr pGeomData, StdRegions::FactorMap factors)
        : Operator(pCollExp, pGeomData, factors), IProductWRTBase_Helper(),
          MatrixFreeBase(pCollExp[0]->GetStdExp()->GetTotPoints(),
                         pCollExp[0]->GetStdExp()->GetNcoeffs(),
                         pCollExp.size())
    {
        // Basis vector
        const auto dim = pCollExp[0]->GetStdExp()->GetShapeDimension();
        std::vector<LibUtilities::BasisSharedPtr> basis(dim);
        for (unsigned int i = 0; i < dim; ++i)
        {
            basis[i] = pCollExp[0]->GetBasis(i);
        }

        // Get shape type
        auto shapeType = pCollExp[0]->GetStdExp()->DetShapeType();

        // Generate operator string and create operator.
        std::string op_string = "IProduct";
        op_string += MatrixFree::GetOpstring(shapeType, m_isDeformed);
        auto oper = MatrixFree::GetOperatorFactory().CreateInstance(
            op_string, basis, pCollExp.size());

        oper->SetUpBdata(basis);
        oper->SetUpZW(basis);

        // Set Jacobian
        oper->SetJac(pGeomData->GetJacInterLeave(pCollExp, m_nElmtPad));

        m_oper = std::dynamic_pointer_cast<MatrixFree::IProduct>(oper);
        ASSERTL0(m_oper, "Failed to cast pointer.");
    }
};

/// Factory initialisation for the IProductWRTBase_MatrixFree operators
OperatorKey IProductWRTBase_MatrixFree::m_typeArr[] = {
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eSegment, eIProductWRTBase, eMatrixFree, false),
        IProductWRTBase_MatrixFree::create, "IProductWRTBase_MatrixFree_Seg"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eQuadrilateral, eIProductWRTBase, eMatrixFree, false),
        IProductWRTBase_MatrixFree::create, "IProductWRTBase_MatrixFree_Quad"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTriangle, eIProductWRTBase, eMatrixFree, false),
        IProductWRTBase_MatrixFree::create, "IProductWRTBase_MatrixFree_Tri"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eHexahedron, eIProductWRTBase, eMatrixFree, false),
        IProductWRTBase_MatrixFree::create, "IProductWRTBase_MatrixFree_Hex"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePrism, eIProductWRTBase, eMatrixFree, false),
        IProductWRTBase_MatrixFree::create, "IProductWRTBase_MatrixFree_Prism"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePyramid, eIProductWRTBase, eMatrixFree, false),
        IProductWRTBase_MatrixFree::create, "IProductWRTBase_MatrixFree_Pyr"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTetrahedron, eIProductWRTBase, eMatrixFree, false),
        IProductWRTBase_MatrixFree::create, "IProductWRTBase_MatrixFree_Tet")};

/**
 * @brief Inner product operator using element-wise operation
 */
class IProductWRTBase_IterPerExp final : virtual public Operator,
                                         virtual public IProductWRTBase_Helper
{
public:
    OPERATOR_CREATE(IProductWRTBase_IterPerExp)

    ~IProductWRTBase_IterPerExp() final = default;

    void operator()(const Array<OneD, const NekDouble> &input,
                    Array<OneD, NekDouble> &output,
                    [[maybe_unused]] Array<OneD, NekDouble> &output1,
                    [[maybe_unused]] Array<OneD, NekDouble> &output2,
                    Array<OneD, NekDouble> &wsp) final
    {
        ASSERTL1(wsp.size() == m_wspSize, "Incorrect workspace size");

        const int nCoeffs = m_stdExp->GetNcoeffs();
        const int nPhys   = m_stdExp->GetTotPoints();
        Array<OneD, NekDouble> tmp;

        Vmath::Vmul(m_jacWStdW.size(), m_jacWStdW, 1, input, 1, wsp, 1);

        for (int i = 0; i < m_numElmt; ++i)
        {
            m_stdExp->IProductWRTBase_SumFac(wsp + i * nPhys,
                                             tmp = output + i * nCoeffs, false);
        }
    }

    void operator()([[maybe_unused]] int dir,
                    [[maybe_unused]] const Array<OneD, const NekDouble> &input,
                    [[maybe_unused]] Array<OneD, NekDouble> &output,
                    [[maybe_unused]] Array<OneD, NekDouble> &wsp) final
    {
        NEKERROR(ErrorUtil::efatal, "Not valid for this operator.");
    }

protected:
    Array<OneD, NekDouble> m_jacWStdW;

private:
    IProductWRTBase_IterPerExp(
        vector<StdRegions::StdExpansionSharedPtr> pCollExp,
        CoalescedGeomDataSharedPtr pGeomData, StdRegions::FactorMap factors)
        : Operator(pCollExp, pGeomData, factors), IProductWRTBase_Helper()
    {
        int nqtot = pCollExp[0]->GetTotPoints();

        m_jacWStdW = pGeomData->GetJacWithStdWeights(pCollExp);

        m_wspSize = nqtot * m_numElmt;
    }
};

/// Factory initialisation for the IProductWRTBase_IterPerExp operators
OperatorKey IProductWRTBase_IterPerExp::m_typeArr[] = {
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eSegment, eIProductWRTBase, eIterPerExp, false),
        IProductWRTBase_IterPerExp::create, "IProductWRTBase_IterPerExp_Seg"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTriangle, eIProductWRTBase, eIterPerExp, false),
        IProductWRTBase_IterPerExp::create, "IProductWRTBase_IterPerExp_Tri"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTriangle, eIProductWRTBase, eIterPerExp, true),
        IProductWRTBase_IterPerExp::create,
        "IProductWRTBase_IterPerExp_NodalTri"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eQuadrilateral, eIProductWRTBase, eIterPerExp, false),
        IProductWRTBase_IterPerExp::create, "IProductWRTBase_IterPerExp_Quad"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTetrahedron, eIProductWRTBase, eIterPerExp, false),
        IProductWRTBase_IterPerExp::create, "IProductWRTBase_IterPerExp_Tet"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTetrahedron, eIProductWRTBase, eIterPerExp, true),
        IProductWRTBase_IterPerExp::create,
        "IProductWRTBase_IterPerExp_NodalTet"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePyramid, eIProductWRTBase, eIterPerExp, false),
        IProductWRTBase_IterPerExp::create, "IProductWRTBase_IterPerExp_Pyr"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePrism, eIProductWRTBase, eIterPerExp, false),
        IProductWRTBase_IterPerExp::create, "IProductWRTBase_IterPerExp_Prism"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePrism, eIProductWRTBase, eIterPerExp, true),
        IProductWRTBase_IterPerExp::create,
        "IProductWRTBase_IterPerExp_NodalPrism"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eHexahedron, eIProductWRTBase, eIterPerExp, false),
        IProductWRTBase_IterPerExp::create, "IProductWRTBase_IterPerExp_Hex"),
};

/**
 * @brief Inner product operator using original MultiRegions implementation.
 */
class IProductWRTBase_NoCollection final : virtual public Operator,
                                           virtual public IProductWRTBase_Helper
{
public:
    OPERATOR_CREATE(IProductWRTBase_NoCollection)

    ~IProductWRTBase_NoCollection() final = default;

    void operator()(const Array<OneD, const NekDouble> &input,
                    Array<OneD, NekDouble> &output,
                    [[maybe_unused]] Array<OneD, NekDouble> &output1,
                    [[maybe_unused]] Array<OneD, NekDouble> &output2,
                    [[maybe_unused]] Array<OneD, NekDouble> &wsp) final
    {
        const int nCoeffs = m_expList[0]->GetNcoeffs();
        const int nPhys   = m_expList[0]->GetTotPoints();
        Array<OneD, NekDouble> tmp;

        for (int i = 0; i < m_numElmt; ++i)
        {
            m_expList[i]->IProductWRTBase(input + i * nPhys,
                                          tmp = output + i * nCoeffs);
        }
    }

    void operator()([[maybe_unused]] int dir,
                    [[maybe_unused]] const Array<OneD, const NekDouble> &input,
                    [[maybe_unused]] Array<OneD, NekDouble> &output,
                    [[maybe_unused]] Array<OneD, NekDouble> &wsp) final
    {
        NEKERROR(ErrorUtil::efatal, "Not valid for this operator.");
    }

protected:
    vector<StdRegions::StdExpansionSharedPtr> m_expList;

private:
    IProductWRTBase_NoCollection(
        vector<StdRegions::StdExpansionSharedPtr> pCollExp,
        CoalescedGeomDataSharedPtr pGeomData, StdRegions::FactorMap factors)
        : Operator(pCollExp, pGeomData, factors), IProductWRTBase_Helper()
    {
        m_expList = pCollExp;
    }
};

/// Factory initialisation for the IProductWRTBase_NoCollection operators
OperatorKey IProductWRTBase_NoCollection::m_typeArr[] = {
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eSegment, eIProductWRTBase, eNoCollection, false),
        IProductWRTBase_NoCollection::create,
        "IProductWRTBase_NoCollection_Seg"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTriangle, eIProductWRTBase, eNoCollection, false),
        IProductWRTBase_NoCollection::create,
        "IProductWRTBase_NoCollection_Tri"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTriangle, eIProductWRTBase, eNoCollection, true),
        IProductWRTBase_NoCollection::create,
        "IProductWRTBase_NoCollection_NodalTri"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eQuadrilateral, eIProductWRTBase, eNoCollection, false),
        IProductWRTBase_NoCollection::create,
        "IProductWRTBase_NoCollection_Quad"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTetrahedron, eIProductWRTBase, eNoCollection, false),
        IProductWRTBase_NoCollection::create,
        "IProductWRTBase_NoCollection_Tet"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTetrahedron, eIProductWRTBase, eNoCollection, true),
        IProductWRTBase_NoCollection::create,
        "IProductWRTBase_NoCollection_NodalTet"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePyramid, eIProductWRTBase, eNoCollection, false),
        IProductWRTBase_NoCollection::create,
        "IProductWRTBase_NoCollection_Pyr"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePrism, eIProductWRTBase, eNoCollection, false),
        IProductWRTBase_NoCollection::create,
        "IProductWRTBase_NoCollection_Prism"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePrism, eIProductWRTBase, eNoCollection, true),
        IProductWRTBase_NoCollection::create,
        "IProductWRTBase_NoCollection_NodalPrism"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eHexahedron, eIProductWRTBase, eNoCollection, false),
        IProductWRTBase_NoCollection::create,
        "IProductWRTBase_NoCollection_Hex"),
};

/**
 * @brief Inner product operator using sum-factorisation (Segment)
 */
class IProductWRTBase_SumFac_Seg final : virtual public Operator,
                                         virtual public IProductWRTBase_Helper
{
public:
    OPERATOR_CREATE(IProductWRTBase_SumFac_Seg)

    ~IProductWRTBase_SumFac_Seg() final = default;

    void operator()(const Array<OneD, const NekDouble> &input,
                    Array<OneD, NekDouble> &output,
                    [[maybe_unused]] Array<OneD, NekDouble> &output1,
                    [[maybe_unused]] Array<OneD, NekDouble> &output2,
                    Array<OneD, NekDouble> &wsp) final
    {
        if (m_colldir0)
        {
            Vmath::Vmul(m_numElmt * m_nquad0, m_jacWStdW, 1, input, 1, output,
                        1);
        }
        else
        {
            Vmath::Vmul(m_numElmt * m_nquad0, m_jacWStdW, 1, input, 1, wsp, 1);

            // out = B0*in;
            Blas::Dgemm('T', 'N', m_nmodes0, m_numElmt, m_nquad0, 1.0,
                        m_base0.data(), m_nquad0, &wsp[0], m_nquad0, 0.0,
                        &output[0], m_nmodes0);
        }
    }

    void operator()([[maybe_unused]] int dir,
                    [[maybe_unused]] const Array<OneD, const NekDouble> &input,
                    [[maybe_unused]] Array<OneD, NekDouble> &output,
                    [[maybe_unused]] Array<OneD, NekDouble> &wsp) final
    {
        NEKERROR(ErrorUtil::efatal, "Not valid for this operator.");
    }

protected:
    const int m_nquad0;
    const int m_nmodes0;
    const bool m_colldir0;
    Array<OneD, const NekDouble> m_jacWStdW;
    Array<OneD, const NekDouble> m_base0;

private:
    IProductWRTBase_SumFac_Seg(
        vector<StdRegions::StdExpansionSharedPtr> pCollExp,
        CoalescedGeomDataSharedPtr pGeomData, StdRegions::FactorMap factors)
        : Operator(pCollExp, pGeomData, factors), IProductWRTBase_Helper(),
          m_nquad0(m_stdExp->GetNumPoints(0)),
          m_nmodes0(m_stdExp->GetBasisNumModes(0)),
          m_colldir0(m_stdExp->GetBasis(0)->Collocation()),
          m_base0(m_stdExp->GetBasis(0)->GetBdata())
    {
        m_wspSize  = m_numElmt * m_nquad0;
        m_jacWStdW = pGeomData->GetJacWithStdWeights(pCollExp);
    }
};

/// Factory initialisation for the IProductWRTBase_SumFac_Seg operator
OperatorKey IProductWRTBase_SumFac_Seg::m_type =
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eSegment, eIProductWRTBase, eSumFac, false),
        IProductWRTBase_SumFac_Seg::create, "IProductWRTBase_SumFac_Seg");

/**
 * @brief Inner product operator using sum-factorisation (Quad)
 */
class IProductWRTBase_SumFac_Quad final : virtual public Operator,
                                          virtual public IProductWRTBase_Helper
{
public:
    OPERATOR_CREATE(IProductWRTBase_SumFac_Quad)

    ~IProductWRTBase_SumFac_Quad() final = default;

    void operator()(const Array<OneD, const NekDouble> &input,
                    Array<OneD, NekDouble> &output,
                    [[maybe_unused]] Array<OneD, NekDouble> &output1,
                    [[maybe_unused]] Array<OneD, NekDouble> &output2,
                    Array<OneD, NekDouble> &wsp) final
    {
        ASSERTL1(wsp.size() == m_wspSize, "Incorrect workspace size");

        QuadIProduct(m_colldir0, m_colldir1, m_numElmt, m_nquad0, m_nquad1,
                     m_nmodes0, m_nmodes1, m_base0, m_base1, m_jacWStdW, input,
                     output, wsp);
    }

    void operator()([[maybe_unused]] int dir,
                    [[maybe_unused]] const Array<OneD, const NekDouble> &input,
                    [[maybe_unused]] Array<OneD, NekDouble> &output,
                    [[maybe_unused]] Array<OneD, NekDouble> &wsp) final
    {
        NEKERROR(ErrorUtil::efatal, "Not valid for this operator.");
    }

protected:
    const int m_nquad0;
    const int m_nquad1;
    const int m_nmodes0;
    const int m_nmodes1;
    const bool m_colldir0;
    const bool m_colldir1;
    Array<OneD, const NekDouble> m_jacWStdW;
    Array<OneD, const NekDouble> m_base0;
    Array<OneD, const NekDouble> m_base1;

private:
    IProductWRTBase_SumFac_Quad(
        vector<StdRegions::StdExpansionSharedPtr> pCollExp,
        CoalescedGeomDataSharedPtr pGeomData, StdRegions::FactorMap factors)
        : Operator(pCollExp, pGeomData, factors), IProductWRTBase_Helper(),
          m_nquad0(m_stdExp->GetNumPoints(0)),
          m_nquad1(m_stdExp->GetNumPoints(1)),
          m_nmodes0(m_stdExp->GetBasisNumModes(0)),
          m_nmodes1(m_stdExp->GetBasisNumModes(1)),
          m_colldir0(m_stdExp->GetBasis(0)->Collocation()),
          m_colldir1(m_stdExp->GetBasis(1)->Collocation()),
          m_base0(m_stdExp->GetBasis(0)->GetBdata()),
          m_base1(m_stdExp->GetBasis(1)->GetBdata())
    {
        m_jacWStdW = pGeomData->GetJacWithStdWeights(pCollExp);
        m_wspSize =
            2 * m_numElmt * (max(m_nquad0 * m_nquad1, m_nmodes0 * m_nmodes1));
    }
};

/// Factory initialisation for the IProductWRTBase_SumFac_Quad operator
OperatorKey IProductWRTBase_SumFac_Quad::m_type =
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eQuadrilateral, eIProductWRTBase, eSumFac, false),
        IProductWRTBase_SumFac_Quad::create, "IProductWRTBase_SumFac_Quad");

/**
 * @brief Inner product operator using sum-factorisation (Tri)
 */
class IProductWRTBase_SumFac_Tri final : virtual public Operator,
                                         virtual public IProductWRTBase_Helper
{
public:
    OPERATOR_CREATE(IProductWRTBase_SumFac_Tri)

    ~IProductWRTBase_SumFac_Tri() final = default;

    void operator()(const Array<OneD, const NekDouble> &input,
                    Array<OneD, NekDouble> &output,
                    [[maybe_unused]] Array<OneD, NekDouble> &output1,
                    [[maybe_unused]] Array<OneD, NekDouble> &output2,
                    Array<OneD, NekDouble> &wsp) final
    {
        ASSERTL1(wsp.size() == m_wspSize, "Incorrect workspace size");

        TriIProduct(m_sortTopVertex, m_numElmt, m_nquad0, m_nquad1, m_nmodes0,
                    m_nmodes1, m_base0, m_base1, m_jacWStdW, input, output,
                    wsp);
    }

    void operator()([[maybe_unused]] int dir,
                    [[maybe_unused]] const Array<OneD, const NekDouble> &input,
                    [[maybe_unused]] Array<OneD, NekDouble> &output,
                    [[maybe_unused]] Array<OneD, NekDouble> &wsp) final
    {
        NEKERROR(ErrorUtil::efatal, "Not valid for this operator.");
    }

protected:
    const int m_nquad0;
    const int m_nquad1;
    const int m_nmodes0;
    const int m_nmodes1;
    Array<OneD, const NekDouble> m_jacWStdW;
    Array<OneD, const NekDouble> m_base0;
    Array<OneD, const NekDouble> m_base1;
    bool m_sortTopVertex;

private:
    IProductWRTBase_SumFac_Tri(
        vector<StdRegions::StdExpansionSharedPtr> pCollExp,
        CoalescedGeomDataSharedPtr pGeomData, StdRegions::FactorMap factors)
        : Operator(pCollExp, pGeomData, factors), IProductWRTBase_Helper(),
          m_nquad0(m_stdExp->GetNumPoints(0)),
          m_nquad1(m_stdExp->GetNumPoints(1)),
          m_nmodes0(m_stdExp->GetBasisNumModes(0)),
          m_nmodes1(m_stdExp->GetBasisNumModes(1)),
          m_base0(m_stdExp->GetBasis(0)->GetBdata()),
          m_base1(m_stdExp->GetBasis(1)->GetBdata())
    {
        m_jacWStdW = pGeomData->GetJacWithStdWeights(pCollExp);
        m_wspSize =
            2 * m_numElmt * (max(m_nquad0 * m_nquad1, m_nmodes0 * m_nmodes1));
        if (m_stdExp->GetBasis(0)->GetBasisType() == LibUtilities::eModified_A)
        {
            m_sortTopVertex = true;
        }
        else
        {
            m_sortTopVertex = false;
        }
    }
};

/// Factory initialisation for the IProductWRTBase_SumFac_Tri operator
OperatorKey IProductWRTBase_SumFac_Tri::m_type =
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTriangle, eIProductWRTBase, eSumFac, false),
        IProductWRTBase_SumFac_Tri::create, "IProductWRTBase_SumFac_Tri");

/**
 * @brief Inner Product operator using sum-factorisation (Hex)
 */
class IProductWRTBase_SumFac_Hex final : virtual public Operator,
                                         virtual public IProductWRTBase_Helper
{
public:
    OPERATOR_CREATE(IProductWRTBase_SumFac_Hex)

    ~IProductWRTBase_SumFac_Hex() final = default;

    void operator()(const Array<OneD, const NekDouble> &input,
                    Array<OneD, NekDouble> &output,
                    [[maybe_unused]] Array<OneD, NekDouble> &output1,
                    [[maybe_unused]] Array<OneD, NekDouble> &output2,
                    Array<OneD, NekDouble> &wsp) final
    {
        ASSERTL1(wsp.size() == m_wspSize, "Incorrect workspace size");

        HexIProduct(m_colldir0, m_colldir1, m_colldir2, m_numElmt, m_nquad0,
                    m_nquad1, m_nquad2, m_nmodes0, m_nmodes1, m_nmodes2,
                    m_base0, m_base1, m_base2, m_jacWStdW, input, output, wsp);
    }

    void operator()([[maybe_unused]] int dir,
                    [[maybe_unused]] const Array<OneD, const NekDouble> &input,
                    [[maybe_unused]] Array<OneD, NekDouble> &output,
                    [[maybe_unused]] Array<OneD, NekDouble> &wsp) final
    {
        NEKERROR(ErrorUtil::efatal, "Not valid for this operator.");
    }

protected:
    const int m_nquad0;
    const int m_nquad1;
    const int m_nquad2;
    const int m_nmodes0;
    const int m_nmodes1;
    const int m_nmodes2;
    const bool m_colldir0;
    const bool m_colldir1;
    const bool m_colldir2;
    Array<OneD, const NekDouble> m_jacWStdW;
    Array<OneD, const NekDouble> m_base0;
    Array<OneD, const NekDouble> m_base1;
    Array<OneD, const NekDouble> m_base2;

private:
    IProductWRTBase_SumFac_Hex(
        vector<StdRegions::StdExpansionSharedPtr> pCollExp,
        CoalescedGeomDataSharedPtr pGeomData, StdRegions::FactorMap factors)
        : Operator(pCollExp, pGeomData, factors), IProductWRTBase_Helper(),
          m_nquad0(m_stdExp->GetNumPoints(0)),
          m_nquad1(m_stdExp->GetNumPoints(1)),
          m_nquad2(m_stdExp->GetNumPoints(2)),
          m_nmodes0(m_stdExp->GetBasisNumModes(0)),
          m_nmodes1(m_stdExp->GetBasisNumModes(1)),
          m_nmodes2(m_stdExp->GetBasisNumModes(2)),
          m_colldir0(m_stdExp->GetBasis(0)->Collocation()),
          m_colldir1(m_stdExp->GetBasis(1)->Collocation()),
          m_colldir2(m_stdExp->GetBasis(2)->Collocation()),
          m_base0(m_stdExp->GetBasis(0)->GetBdata()),
          m_base1(m_stdExp->GetBasis(1)->GetBdata()),
          m_base2(m_stdExp->GetBasis(2)->GetBdata())

    {
        m_jacWStdW = pGeomData->GetJacWithStdWeights(pCollExp);
        m_wspSize  = 3 * m_numElmt *
                    (max(m_nquad0 * m_nquad1 * m_nquad2,
                         m_nmodes0 * m_nmodes1 * m_nmodes2));
    }
};

/// Factory initialisation for the IProductWRTBase_SumFac_Hex operator
OperatorKey IProductWRTBase_SumFac_Hex::m_type =
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eHexahedron, eIProductWRTBase, eSumFac, false),
        IProductWRTBase_SumFac_Hex::create, "IProductWRTBase_SumFac_Hex");

/**
 * @brief Inner product operator using sum-factorisation (Tet)
 */
class IProductWRTBase_SumFac_Tet final : virtual public Operator,
                                         virtual public IProductWRTBase_Helper
{
public:
    OPERATOR_CREATE(IProductWRTBase_SumFac_Tet)

    ~IProductWRTBase_SumFac_Tet() final = default;

    void operator()(const Array<OneD, const NekDouble> &input,
                    Array<OneD, NekDouble> &output,
                    [[maybe_unused]] Array<OneD, NekDouble> &output1,
                    [[maybe_unused]] Array<OneD, NekDouble> &output2,
                    Array<OneD, NekDouble> &wsp) final
    {
        ASSERTL1(wsp.size() == m_wspSize, "Incorrect workspace size");

        TetIProduct(m_sortTopEdge, m_numElmt, m_nquad0, m_nquad1, m_nquad2,
                    m_nmodes0, m_nmodes1, m_nmodes2, m_base0, m_base1, m_base2,
                    m_jacWStdW, input, output, wsp);
    }

    void operator()([[maybe_unused]] int dir,
                    [[maybe_unused]] const Array<OneD, const NekDouble> &input,
                    [[maybe_unused]] Array<OneD, NekDouble> &output,
                    [[maybe_unused]] Array<OneD, NekDouble> &wsp) final
    {
        NEKERROR(ErrorUtil::efatal, "Not valid for this operator.");
    }

protected:
    const int m_nquad0;
    const int m_nquad1;
    const int m_nquad2;
    const int m_nmodes0;
    const int m_nmodes1;
    const int m_nmodes2;
    Array<OneD, const NekDouble> m_jacWStdW;
    Array<OneD, const NekDouble> m_base0;
    Array<OneD, const NekDouble> m_base1;
    Array<OneD, const NekDouble> m_base2;
    bool m_sortTopEdge;

private:
    IProductWRTBase_SumFac_Tet(
        vector<StdRegions::StdExpansionSharedPtr> pCollExp,
        CoalescedGeomDataSharedPtr pGeomData, StdRegions::FactorMap factors)
        : Operator(pCollExp, pGeomData, factors), IProductWRTBase_Helper(),
          m_nquad0(m_stdExp->GetNumPoints(0)),
          m_nquad1(m_stdExp->GetNumPoints(1)),
          m_nquad2(m_stdExp->GetNumPoints(2)),
          m_nmodes0(m_stdExp->GetBasisNumModes(0)),
          m_nmodes1(m_stdExp->GetBasisNumModes(1)),
          m_nmodes2(m_stdExp->GetBasisNumModes(2)),
          m_base0(m_stdExp->GetBasis(0)->GetBdata()),
          m_base1(m_stdExp->GetBasis(1)->GetBdata()),
          m_base2(m_stdExp->GetBasis(2)->GetBdata())
    {
        m_jacWStdW = pGeomData->GetJacWithStdWeights(pCollExp);
        m_wspSize =
            m_numElmt *
            (max(m_nquad0 * m_nquad1 * m_nquad2,
                 m_nquad2 * m_nmodes0 * (2 * m_nmodes1 - m_nmodes0 + 1) / 2) +
             m_nquad2 * m_nquad1 * m_nmodes0);

        if (m_stdExp->GetBasis(0)->GetBasisType() == LibUtilities::eModified_A)
        {
            m_sortTopEdge = true;
        }
        else
        {
            m_sortTopEdge = false;
        }
    }
};

/// Factory initialisation for the IProductWRTBase_SumFac_Tet operator
OperatorKey IProductWRTBase_SumFac_Tet::m_type =
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTetrahedron, eIProductWRTBase, eSumFac, false),
        IProductWRTBase_SumFac_Tet::create, "IProductWRTBase_SumFac_Tet");

/**
 * @brief Inner Product operator using sum-factorisation (Prism)
 */
class IProductWRTBase_SumFac_Prism final : virtual public Operator,
                                           IProductWRTBase_Helper
{
public:
    OPERATOR_CREATE(IProductWRTBase_SumFac_Prism)

    ~IProductWRTBase_SumFac_Prism() final = default;

    void operator()(const Array<OneD, const NekDouble> &input,
                    Array<OneD, NekDouble> &output,
                    [[maybe_unused]] Array<OneD, NekDouble> &output1,
                    [[maybe_unused]] Array<OneD, NekDouble> &output2,
                    Array<OneD, NekDouble> &wsp) final
    {
        ASSERTL1(wsp.size() == m_wspSize, "Incorrect workspace size");

        PrismIProduct(m_sortTopVertex, m_numElmt, m_nquad0, m_nquad1, m_nquad2,
                      m_nmodes0, m_nmodes1, m_nmodes2, m_base0, m_base1,
                      m_base2, m_jacWStdW, input, output, wsp);
    }

    void operator()([[maybe_unused]] int dir,
                    [[maybe_unused]] const Array<OneD, const NekDouble> &input,
                    [[maybe_unused]] Array<OneD, NekDouble> &output,
                    [[maybe_unused]] Array<OneD, NekDouble> &wsp) final
    {
        NEKERROR(ErrorUtil::efatal, "Not valid for this operator.");
    }

protected:
    const int m_nquad0;
    const int m_nquad1;
    const int m_nquad2;
    const int m_nmodes0;
    const int m_nmodes1;
    const int m_nmodes2;
    Array<OneD, const NekDouble> m_jacWStdW;
    Array<OneD, const NekDouble> m_base0;
    Array<OneD, const NekDouble> m_base1;
    Array<OneD, const NekDouble> m_base2;
    bool m_sortTopVertex;

private:
    IProductWRTBase_SumFac_Prism(
        vector<StdRegions::StdExpansionSharedPtr> pCollExp,
        CoalescedGeomDataSharedPtr pGeomData, StdRegions::FactorMap factors)
        : Operator(pCollExp, pGeomData, factors), IProductWRTBase_Helper(),
          m_nquad0(m_stdExp->GetNumPoints(0)),
          m_nquad1(m_stdExp->GetNumPoints(1)),
          m_nquad2(m_stdExp->GetNumPoints(2)),
          m_nmodes0(m_stdExp->GetBasisNumModes(0)),
          m_nmodes1(m_stdExp->GetBasisNumModes(1)),
          m_nmodes2(m_stdExp->GetBasisNumModes(2)),
          m_base0(m_stdExp->GetBasis(0)->GetBdata()),
          m_base1(m_stdExp->GetBasis(1)->GetBdata()),
          m_base2(m_stdExp->GetBasis(2)->GetBdata())

    {
        m_jacWStdW = pGeomData->GetJacWithStdWeights(pCollExp);

        m_wspSize = m_numElmt * m_nquad2 *
                        (max(m_nquad0 * m_nquad1, m_nmodes0 * m_nmodes1)) +
                    m_nquad1 * m_nquad2 * m_numElmt * m_nmodes0;

        if (m_stdExp->GetBasis(0)->GetBasisType() == LibUtilities::eModified_A)
        {
            m_sortTopVertex = true;
        }
        else
        {
            m_sortTopVertex = false;
        }
    }
};

/// Factory initialisation for the IProductWRTBase_SumFac_Prism operator
OperatorKey IProductWRTBase_SumFac_Prism::m_type =
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePrism, eIProductWRTBase, eSumFac, false),
        IProductWRTBase_SumFac_Prism::create, "IProductWRTBase_SumFac_Prism");

/**
 * @brief Inner Product operator using sum-factorisation (Pyr)
 */
class IProductWRTBase_SumFac_Pyr final : virtual public Operator,
                                         IProductWRTBase_Helper
{
public:
    OPERATOR_CREATE(IProductWRTBase_SumFac_Pyr)

    ~IProductWRTBase_SumFac_Pyr() final = default;

    void operator()(const Array<OneD, const NekDouble> &input,
                    Array<OneD, NekDouble> &output,
                    [[maybe_unused]] Array<OneD, NekDouble> &output1,
                    [[maybe_unused]] Array<OneD, NekDouble> &output2,
                    Array<OneD, NekDouble> &wsp) final
    {
        ASSERTL1(wsp.size() == m_wspSize, "Incorrect workspace size");

        PyrIProduct(m_sortTopVertex, m_numElmt, m_nquad0, m_nquad1, m_nquad2,
                    m_nmodes0, m_nmodes1, m_nmodes2, m_base0, m_base1, m_base2,
                    m_jacWStdW, input, output, wsp);
    }

    void operator()([[maybe_unused]] int dir,
                    [[maybe_unused]] const Array<OneD, const NekDouble> &input,
                    [[maybe_unused]] Array<OneD, NekDouble> &output,
                    [[maybe_unused]] Array<OneD, NekDouble> &wsp) final
    {
        NEKERROR(ErrorUtil::efatal, "Not valid for this operator.");
    }

protected:
    const int m_nquad0;
    const int m_nquad1;
    const int m_nquad2;
    const int m_nmodes0;
    const int m_nmodes1;
    const int m_nmodes2;
    Array<OneD, const NekDouble> m_jacWStdW;
    Array<OneD, const NekDouble> m_base0;
    Array<OneD, const NekDouble> m_base1;
    Array<OneD, const NekDouble> m_base2;
    bool m_sortTopVertex;

private:
    IProductWRTBase_SumFac_Pyr(
        vector<StdRegions::StdExpansionSharedPtr> pCollExp,
        CoalescedGeomDataSharedPtr pGeomData, StdRegions::FactorMap factors)
        : Operator(pCollExp, pGeomData, factors), IProductWRTBase_Helper(),
          m_nquad0(m_stdExp->GetNumPoints(0)),
          m_nquad1(m_stdExp->GetNumPoints(1)),
          m_nquad2(m_stdExp->GetNumPoints(2)),
          m_nmodes0(m_stdExp->GetBasisNumModes(0)),
          m_nmodes1(m_stdExp->GetBasisNumModes(1)),
          m_nmodes2(m_stdExp->GetBasisNumModes(2)),
          m_base0(m_stdExp->GetBasis(0)->GetBdata()),
          m_base1(m_stdExp->GetBasis(1)->GetBdata()),
          m_base2(m_stdExp->GetBasis(2)->GetBdata())

    {
        m_jacWStdW = pGeomData->GetJacWithStdWeights(pCollExp);

        m_wspSize = m_numElmt * m_nquad2 *
                        (max(m_nquad0 * m_nquad1, m_nmodes0 * m_nmodes1)) +
                    m_nquad1 * m_nquad2 * m_numElmt * m_nmodes0;

        if (m_stdExp->GetBasis(0)->GetBasisType() == LibUtilities::eModified_A)
        {
            m_sortTopVertex = true;
        }
        else
        {
            m_sortTopVertex = false;
        }
    }
};

/// Factory initialisation for the IProductWRTBase_SumFac_Pyr operator
OperatorKey IProductWRTBase_SumFac_Pyr::m_type =
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePyramid, eIProductWRTBase, eSumFac, false),
        IProductWRTBase_SumFac_Pyr::create, "IProductWRTBase_SumFac_Pyr");

} // namespace Nektar::Collections
