///////////////////////////////////////////////////////////////////////////////
//
// File: IProduct.h
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
// Description:
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBRARY_MF_IPRODUCT_H
#define NEKTAR_LIBRARY_MF_IPRODUCT_H

#include <LibUtilities/BasicUtils/ShapeType.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Foundations/Basis.h>

#include "Operator.hpp"

#include "IProductKernels.hpp"

namespace Nektar::MatrixFree
{

// As each operator has seven shapes over three dimension to get to the
// "work" each operator uses a series of preprocessor directives based
// on the dimension and shape type so to limit the code
// inclusion. This keeps the library size as minimal as possible while
// removing runtime conditionals.

// The preprocessor directives, SHAPE_DIMENSION_?D and SHAPE_TYPE_*
// are constructed by CMake in the CMakeLists.txt file which uses an
// implementation file.  See the CMakeLists.txt files for more
// details.
template <LibUtilities::ShapeType SHAPE_TYPE, bool DEFORMED = false>
struct IProductTemplate : public IProduct, public Helper<SHAPE_TYPE, DEFORMED>
{
    IProductTemplate(std::vector<LibUtilities::BasisSharedPtr> basis, int nElmt)
        : IProduct(basis, nElmt), Helper<SHAPE_TYPE, DEFORMED>(basis, nElmt)
    {
        constexpr auto DIM = LibUtilities::ShapeTypeDimMap[SHAPE_TYPE];

        if (DIM == 1)
        {
            m_nmTot = LibUtilities::GetNumberOfCoefficients(SHAPE_TYPE,
                                                            this->m_nm[0]);
        }
        else if (DIM == 2)
        {
            m_nmTot = LibUtilities::GetNumberOfCoefficients(
                SHAPE_TYPE, this->m_nm[0], this->m_nm[1]);
        }
        else if (DIM == 3)
        {
            m_nmTot = LibUtilities::GetNumberOfCoefficients(
                SHAPE_TYPE, this->m_nm[0], this->m_nm[1], this->m_nm[2]);
        }
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis, int nElmt)
    {
        return std::make_shared<IProductTemplate<SHAPE_TYPE, DEFORMED>>(basis,
                                                                        nElmt);
    }

    void operator()(const Array<OneD, const NekDouble> &input,
                    Array<OneD, NekDouble> &output) final
    {
        const int nm0 = this->m_nm[0];
        const int nq0 = this->m_nq[0];
#if defined(SHAPE_DIMENSION_2D)
        const int nm1 = this->m_nm[1];
        const int nq1 = this->m_nq[1];
#elif defined(SHAPE_DIMENSION_3D)
        const int nm1 = this->m_nm[1];
        const int nm2 = this->m_nm[2];
        const int nq1 = this->m_nq[1];
        const int nq2 = this->m_nq[2];
#endif
#include "SwitchNodesPoints.h"
    }

    // There must be separate 1D, 2D, and 3D operators because the
    // helper base class is based on the shape type dim. Thus indices
    // for the basis must be respected.

    // Further there are duplicate 1D, 2D, and 3D operators so to
    // allow for compile time array sizes to be part of the template
    // and thus gain loop unrolling. The only difference is the
    // location of the size value, in the template or the function:
    // foo<int bar>() {} vs foo() { int bar = ... }

#if defined(SHAPE_DIMENSION_1D)

    // Non-size based operator.
    void operator1D(const Array<OneD, const NekDouble> &input,
                    Array<OneD, NekDouble> &output)
    {
        const auto nm0 = this->m_nm[0];
        const auto nq0 = this->m_nq[0];

        const auto nqTot    = nq0;
        const auto nqBlocks = nqTot * vec_t::width;
        const auto nmBlocks = m_nmTot * vec_t::width;

        auto *inptr  = &input[0];
        auto *outptr = &output[0];

        std::vector<vec_t, allocator<vec_t>> tmpIn(nqTot), tmpOut(m_nmTot);

        vec_t *jac_ptr = &((*this->m_jac)[0]);

        // temporary aligned storage for local fields
        NekDouble *locField = static_cast<NekDouble *>(::operator new[](
            nqBlocks * sizeof(NekDouble), std::align_val_t(vec_t::alignment)));

        for (int e = 0; e < this->m_nBlocks - 1; ++e)
        {
            // Load data to aligned storage and interleave it
            // std::copy(inptr, inptr + nqBlocks, locField);
            // load_interleave(locField, nqTot, tmpIn);
            load_unalign_interleave(inptr, nqTot, tmpIn);

            IProduct1DKernel<SHAPE_TYPE, false, false, DEFORMED>(
                nm0, nq0, tmpIn, this->m_bdata[0], this->m_w[0], jac_ptr,
                tmpOut);

            // de-interleave and store data
            // deinterleave_store(tmpOut, m_nmTot, locField);
            // std::copy(locField, locField + nmBlocks, outptr);
            deinterleave_unalign_store(tmpOut, m_nmTot, outptr);

            inptr += nqBlocks;
            outptr += nmBlocks;
            if constexpr (DEFORMED)
            {
                jac_ptr += nqTot;
            }
            else
            {
                ++jac_ptr;
            }
        }
        // last block
        {
            int acturalSize = nqBlocks - this->m_nPads * nqTot;
            std::copy(inptr, inptr + acturalSize, locField);
            load_interleave(locField, nqTot, tmpIn);

            IProduct1DKernel<SHAPE_TYPE, false, false, DEFORMED>(
                nm0, nq0, tmpIn, this->m_bdata[0], this->m_w[0], jac_ptr,
                tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, locField);
            acturalSize = nmBlocks - this->m_nPads * m_nmTot;
            std::copy(locField, locField + acturalSize, outptr);
        }
        // free aligned memory
        ::operator delete[](locField, std::align_val_t(vec_t::alignment));
    }

    // Size based template version.
    template <int nm0, int nq0>
    void operator1D(const Array<OneD, const NekDouble> &input,
                    Array<OneD, NekDouble> &output)
    {
        constexpr auto nqTot    = nq0;
        constexpr auto nqBlocks = nqTot * vec_t::width;
        const auto nmBlocks     = m_nmTot * vec_t::width;

        auto *inptr  = &input[0];
        auto *outptr = &output[0];

        std::vector<vec_t, allocator<vec_t>> tmpIn(nqTot), tmpOut(m_nmTot);

        vec_t *jac_ptr = &((*this->m_jac)[0]);

        // temporary aligned storage for local fields
        alignas(vec_t::alignment) NekDouble locField[nqBlocks];

        for (int e = 0; e < this->m_nBlocks - 1; ++e)
        {
            // Load data to aligned storage and interleave it
            // std::copy(inptr, inptr + nqBlocks, locField);
            // load_interleave(locField, nqTot, tmpIn);
            load_unalign_interleave(inptr, nqTot, tmpIn);

            IProduct1DKernel<SHAPE_TYPE, false, false, DEFORMED>(
                nm0, nq0, tmpIn, this->m_bdata[0], this->m_w[0], jac_ptr,
                tmpOut);

            // de-interleave and store data
            // deinterleave_store(tmpOut, m_nmTot, locField);
            // std::copy(locField, locField + nmBlocks, outptr);
            deinterleave_unalign_store(tmpOut, m_nmTot, outptr);

            inptr += nqBlocks;
            outptr += nmBlocks;
            if constexpr (DEFORMED)
            {
                jac_ptr += nqTot;
            }
            else
            {
                ++jac_ptr;
            }
        }
        // last block
        {
            int acturalSize = nqBlocks - this->m_nPads * nqTot;
            std::copy(inptr, inptr + acturalSize, locField);
            load_interleave(locField, nqTot, tmpIn);

            IProduct1DKernel<SHAPE_TYPE, false, false, DEFORMED>(
                nm0, nq0, tmpIn, this->m_bdata[0], this->m_w[0], jac_ptr,
                tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, locField);
            acturalSize = nmBlocks - this->m_nPads * m_nmTot;
            std::copy(locField, locField + acturalSize, outptr);
        }
    }

#elif defined(SHAPE_DIMENSION_2D)

    // Non-size based operator.
    void operator2D(const Array<OneD, const NekDouble> &input,
                    Array<OneD, NekDouble> &output)
    {
        const auto nm0 = this->m_nm[0];
        const auto nm1 = this->m_nm[1];

        const auto nq0 = this->m_nq[0];
        const auto nq1 = this->m_nq[1];

        const auto nqTot    = nq0 * nq1;
        const auto nqBlocks = nqTot * vec_t::width;
        const auto nmBlocks = m_nmTot * vec_t::width;

        auto *inptr  = &input[0];
        auto *outptr = &output[0];

        const bool correct =
            (m_basis[0]->GetBasisType() == LibUtilities::eModified_A);

        // Workspace for kernels - also checks preconditions
        size_t wsp0Size = 0;
        IProduct2DWorkspace<SHAPE_TYPE>(nm0, nm1, nq0, nq1, wsp0Size);

        std::vector<vec_t, allocator<vec_t>> wsp0(wsp0Size), tmpIn(nqTot),
            tmpOut(m_nmTot);

        vec_t *jac_ptr = &((*this->m_jac)[0]);

        // temporary aligned storage for local fields
        NekDouble *locField = static_cast<NekDouble *>(::operator new[](
            nqBlocks * sizeof(NekDouble), std::align_val_t(vec_t::alignment)));

        for (int e = 0; e < this->m_nBlocks - 1; ++e)
        {
            // Load data to aligned storage and interleave it
            // std::copy(inptr, inptr + nqBlocks, locField);
            // load_interleave(locField, nqTot, tmpIn);
            load_unalign_interleave(inptr, nqTot, tmpIn);

            IProduct2DKernel<SHAPE_TYPE, false, false, DEFORMED>(
                nm0, nm1, nq0, nq1, correct, tmpIn, this->m_bdata[0],
                this->m_bdata[1], this->m_w[0], this->m_w[1], jac_ptr, wsp0,
                tmpOut);

            // de-interleave and store data
            // deinterleave_store(tmpOut, m_nmTot, locField);
            // std::copy(locField, locField + nmBlocks, outptr);
            deinterleave_unalign_store(tmpOut, m_nmTot, outptr);

            inptr += nqBlocks;
            outptr += nmBlocks;
            if constexpr (DEFORMED)
            {
                jac_ptr += nqTot;
            }
            else
            {
                ++jac_ptr;
            }
        }
        // last block
        {
            int acturalSize = nqBlocks - this->m_nPads * nqTot;
            std::copy(inptr, inptr + acturalSize, locField);
            load_interleave(locField, nqTot, tmpIn);

            IProduct2DKernel<SHAPE_TYPE, false, false, DEFORMED>(
                nm0, nm1, nq0, nq1, correct, tmpIn, this->m_bdata[0],
                this->m_bdata[1], this->m_w[0], this->m_w[1], jac_ptr, wsp0,
                tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, locField);
            acturalSize = nmBlocks - this->m_nPads * m_nmTot;
            std::copy(locField, locField + acturalSize, outptr);
        }
        // free aligned memory
        ::operator delete[](locField, std::align_val_t(vec_t::alignment));
    }

    // Size based template version.
    template <int nm0, int nm1, int nq0, int nq1>
    void operator2D(const Array<OneD, const NekDouble> &input,
                    Array<OneD, NekDouble> &output)
    {
        constexpr auto nqTot    = nq0 * nq1;
        constexpr auto nqBlocks = nqTot * vec_t::width;
        const auto nmBlocks     = m_nmTot * vec_t::width;

        auto *inptr  = &input[0];
        auto *outptr = &output[0];

        const bool correct =
            (m_basis[0]->GetBasisType() == LibUtilities::eModified_A);

        // Workspace for kernels - also checks preconditions
        size_t wsp0Size = 0;
        IProduct2DWorkspace<SHAPE_TYPE>(nm0, nm1, nq0, nq1, wsp0Size);

        std::vector<vec_t, allocator<vec_t>> wsp0(wsp0Size), tmpIn(nqTot),
            tmpOut(m_nmTot);

        vec_t *jac_ptr = &((*this->m_jac)[0]);

        // temporary aligned storage for local fields
        alignas(vec_t::alignment) NekDouble locField[nqBlocks];

        for (int e = 0; e < this->m_nBlocks - 1; ++e)
        {
            // Load data to aligned storage and interleave it
            // std::copy(inptr, inptr + nqBlocks, locField);
            // load_interleave(locField, nqTot, tmpIn);
            load_unalign_interleave(inptr, nqTot, tmpIn);

            IProduct2DKernel<SHAPE_TYPE, false, false, DEFORMED>(
                nm0, nm1, nq0, nq1, correct, tmpIn, this->m_bdata[0],
                this->m_bdata[1], this->m_w[0], this->m_w[1], jac_ptr, wsp0,
                tmpOut);

            // de-interleave and store data
            // deinterleave_store(tmpOut, m_nmTot, locField);
            // std::copy(locField, locField + nmBlocks, outptr);
            deinterleave_unalign_store(tmpOut, m_nmTot, outptr);

            inptr += nqBlocks;
            outptr += nmBlocks;
            if constexpr (DEFORMED)
            {
                jac_ptr += nqTot;
            }
            else
            {
                ++jac_ptr;
            }
        }
        // last block
        {
            int acturalSize = nqBlocks - this->m_nPads * nqTot;
            std::copy(inptr, inptr + acturalSize, locField);
            load_interleave(locField, nqTot, tmpIn);

            IProduct2DKernel<SHAPE_TYPE, false, false, DEFORMED>(
                nm0, nm1, nq0, nq1, correct, tmpIn, this->m_bdata[0],
                this->m_bdata[1], this->m_w[0], this->m_w[1], jac_ptr, wsp0,
                tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, locField);
            acturalSize = nmBlocks - this->m_nPads * m_nmTot;
            std::copy(locField, locField + acturalSize, outptr);
        }
    }

#elif defined(SHAPE_DIMENSION_3D)

    // Non-size based operator.
    void operator3D(const Array<OneD, const NekDouble> &input,
                    Array<OneD, NekDouble> &output)
    {
        const auto nm0 = this->m_nm[0];
        const auto nm1 = this->m_nm[1];
        const auto nm2 = this->m_nm[2];

        const auto nq0 = this->m_nq[0];
        const auto nq1 = this->m_nq[1];
        const auto nq2 = this->m_nq[2];

        const auto nqTot    = nq0 * nq1 * nq2;
        const auto nqBlocks = nqTot * vec_t::width;
        const auto nmBlocks = m_nmTot * vec_t::width;

        auto *inptr  = &input[0];
        auto *outptr = &output[0];

        const bool correct =
            (m_basis[0]->GetBasisType() == LibUtilities::eModified_A);

        // Workspace for kernels - also checks preconditions
        size_t wsp0Size = 0, wsp1Size = 0, wsp2Size = 0;
        IProduct3DWorkspace<SHAPE_TYPE>(nm0, nm1, nm2, nq0, nq1, nq2, wsp0Size,
                                        wsp1Size, wsp2Size);

        std::vector<vec_t, allocator<vec_t>> wsp0(wsp0Size), wsp1(wsp1Size),
            wsp2(wsp2Size), tmpIn(nqTot), tmpOut(m_nmTot);

        vec_t *jac_ptr = &((*this->m_jac)[0]);

        // temporary aligned storage for local fields
        NekDouble *locField = static_cast<NekDouble *>(::operator new[](
            nqBlocks * sizeof(NekDouble), std::align_val_t(vec_t::alignment)));

        for (int e = 0; e < this->m_nBlocks - 1; ++e)
        {
            // Load data to aligned storage and interleave it
            // std::copy(inptr, inptr + nqBlocks, locField);
            // load_interleave(locField, nqTot, tmpIn);
            load_unalign_interleave(inptr, nqTot, tmpIn);

            IProduct3DKernel<SHAPE_TYPE, false, false, DEFORMED>(
                nm0, nm1, nm2, nq0, nq1, nq2, correct, tmpIn, this->m_bdata[0],
                this->m_bdata[1], this->m_bdata[2], this->m_w[0], this->m_w[1],
                this->m_w[2], jac_ptr, wsp0, wsp1, wsp2, tmpOut);

            // de-interleave and store data
            // deinterleave_store(tmpOut, m_nmTot, locField);
            // std::copy(locField, locField + nmBlocks, outptr);
            deinterleave_unalign_store(tmpOut, m_nmTot, outptr);

            inptr += nqBlocks;
            outptr += nmBlocks;
            if constexpr (DEFORMED)
            {
                jac_ptr += nqTot;
            }
            else
            {
                ++jac_ptr;
            }
        }
        // last block
        {
            int acturalSize = nqBlocks - this->m_nPads * nqTot;
            std::copy(inptr, inptr + acturalSize, locField);
            load_interleave(locField, nqTot, tmpIn);

            IProduct3DKernel<SHAPE_TYPE, false, false, DEFORMED>(
                nm0, nm1, nm2, nq0, nq1, nq2, correct, tmpIn, this->m_bdata[0],
                this->m_bdata[1], this->m_bdata[2], this->m_w[0], this->m_w[1],
                this->m_w[2], jac_ptr, wsp0, wsp1, wsp2, tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, locField);
            acturalSize = nmBlocks - this->m_nPads * m_nmTot;
            std::copy(locField, locField + acturalSize, outptr);
        }
        // free aligned memory
        ::operator delete[](locField, std::align_val_t(vec_t::alignment));
    }

    // Size based template version.
    template <int nm0, int nm1, int nm2, int nq0, int nq1, int nq2>
    void operator3D(const Array<OneD, const NekDouble> &input,
                    Array<OneD, NekDouble> &output)
    {
        constexpr auto nqTot    = nq0 * nq1 * nq2;
        constexpr auto nqBlocks = nqTot * vec_t::width;
        const auto nmBlocks     = m_nmTot * vec_t::width;

        auto *inptr  = &input[0];
        auto *outptr = &output[0];

        const bool correct =
            (m_basis[0]->GetBasisType() == LibUtilities::eModified_A);

        // Workspace for kernels - also checks preconditions
        size_t wsp0Size = 0, wsp1Size = 0, wsp2Size = 0;
        IProduct3DWorkspace<SHAPE_TYPE>(nm0, nm1, nm2, nq0, nq1, nq2, wsp0Size,
                                        wsp1Size, wsp2Size);

        std::vector<vec_t, allocator<vec_t>> wsp0(wsp0Size), wsp1(wsp1Size),
            wsp2(wsp2Size), tmpIn(nqTot), tmpOut(m_nmTot);

        vec_t *jac_ptr = &((*this->m_jac)[0]);

        // temporary aligned storage for local fields
        alignas(vec_t::alignment) NekDouble locField[nqBlocks];

        for (int e = 0; e < this->m_nBlocks - 1; ++e)
        {
            // Load data to aligned storage and interleave it
            // std::copy(inptr, inptr + nqBlocks, locField);
            // load_interleave(locField, nqTot, tmpIn);
            load_unalign_interleave(inptr, nqTot, tmpIn);

            IProduct3DKernel<SHAPE_TYPE, false, false, DEFORMED>(
                nm0, nm1, nm2, nq0, nq1, nq2, correct, tmpIn, this->m_bdata[0],
                this->m_bdata[1], this->m_bdata[2], this->m_w[0], this->m_w[1],
                this->m_w[2], jac_ptr, wsp0, wsp1, wsp2, tmpOut);

            // de-interleave and store data
            // deinterleave_store(tmpOut, m_nmTot, locField);
            // std::copy(locField, locField + nmBlocks, outptr);
            deinterleave_unalign_store(tmpOut, m_nmTot, outptr);

            inptr += nqBlocks;
            outptr += nmBlocks;
            if constexpr (DEFORMED)
            {
                jac_ptr += nqTot;
            }
            else
            {
                ++jac_ptr;
            }
        }
        // last block
        {
            int acturalSize = nqBlocks - this->m_nPads * nqTot;
            std::copy(inptr, inptr + acturalSize, locField);
            load_interleave(locField, nqTot, tmpIn);

            IProduct3DKernel<SHAPE_TYPE, false, false, DEFORMED>(
                nm0, nm1, nm2, nq0, nq1, nq2, correct, tmpIn, this->m_bdata[0],
                this->m_bdata[1], this->m_bdata[2], this->m_w[0], this->m_w[1],
                this->m_w[2], jac_ptr, wsp0, wsp1, wsp2, tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, locField);
            acturalSize = nmBlocks - this->m_nPads * m_nmTot;
            std::copy(locField, locField + acturalSize, outptr);
        }
    }

#endif // SHAPE_DIMENSION

private:
    int m_nmTot;
};

} // namespace Nektar::MatrixFree

#endif
