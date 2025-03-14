///////////////////////////////////////////////////////////////////////////////
//
// File: sse2.hpp
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
// Description: Vector type using sse2 extension.
// Note that this is not a full implementation: only the int type is
// implemented as support index type for avx2.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_LIBUTILITES_SIMDLIB_SSE2_H
#define NEKTAR_LIB_LIBUTILITES_SIMDLIB_SSE2_H

#if defined(__x86_64__)
#include <immintrin.h>
#if defined(__INTEL_COMPILER) && !defined(TINYSIMD_HAS_SVML)
#define TINYSIMD_HAS_SVML
#endif
#endif
#include "traits.hpp"
#include <cstdint>

namespace tinysimd::abi
{

template <typename scalarType> struct sse2
{
    using type = void;
};

} // namespace tinysimd::abi

#if defined(__SSE2__) && defined(NEKTAR_ENABLE_SIMD_SSE2)

namespace tinysimd
{

// forward declaration of concrete types
template <typename T> struct sse2Int4;

namespace abi
{

// mapping between abstract types and concrete types
template <> struct sse2<std::int32_t>
{
    using type = sse2Int4<std::int32_t>;
};
template <> struct sse2<std::uint32_t>
{
    using type = sse2Int4<std::uint32_t>;
};

} // namespace abi

// concrete types
template <typename T> struct sse2Int4
{
    static_assert(std::is_integral<T>::value && sizeof(T) == 4,
                  "4 bytes Integral required.");

    static constexpr unsigned int width     = 4;
    static constexpr unsigned int alignment = 16;

    using scalarType  = T;
    using vectorType  = __m128i;
    using scalarArray = scalarType[width];

    // storage
    vectorType _data;

    // ctors
    inline sse2Int4()                    = default;
    inline sse2Int4(const sse2Int4 &rhs) = default;
    inline sse2Int4(const vectorType &rhs) : _data(rhs)
    {
    }
    inline sse2Int4(const scalarType rhs)
    {
        _data = _mm_set1_epi32(rhs);
    }

    // store
    inline void store(scalarType *p) const
    {
        _mm_store_si128(reinterpret_cast<vectorType *>(p), _data);
    }

    template <class flag,
              typename std::enable_if<is_requiring_alignment<flag>::value &&
                                          !is_streaming<flag>::value,
                                      bool>::type = 0>
    inline void store(scalarType *p, flag) const
    {
        _mm_store_si128(reinterpret_cast<vectorType *>(p), _data);
    }

    template <class flag,
              typename std::enable_if<!is_requiring_alignment<flag>::value,
                                      bool>::type = 0>
    inline void store(scalarType *p, flag) const
    {
        _mm_storeu_si128(reinterpret_cast<vectorType *>(p), _data);
    }

    inline void load(const scalarType *p)
    {
        _data = _mm_load_si128(reinterpret_cast<const vectorType *>(p));
    }

    template <class flag,
              typename std::enable_if<is_requiring_alignment<flag>::value &&
                                          !is_streaming<flag>::value,
                                      bool>::type = 0>
    inline void load(const scalarType *p, flag)
    {
        _data = _mm_load_si128(reinterpret_cast<const vectorType *>(p));
    }

    template <class flag,
              typename std::enable_if<!is_requiring_alignment<flag>::value,
                                      bool>::type = 0>
    inline void load(const scalarType *p, flag)
    {
        _data = _mm_loadu_si128(reinterpret_cast<const vectorType *>(p));
    }

    // gather/scatter with sse2
    inline void gather(scalarType const *p, const sse2Int4<T> &indices)
    {
        _data = _mm_i32gather_epi32(p, indices._data, 8);
    }

    inline void scatter(scalarType *out, const sse2Int4<T> &indices) const
    {
        // no scatter intrinsics for sse2
        alignas(alignment) scalarArray tmp;
        _mm_store_epi32(tmp, _data);

        out[_mm_extract_epi32(indices._data, 0)] = tmp[0]; // SSE4.1
        out[_mm_extract_epi32(indices._data, 1)] = tmp[1];
    }

    inline void broadcast(const scalarType rhs)
    {
        _data = _mm_set1_epi32(rhs);
    }

    // subscript
    // subscript operators are convienient but expensive
    // should not be used in optimized kernels
    inline scalarType operator[](size_t i) const
    {
        alignas(alignment) scalarArray tmp;
        store(tmp, is_aligned);
        return tmp[i];
    }
};

} // namespace tinysimd

#endif // defined(__SSE2__) && defined(NEKTAR_ENABLE_SIMD_SSE2)

#endif
