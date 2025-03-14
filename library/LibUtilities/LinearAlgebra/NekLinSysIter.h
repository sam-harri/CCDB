///////////////////////////////////////////////////////////////////////////////
//
// File: NekLinSysIter.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
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
// Description: NekLinSysIter header
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_LINSYS_ITERAT_H
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_LINSYS_ITERAT_H

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/LinearAlgebra/NekSys.h>

namespace Nektar::LibUtilities
{

class NekLinSysIter;

typedef std::shared_ptr<NekLinSysIter> NekLinSysIterSharedPtr;

typedef LibUtilities::NekFactory<
    std::string, NekLinSysIter, const LibUtilities::SessionReaderSharedPtr &,
    const LibUtilities::CommSharedPtr &, const int, const NekSysKey &>
    NekLinSysIterFactory;
LIB_UTILITIES_EXPORT NekLinSysIterFactory &GetNekLinSysIterFactory();

class NekLinSysIter : public NekSys
{
public:
    /// Support creation through MemoryManager.
    friend class MemoryManager<NekLinSysIter>;

    LIB_UTILITIES_EXPORT static NekLinSysIterSharedPtr CreateInstance(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const LibUtilities::CommSharedPtr &vRowComm, const int nDimen,
        const NekSysKey &pKey)
    {
        NekLinSysIterSharedPtr p =
            MemoryManager<NekLinSysIter>::AllocateSharedPtr(pSession, vRowComm,
                                                            nDimen, pKey);
        return p;
    }

    LIB_UTILITIES_EXPORT NekLinSysIter(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const LibUtilities::CommSharedPtr &vRowComm, const int nDimen,
        const NekSysKey &pKey);
    LIB_UTILITIES_EXPORT ~NekLinSysIter() override = default;

    LIB_UTILITIES_EXPORT void SetUniversalUniqueMap(
        const Array<OneD, const int> &map);

    LIB_UTILITIES_EXPORT void SetNekLinSysTolerance(const NekDouble in)
    {
        m_NekLinSysTolerance = fmax(1.0E-16, in);
    }

    LIB_UTILITIES_EXPORT void SetNekLinSysMaxIterations(const unsigned int in)
    {
        m_NekLinSysMaxIterations = in;
    }

    LIB_UTILITIES_EXPORT int GetNekLinSysTolerance()
    {
        return m_NekLinSysTolerance;
    }

    LIB_UTILITIES_EXPORT bool IsLocal()
    {
        return m_isLocal;
    }

protected:
    /// Global to universal unique map
    Array<OneD, int> m_map;

    NekDouble m_NekLinSysTolerance;
    int m_NekLinSysMaxIterations;
    int m_totalIterations = 0;

    // Boolean to identify if iteration acts on local storage
    bool m_isLocal;

    void v_InitObject() override;

    void SetUniversalUniqueMap();

    void Set_Rhs_Magnitude(const Array<OneD, NekDouble> &pIn);

    void ConvergenceCheck(const Array<OneD, const NekDouble> &Residual);

private:
};

} // namespace Nektar::LibUtilities

#endif
