///////////////////////////////////////////////////////////////////////////////
//
// File: NekLinSysIterCGLoc.h
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
// Description: NekLinSysIterCGLoc header
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_LINSYS_ITERAT_CG_H
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_LINSYS_ITERAT_CG_H

#include <LibUtilities/LinearAlgebra/NekLinSysIter.h>

namespace Nektar::LibUtilities
{
/// A global linear system.
class NekLinSysIterCGLoc;

typedef std::shared_ptr<NekLinSysIterCGLoc> NekLinSysIterCGLocSharedPtr;

class NekLinSysIterCGLoc : public NekLinSysIter
{
public:
    /// Support creation through MemoryManager.
    friend class MemoryManager<NekLinSysIterCGLoc>;

    LIB_UTILITIES_EXPORT static NekLinSysIterSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const LibUtilities::CommSharedPtr &vRowComm, const int nDimen,
        const NekSysKey &pKey)
    {
        NekLinSysIterCGLocSharedPtr p =
            MemoryManager<NekLinSysIterCGLoc>::AllocateSharedPtr(
                pSession, vRowComm, nDimen, pKey);
        p->InitObject();
        return p;
    }

    static std::string className;

    /// Constructor for full direct matrix solve.
    LIB_UTILITIES_EXPORT NekLinSysIterCGLoc(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const LibUtilities::CommSharedPtr &vRowComm, const int nDimen,
        const NekSysKey &pKey);
    LIB_UTILITIES_EXPORT ~NekLinSysIterCGLoc() override = default;

protected:
    void v_InitObject() override;

    int v_SolveSystem(const int nGlobal,
                      const Array<OneD, const NekDouble> &pInput,
                      Array<OneD, NekDouble> &pOutput, const int nDir) override;

private:
    /// Actual iterative solve
    void DoConjugateGradient(const int pNumRows,
                             const Array<OneD, const NekDouble> &pInput,
                             Array<OneD, NekDouble> &pOutput);
};
} // namespace Nektar::LibUtilities

#endif
