///////////////////////////////////////////////////////////////////////////////
//
// File: PreconditionerBlock.h
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
// Description: Block Preconditioner header
//
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_LIB_MULTIREGIONS_PRECONDITIONERBLOCK_H
#define NEKTAR_LIB_MULTIREGIONS_PRECONDITIONERBLOCK_H
#include <LocalRegions/PrismExp.h>
#include <LocalRegions/TetExp.h>
#include <MultiRegions/AssemblyMap/AssemblyMapCG.h>
#include <MultiRegions/GlobalLinSys.h>
#include <MultiRegions/MultiRegionsDeclspec.h>
#include <MultiRegions/Preconditioner.h>

namespace Nektar::MultiRegions
{
class PreconditionerBlock;
typedef std::shared_ptr<PreconditionerBlock> PreconditionerBlockSharedPtr;

class PreconditionerBlock : public Preconditioner
{
public:
    /// Creates an instance of this class
    static PreconditionerSharedPtr create(
        const std::shared_ptr<GlobalLinSys> &plinsys,
        const std::shared_ptr<AssemblyMap> &pLocToGloMap)
    {
        PreconditionerSharedPtr p =
            MemoryManager<PreconditionerBlock>::AllocateSharedPtr(plinsys,
                                                                  pLocToGloMap);
        p->InitObject();
        return p;
    }

    /// Name of class
    static std::string className;

    MULTI_REGIONS_EXPORT PreconditionerBlock(
        const std::shared_ptr<GlobalLinSys> &plinsys,
        const AssemblyMapSharedPtr &pLocToGloMap);

    MULTI_REGIONS_EXPORT
    ~PreconditionerBlock() override
    {
        Gs::Free(m_gs_mapping);
    }

protected:
    bool m_isFull;
    DNekBlkMatSharedPtr m_blkMat;

    // Hold gs mapping as member to hold and delete only at descrution
    Gs::gs_data *m_gs_mapping;

    void v_InitObject() override;
    void v_DoPreconditioner(const Array<OneD, NekDouble> &pInput,
                            Array<OneD, NekDouble> &pOutput,
                            const bool &isLocal = false) override;

    void v_BuildPreconditioner() override;

private:
    void BlockPreconditionerCG(void);
    void BlockPreconditionerHDG(void);
};
} // namespace Nektar::MultiRegions

#endif
