///////////////////////////////////////////////////////////////////////////////
//
// File: PreconditionerDiagonal.h
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
// Description: PreconditionerDiagonal header
//
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_LIB_MULTIREGIONS_PRECONDITIONERDIAGONAL_H
#define NEKTAR_LIB_MULTIREGIONS_PRECONDITIONERDIAGONAL_H
#include <MultiRegions/AssemblyMap/AssemblyMapCG.h>
#include <MultiRegions/GlobalLinSys.h>
#include <MultiRegions/MultiRegionsDeclspec.h>
#include <MultiRegions/Preconditioner.h>

namespace Nektar::MultiRegions
{
class PreconditionerDiagonal;
typedef std::shared_ptr<PreconditionerDiagonal> PreconditionerDiagonalSharedPtr;

class PreconditionerDiagonal : public Preconditioner
{
public:
    /// Creates an instance of this class
    static PreconditionerSharedPtr create(
        const std::shared_ptr<GlobalLinSys> &plinsys,
        const std::shared_ptr<AssemblyMap> &pLocToGloMap)
    {
        PreconditionerSharedPtr p =
            MemoryManager<PreconditionerDiagonal>::AllocateSharedPtr(
                plinsys, pLocToGloMap);
        p->InitObject();
        return p;
    }

    /// Name of class
    static std::string className;
    // static std::string className1;

    MULTI_REGIONS_EXPORT PreconditionerDiagonal(
        const std::shared_ptr<GlobalLinSys> &plinsys,
        const AssemblyMapSharedPtr &pLocToGloMap);

    MULTI_REGIONS_EXPORT
    ~PreconditionerDiagonal() override
    {
    }

protected:
    Array<OneD, NekDouble> m_diagonals;

    void v_InitObject() override;

    void v_DoPreconditioner(const Array<OneD, NekDouble> &pInput,
                            Array<OneD, NekDouble> &pOutput,
                            const bool &IsLocal = false) override;

    void v_BuildPreconditioner() override;

private:
    void DiagonalPreconditionerSum(void);

    void StaticCondDiagonalPreconditionerSum(void);

    static std::string lookupIds[];
    static std::string def;
};

class PreconditionerNull;
typedef std::shared_ptr<PreconditionerNull> PreconditionerNullSharedPtr;

class PreconditionerNull : public Preconditioner
{
public:
    /// Creates an instance of this class
    static PreconditionerSharedPtr create(
        const std::shared_ptr<GlobalLinSys> &plinsys,
        const std::shared_ptr<AssemblyMap> &pLocToGloMap)
    {
        PreconditionerSharedPtr p =
            MemoryManager<PreconditionerNull>::AllocateSharedPtr(plinsys,
                                                                 pLocToGloMap);
        p->InitObject();
        return p;
    }

    /// Name of class
    static std::string className;

    MULTI_REGIONS_EXPORT PreconditionerNull(
        const std::shared_ptr<GlobalLinSys> &plinsys,
        const AssemblyMapSharedPtr &pLocToGloMap);

    MULTI_REGIONS_EXPORT
    ~PreconditionerNull() override
    {
    }

    void v_InitObject() override;

    void v_DoPreconditioner(const Array<OneD, NekDouble> &pInput,
                            Array<OneD, NekDouble> &pOutput,
                            const bool &isLocal = false) override;

    void v_BuildPreconditioner() override;

private:
    static std::string lookupIds[];
    static std::string def;
};

class PreconditionerJacobi;
typedef std::shared_ptr<PreconditionerJacobi> PreconditionerJacobiSharedPtr;

class PreconditionerJacobi : public PreconditionerDiagonal
{
public:
    /// Creates an instance of this class
    static PreconditionerSharedPtr create(
        const std::shared_ptr<GlobalLinSys> &plinsys,
        const std::shared_ptr<AssemblyMap> &pLocToGloMap)
    {
        PreconditionerSharedPtr p =
            MemoryManager<PreconditionerJacobi>::AllocateSharedPtr(
                plinsys, pLocToGloMap);
        p->InitObject();
        return p;
    }

    /// Name of class
    static std::string className;
    // static std::string className1;

    MULTI_REGIONS_EXPORT PreconditionerJacobi(
        const std::shared_ptr<GlobalLinSys> &plinsys,
        const AssemblyMapSharedPtr &pLocToGloMap);

    MULTI_REGIONS_EXPORT
    ~PreconditionerJacobi() override
    {
    }

protected:
    void v_InitObject() override;

    void v_BuildPreconditioner() override;

    void v_DoPreconditioner(const Array<OneD, NekDouble> &pInput,
                            Array<OneD, NekDouble> &pOutput,
                            const bool &IsLocal = false) override;

private:
    int m_niter;

    static std::string lookupIds[];
    static std::string def;
};

} // namespace Nektar::MultiRegions

#endif
