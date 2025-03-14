///////////////////////////////////////////////////////////////////////////////
//
// File: ImplicitExtrapolate.h
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
// Description: Abstract base class for ImplicitExtrapolate.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_IMPLICITEXTRAPOLATE_H
#define NEKTAR_SOLVERS_IMPLICITEXTRAPOLATE_H

#include <IncNavierStokesSolver/EquationSystems/WeakPressureExtrapolate.h>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <MultiRegions/ExpList.h>

namespace Nektar
{
//--------
// Implicit Extrapolate
// --------

class ImplicitExtrapolate;

typedef std::shared_ptr<ImplicitExtrapolate> ImplicitExtrapolateSharedPtr;

class ImplicitExtrapolate : public WeakPressureExtrapolate
{
public:
    /// Creates an instance of this class
    static ExtrapolateSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        MultiRegions::ExpListSharedPtr &pPressure, const Array<OneD, int> &pVel,
        const SolverUtils::AdvectionSharedPtr &advObject)
    {
        ExtrapolateSharedPtr p =
            MemoryManager<ImplicitExtrapolate>::AllocateSharedPtr(
                pSession, pFields, pPressure, pVel, advObject);
        return p;
    }

    /// Name of class
    static std::string className;

    ImplicitExtrapolate(const LibUtilities::SessionReaderSharedPtr pSession,
                        Array<OneD, MultiRegions::ExpListSharedPtr> pFields,
                        MultiRegions::ExpListSharedPtr pPressure,
                        const Array<OneD, int> pVel,
                        const SolverUtils::AdvectionSharedPtr advObject);

    ~ImplicitExtrapolate() override;

protected:
    void v_EvaluatePressureBCs(
        const Array<OneD, const Array<OneD, NekDouble>> &fields,
        const Array<OneD, const Array<OneD, NekDouble>> &N,
        NekDouble kinvis) override;

    static std::string solverTypeLookupId;
};

} // namespace Nektar

#endif
