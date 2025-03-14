///////////////////////////////////////////////////////////////////////////////
//
// File: TimeDependentBC.cpp
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
// Description: Time dependent boundary condition
//
///////////////////////////////////////////////////////////////////////////////

#include "TimeDependentBC.h"

using namespace std;

namespace Nektar
{

std::string TimeDependentBC::className =
    GetCFSBndCondFactory().RegisterCreatorFunction(
        "TimeDependent", TimeDependentBC::create,
        "Time dependent boundary condition.");

TimeDependentBC::TimeDependentBC(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    const Array<OneD, Array<OneD, NekDouble>> &pTraceNormals,
    const Array<OneD, Array<OneD, NekDouble>> &pGridVelocity,
    const int pSpaceDim, const int bcRegion, const int cnt)
    : CFSBndCond(pSession, pFields, pTraceNormals, pGridVelocity, pSpaceDim,
                 bcRegion, cnt)
{
}

void TimeDependentBC::v_Apply(
    [[maybe_unused]] Array<OneD, Array<OneD, NekDouble>> &Fwd,
    Array<OneD, Array<OneD, NekDouble>> &physarray, const NekDouble &time)
{
    int nvariables = physarray.size();
    std::string varName;
    for (int i = 0; i < nvariables; ++i)
    {
        varName = m_session->GetVariable(i);
        m_fields[i]->EvaluateBoundaryConditions(time, varName);
    }
}

} // namespace Nektar
