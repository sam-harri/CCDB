///////////////////////////////////////////////////////////////////////////////
//
// File: EulerImplicitCFE.cpp
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
// Description: Euler Implicit equations in consƒervative variables without
// artificial diffusion
//
///////////////////////////////////////////////////////////////////////////////

#include <CompressibleFlowSolver/EquationSystems/EulerImplicitCFE.h>

using namespace std;

namespace Nektar
{
string EulerImplicitCFE::className =
    SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
        "EulerImplicitCFE", EulerImplicitCFE::create,
        "Euler Implicit equations in conservative variables.");

/**
 *
 */
EulerImplicitCFE::EulerImplicitCFE(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : UnsteadySystem(pSession, pGraph),
      CompressibleFlowSystem(pSession, pGraph), CFSImplicit(pSession, pGraph)
{
}

/**
 *
 */
void EulerImplicitCFE::v_InitObject(bool DeclareFields)
{
    CFSImplicit::v_InitObject(DeclareFields);
    m_viscousJacFlag = false;
}

/**
 *
 */
bool EulerImplicitCFE::v_SupportsShockCaptType(const std::string type) const
{
    return type == "Off";
}
} // namespace Nektar
