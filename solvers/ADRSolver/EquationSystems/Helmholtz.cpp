///////////////////////////////////////////////////////////////////////////////
//
// File: Helmholtz.cpp
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
// Description: Helmholtz solve routines
//
///////////////////////////////////////////////////////////////////////////////

#include <ADRSolver/EquationSystems/Helmholtz.h>

namespace Nektar
{
std::string Helmholtz::className =
    GetEquationSystemFactory().RegisterCreatorFunction("Helmholtz",
                                                       Helmholtz::create);

Helmholtz::Helmholtz(const LibUtilities::SessionReaderSharedPtr &pSession,
                     const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : Poisson(pSession, pGraph)
{
    if (pSession->DefinesParameter("Lambda"))
    {
        m_factors[StdRegions::eFactorLambda] =
            m_session->GetParameter("Lambda");
    }
}

void Helmholtz::v_InitObject(bool DeclareFields)
{
    Poisson::v_InitObject(DeclareFields);
}

void Helmholtz::v_GenerateSummary(SolverUtils::SummaryList &s)
{
    Poisson::v_GenerateSummary(s);
    SolverUtils::AddSummaryItem(s, "Lambda",
                                m_factors[StdRegions::eFactorLambda]);
}

Array<OneD, bool> Helmholtz::v_GetSystemSingularChecks()
{
    if (m_factors[StdRegions::eFactorLambda] == 0.0)
    {
        return Array<OneD, bool>(m_session->GetVariables().size(), true);
    }
    else
    {
        return Array<OneD, bool>(m_session->GetVariables().size(), false);
    }
}
} // namespace Nektar
