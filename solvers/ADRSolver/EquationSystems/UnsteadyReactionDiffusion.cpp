///////////////////////////////////////////////////////////////////////////////
//
// File: UnsteadyReactionDiffusion.cpp
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
// Description: Unsteady reaction-diffusion solve routines
//
///////////////////////////////////////////////////////////////////////////////

#include <ADRSolver/EquationSystems/UnsteadyReactionDiffusion.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <iomanip>
#include <iostream>

using namespace std;

namespace Nektar
{
string UnsteadyReactionDiffusion::className =
    GetEquationSystemFactory().RegisterCreatorFunction(
        "UnsteadyReactionDiffusion", UnsteadyReactionDiffusion::create);

UnsteadyReactionDiffusion::UnsteadyReactionDiffusion(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : UnsteadyDiffusion(pSession, pGraph)
{
}

/**
 * @brief Initialisation object for the unsteady reaction-diffusion problem.
 */
void UnsteadyReactionDiffusion::v_InitObject(bool DeclareFields)
{
    UnsteadySystem::v_InitObject(DeclareFields);

    ASSERTL0(m_intScheme->GetIntegrationSchemeType() == LibUtilities::eIMEX,
             "Reaction-diffusion requires an implicit-explicit timestepping"
             " (e.g. IMEXOrder2)");

    // Load diffusion parameter
    m_session->LoadParameter("epsilon", m_epsilon, 0.0);

    m_session->MatchSolverInfo("SpectralVanishingViscosity", "True",
                               m_useSpecVanVisc, false);

    if (m_useSpecVanVisc)
    {
        m_session->LoadParameter("SVVCutoffRatio", m_sVVCutoffRatio, 0.75);
        m_session->LoadParameter("SVVDiffCoeff", m_sVVDiffCoeff, 0.1);
    }

    int npoints = m_fields[0]->GetNpoints();

    if (m_session->DefinesParameter("d00"))
    {
        m_d00 = m_session->GetParameter("d00");
        m_varcoeff[StdRegions::eVarCoeffD00] =
            Array<OneD, NekDouble>(npoints, m_session->GetParameter("d00"));
    }
    if (m_session->DefinesParameter("d11"))
    {
        m_d11 = m_session->GetParameter("d11");
        m_varcoeff[StdRegions::eVarCoeffD11] =
            Array<OneD, NekDouble>(npoints, m_session->GetParameter("d11"));
    }
    if (m_session->DefinesParameter("d22"))
    {
        m_d22 = m_session->GetParameter("d22");
        m_varcoeff[StdRegions::eVarCoeffD22] =
            Array<OneD, NekDouble>(npoints, m_session->GetParameter("d22"));
    }

    // Forcing terms
    m_forcing = SolverUtils::Forcing::Load(m_session, shared_from_this(),
                                           m_fields, m_fields.size());

    m_ode.DefineOdeRhs(&UnsteadyReactionDiffusion::DoOdeRhs, this);
    m_ode.DefineProjection(&UnsteadyReactionDiffusion::DoOdeProjection, this);
    m_ode.DefineImplicitSolve(&UnsteadyReactionDiffusion::DoImplicitSolve,
                              this);
}

/**
 * @brief Unsteady diffusion problem destructor.
 */
UnsteadyReactionDiffusion::~UnsteadyReactionDiffusion()
{
}

/**
 * @brief Compute the right-hand side for the unsteady reaction diffusion
 * problem.
 *
 * @param inarray    Given fields.
 * @param outarray   Calculated solution.
 * @param time       Time.
 */
void UnsteadyReactionDiffusion::DoOdeRhs(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    // RHS should be set to zero.
    for (int i = 0; i < outarray.size(); ++i)
    {
        Vmath::Zero(outarray[i].size(), &outarray[i][0], 1);
    }

    // Add forcing terms for reaction.
    for (auto &x : m_forcing)
    {
        // set up non-linear terms
        x->Apply(m_fields, inarray, outarray, time);
    }
}

/**
 * @brief Compute the projection for the unsteady diffusion problem.
 *
 * @param inarray    Given fields.
 * @param outarray   Calculated solution.
 * @param time       Time.
 */
void UnsteadyReactionDiffusion::DoOdeProjection(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    UnsteadyDiffusion::DoOdeProjection(inarray, outarray, time);
}

/**
 * @brief Implicit solution of the unsteady diffusion problem.
 */
void UnsteadyReactionDiffusion::DoImplicitSolve(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
    const NekDouble lambda)
{
    UnsteadyDiffusion::DoImplicitSolve(inarray, outarray, time, lambda);
}

} // namespace Nektar
