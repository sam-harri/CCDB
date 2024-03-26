///////////////////////////////////////////////////////////////////////////////
//
// File: ShallowWaterSystem.cpp
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
// Description: Generic timestepping for shallow water solvers
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <ShallowWaterSolver/EquationSystems/ShallowWaterSystem.h>

using namespace std;

namespace Nektar
{
/**
 * @class ShallowWaterSystem
 *
 * Provides the underlying timestepping framework for shallow water flow solvers
 * including the general timestepping routines. This class is not intended
 * to be directly instantiated, but rather is a base class on which to
 * define shallow water solvers, e.g. SWE, Boussinesq, linear and nonlinear
 * versions.
 *
 * For details on implementing unsteady solvers see
 * \ref sectionADRSolverModuleImplementation
 */

/**
 * Processes SolverInfo parameters from the session file and sets up
 * timestepping-specific code.
 * @param   pSession        Session object to read parameters from.
 */

string ShallowWaterSystem::className =
    SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
        "ShallowWaterSystem", ShallowWaterSystem::create,
        "Auxiliary functions for the shallow water system.");

ShallowWaterSystem::ShallowWaterSystem(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : UnsteadySystem(pSession, pGraph)
{
}

void ShallowWaterSystem::v_InitObject(bool DeclareFields)
{
    UnsteadySystem::v_InitObject(DeclareFields);

    // if discontinuous Galerkin determine numerical flux to use
    if (m_projectionType == MultiRegions::eDiscontinuous)
    {
        ASSERTL0(m_session->DefinesSolverInfo("UPWINDTYPE"),
                 "No UPWINDTYPE defined in session.");
    }

    // Set up locations of velocity vector.
    m_vecLocs    = Array<OneD, Array<OneD, NekDouble>>(1);
    m_vecLocs[0] = Array<OneD, NekDouble>(m_spacedim);
    for (int i = 0; i < m_spacedim; ++i)
    {
        m_vecLocs[0][i] = 1 + i;
    }

    // Load acceleration of gravity
    m_session->LoadParameter("Gravity", m_g, 9.81);

    // input/output in primitive variables
    m_primitive = true;

    EvaluateWaterDepth();

    m_constantDepth = true;
    NekDouble depth = m_depth[0];
    for (int i = 0; i < GetTotPoints(); ++i)
    {
        if (m_depth[i] != depth)
        {
            m_constantDepth = false;
            break;
        }
    }

    // Compute the bottom slopes
    int nq = GetTotPoints();
    if (m_constantDepth != true)
    {
        m_bottomSlope = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
        for (int i = 0; i < m_spacedim; ++i)
        {
            m_bottomSlope[i] = Array<OneD, NekDouble>(nq);
            m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[i], m_depth,
                                   m_bottomSlope[i]);
            Vmath::Neg(nq, m_bottomSlope[i], 1);
        }
    }

    EvaluateCoriolis();
}

void ShallowWaterSystem::v_GenerateSummary(SolverUtils::SummaryList &s)
{
    UnsteadySystem::v_GenerateSummary(s);
    if (m_constantDepth == true)
    {
        SolverUtils::AddSummaryItem(s, "Depth", "constant");
    }
    else
    {
        SolverUtils::AddSummaryItem(s, "Depth", "variable");
    }
}

void ShallowWaterSystem::DoOdeProjection(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    int i;
    int nvariables = inarray.size();

    switch (m_projectionType)
    {
        case MultiRegions::eDiscontinuous:
        {

            // Just copy over array
            if (inarray != outarray)
            {
                int npoints = GetNpoints();

                for (i = 0; i < nvariables; ++i)
                {
                    Vmath::Vcopy(npoints, inarray[i], 1, outarray[i], 1);
                }
            }

            SetBoundaryConditions(outarray, time);
            break;
        }
        case MultiRegions::eGalerkin:
        case MultiRegions::eMixed_CG_Discontinuous:
        {

            EquationSystem::SetBoundaryConditions(time);
            Array<OneD, NekDouble> coeffs(m_fields[0]->GetNcoeffs(), 0.0);

            for (i = 0; i < nvariables; ++i)
            {
                m_fields[i]->FwdTrans(inarray[i], coeffs);
                m_fields[i]->BwdTrans(coeffs, outarray[i]);
            }
            break;
        }
        default:
            ASSERTL0(false, "Unknown projection scheme");
            break;
    }
}

//----------------------------------------------------
void ShallowWaterSystem::SetBoundaryConditions(
    Array<OneD, Array<OneD, NekDouble>> &inarray, NekDouble time)
{
    std::string varName;
    int nvariables = 3;
    int cnt        = 0;
    int nTracePts  = GetTraceTotPoints();

    // Extract trace for boundaries. Needs to be done on all processors to avoid
    // deadlock.
    Array<OneD, Array<OneD, NekDouble>> Fwd(nvariables);
    for (int i = 0; i < nvariables; ++i)
    {
        Fwd[i] = Array<OneD, NekDouble>(nTracePts);
        m_fields[i]->ExtractTracePhys(inarray[i], Fwd[i]);
    }

    // Loop over Boundary Regions
    for (int n = 0; n < m_fields[0]->GetBndConditions().size(); ++n)
    {
        if (m_fields[0]->GetBndConditions()[n]->GetBoundaryConditionType() ==
            SpatialDomains::ePeriodic)
        {
            continue;
        }

        // Wall Boundary Condition
        if (boost::iequals(m_fields[0]->GetBndConditions()[n]->GetUserDefined(),
                           "Wall"))
        {
            WallBoundary2D(n, cnt, Fwd);
        }

        // Time Dependent Boundary Condition (specified in meshfile)
        if (m_fields[0]->GetBndConditions()[n]->IsTimeDependent())
        {
            for (int i = 0; i < nvariables; ++i)
            {
                varName = m_session->GetVariable(i);
                m_fields[i]->EvaluateBoundaryConditions(time, varName);
            }
        }
        cnt += m_fields[0]->GetBndCondExpansions()[n]->GetExpSize();
    }
}

void ShallowWaterSystem::WallBoundary2D(
    int bcRegion, int cnt, Array<OneD, Array<OneD, NekDouble>> &Fwd)
{
    int i;
    int nvariables = 3;

    // Adjust the physical values of the trace to take
    // user defined boundaries into account
    int e, id1, id2, npts;

    for (e = 0; e < m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize();
         ++e)
    {
        npts = m_fields[0]
                   ->GetBndCondExpansions()[bcRegion]
                   ->GetExp(e)
                   ->GetNumPoints(0);
        id1 = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetPhys_Offset(e);
        id2 = m_fields[0]->GetTrace()->GetPhys_Offset(
            m_fields[0]->GetTraceMap()->GetBndCondIDToGlobalTraceID(cnt + e));

        switch (m_expdim)
        {
            case 1:
            {
                // negate the forward flux
                Vmath::Neg(npts, &Fwd[1][id2], 1);
            }
            break;
            case 2:
            {
                Array<OneD, NekDouble> tmp_n(npts);
                Array<OneD, NekDouble> tmp_t(npts);

                Vmath::Vmul(npts, &Fwd[1][id2], 1, &m_traceNormals[0][id2], 1,
                            &tmp_n[0], 1);
                Vmath::Vvtvp(npts, &Fwd[2][id2], 1, &m_traceNormals[1][id2], 1,
                             &tmp_n[0], 1, &tmp_n[0], 1);

                Vmath::Vmul(npts, &Fwd[1][id2], 1, &m_traceNormals[1][id2], 1,
                            &tmp_t[0], 1);
                Vmath::Vvtvm(npts, &Fwd[2][id2], 1, &m_traceNormals[0][id2], 1,
                             &tmp_t[0], 1, &tmp_t[0], 1);

                // negate the normal flux
                Vmath::Neg(npts, tmp_n, 1);

                // rotate back to Cartesian
                Vmath::Vmul(npts, &tmp_t[0], 1, &m_traceNormals[1][id2], 1,
                            &Fwd[1][id2], 1);
                Vmath::Vvtvm(npts, &tmp_n[0], 1, &m_traceNormals[0][id2], 1,
                             &Fwd[1][id2], 1, &Fwd[1][id2], 1);

                Vmath::Vmul(npts, &tmp_t[0], 1, &m_traceNormals[0][id2], 1,
                            &Fwd[2][id2], 1);
                Vmath::Vvtvp(npts, &tmp_n[0], 1, &m_traceNormals[1][id2], 1,
                             &Fwd[2][id2], 1, &Fwd[2][id2], 1);
            }
            break;
            case 3:
                ASSERTL0(false,
                         "3D not implemented for Shallow Water Equations");
                break;
            default:
                ASSERTL0(false, "Illegal expansion dimension");
        }

        // copy boundary adjusted values into the boundary expansion
        for (i = 0; i < nvariables; ++i)
        {
            Vmath::Vcopy(npts, &Fwd[i][id2], 1,
                         &(m_fields[i]
                               ->GetBndCondExpansions()[bcRegion]
                               ->UpdatePhys())[id1],
                         1);
        }
    }
}

void ShallowWaterSystem::WallBoundary(
    int bcRegion, int cnt, Array<OneD, Array<OneD, NekDouble>> &Fwd,
    Array<OneD, Array<OneD, NekDouble>> &physarray)
{
    int i;
    int nvariables = physarray.size();

    // Adjust the physical values of the trace to take
    // user defined boundaries into account
    int e, id1, id2, npts;

    for (e = 0; e < m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize();
         ++e)
    {
        npts = m_fields[0]
                   ->GetBndCondExpansions()[bcRegion]
                   ->GetExp(e)
                   ->GetTotPoints();
        id1 = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetPhys_Offset(e);
        id2 = m_fields[0]->GetTrace()->GetPhys_Offset(
            m_fields[0]->GetTraceMap()->GetBndCondIDToGlobalTraceID(cnt + e));

        // For 2D/3D, define: v* = v - 2(v.n)n
        Array<OneD, NekDouble> tmp(npts, 0.0);

        // Calculate (v.n)
        for (i = 0; i < m_spacedim; ++i)
        {
            Vmath::Vvtvp(npts, &Fwd[1 + i][id2], 1, &m_traceNormals[i][id2], 1,
                         &tmp[0], 1, &tmp[0], 1);
        }

        // Calculate 2.0(v.n)
        Vmath::Smul(npts, -2.0, &tmp[0], 1, &tmp[0], 1);

        // Calculate v* = v - 2.0(v.n)n
        for (i = 0; i < m_spacedim; ++i)
        {
            Vmath::Vvtvp(npts, &tmp[0], 1, &m_traceNormals[i][id2], 1,
                         &Fwd[1 + i][id2], 1, &Fwd[1 + i][id2], 1);
        }

        // copy boundary adjusted values into the boundary expansion
        for (i = 0; i < nvariables; ++i)
        {
            Vmath::Vcopy(npts, &Fwd[i][id2], 1,
                         &(m_fields[i]
                               ->GetBndCondExpansions()[bcRegion]
                               ->UpdatePhys())[id1],
                         1);
        }
    }
}

// physarray contains the conservative variables
void ShallowWaterSystem::AddCoriolis(
    const Array<OneD, const Array<OneD, NekDouble>> &physarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{

    int ncoeffs = GetNcoeffs();
    int nq      = GetTotPoints();

    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> mod(ncoeffs);

    switch (m_projectionType)
    {
        case MultiRegions::eDiscontinuous:
        {
            // add to u equation
            Vmath::Vmul(nq, m_coriolis, 1, physarray[2], 1, tmp, 1);
            m_fields[0]->IProductWRTBase(tmp, mod);
            m_fields[0]->MultiplyByElmtInvMass(mod, mod);
            m_fields[0]->BwdTrans(mod, tmp);
            Vmath::Vadd(nq, tmp, 1, outarray[1], 1, outarray[1], 1);

            // add to v equation
            Vmath::Vmul(nq, m_coriolis, 1, physarray[1], 1, tmp, 1);
            Vmath::Neg(nq, tmp, 1);
            m_fields[0]->IProductWRTBase(tmp, mod);
            m_fields[0]->MultiplyByElmtInvMass(mod, mod);
            m_fields[0]->BwdTrans(mod, tmp);
            Vmath::Vadd(nq, tmp, 1, outarray[2], 1, outarray[2], 1);
        }
        break;
        case MultiRegions::eGalerkin:
        case MultiRegions::eMixed_CG_Discontinuous:
        {
            // add to u equation
            Vmath::Vmul(nq, m_coriolis, 1, physarray[2], 1, tmp, 1);
            Vmath::Vadd(nq, tmp, 1, outarray[1], 1, outarray[1], 1);

            // add to v equation
            Vmath::Vmul(nq, m_coriolis, 1, physarray[1], 1, tmp, 1);
            Vmath::Neg(nq, tmp, 1);
            Vmath::Vadd(nq, tmp, 1, outarray[2], 1, outarray[2], 1);
        }
        break;
        default:
            ASSERTL0(false, "Unknown projection scheme for the NonlinearSWE");
            break;
    }
}

void ShallowWaterSystem::ConservativeToPrimitive()
{
    int nq = GetTotPoints();

    // u = hu/h
    Vmath::Vdiv(nq, m_fields[1]->GetPhys(), 1, m_fields[0]->GetPhys(), 1,
                m_fields[1]->UpdatePhys(), 1);

    // v = hv/ v
    Vmath::Vdiv(nq, m_fields[2]->GetPhys(), 1, m_fields[0]->GetPhys(), 1,
                m_fields[2]->UpdatePhys(), 1);

    // \eta = h - d
    Vmath::Vsub(nq, m_fields[0]->GetPhys(), 1, m_depth, 1,
                m_fields[0]->UpdatePhys(), 1);
}

void ShallowWaterSystem::PrimitiveToConservative()
{
    int nq = GetTotPoints();

    // h = \eta + d
    Vmath::Vadd(nq, m_fields[0]->GetPhys(), 1, m_depth, 1,
                m_fields[0]->UpdatePhys(), 1);

    // hu = h * u
    Vmath::Vmul(nq, m_fields[0]->GetPhys(), 1, m_fields[1]->GetPhys(), 1,
                m_fields[1]->UpdatePhys(), 1);

    // hv = h * v
    Vmath::Vmul(nq, m_fields[0]->GetPhys(), 1, m_fields[2]->GetPhys(), 1,
                m_fields[2]->UpdatePhys(), 1);
}

void ShallowWaterSystem::EvaluateWaterDepth(void)
{
    GetFunction("WaterDepth")->Evaluate("d", m_depth);
}

void ShallowWaterSystem::EvaluateCoriolis(void)
{
    GetFunction("Coriolis")->Evaluate("f", m_coriolis);
}

void ShallowWaterSystem::CopyBoundaryTrace(const Array<OneD, NekDouble> &Fwd,
                                           Array<OneD, NekDouble> &Bwd)
{

    int cnt = 0;
    // loop over Boundary Regions
    for (int bcRegion = 0; bcRegion < m_fields[0]->GetBndConditions().size();
         ++bcRegion)
    {
        if (m_fields[0]
                ->GetBndConditions()[bcRegion]
                ->GetBoundaryConditionType() == SpatialDomains::ePeriodic)
        {
            continue;
        }

        // Copy the forward trace of the field to the backward trace
        int e, id2, npts;

        for (e = 0;
             e < m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize();
             ++e)
        {
            npts = m_fields[0]
                       ->GetBndCondExpansions()[bcRegion]
                       ->GetExp(e)
                       ->GetTotPoints();
            id2 = m_fields[0]->GetTrace()->GetPhys_Offset(
                m_fields[0]->GetTraceMap()->GetBndCondIDToGlobalTraceID(cnt +
                                                                        e));

            Vmath::Vcopy(npts, &Fwd[id2], 1, &Bwd[id2], 1);
        }

        cnt += m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize();
    }
}

} // namespace Nektar
