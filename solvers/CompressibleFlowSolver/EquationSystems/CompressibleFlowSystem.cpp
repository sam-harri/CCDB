///////////////////////////////////////////////////////////////////////////////
//
// File: CompressibleFlowSystem.cpp
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
// Description: Compressible flow system base class with auxiliary functions
//
///////////////////////////////////////////////////////////////////////////////

#include <CompressibleFlowSolver/EquationSystems/CompressibleFlowSystem.h>

#include <LibUtilities/BasicUtils/Timer.h>
#include <SolverUtils/Advection/AdvectionWeakDG.h>

using namespace std;

namespace Nektar
{
CompressibleFlowSystem::CompressibleFlowSystem(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : UnsteadySystem(pSession, pGraph), AdvectionSystem(pSession, pGraph)
{
}

/**
 * @brief Initialization object for CompressibleFlowSystem class.
 */
void CompressibleFlowSystem::v_InitObject(bool DeclareFields)
{
    AdvectionSystem::v_InitObject(DeclareFields);

    for (size_t i = 0; i < m_fields.size(); i++)
    {
        // Use BwdTrans to make sure initial condition is in solution space
        m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(),
                              m_fields[i]->UpdatePhys());
    }

    m_varConv = MemoryManager<VariableConverter>::AllocateSharedPtr(
        m_session, m_spacedim, m_graph);

    ASSERTL0(m_session->DefinesSolverInfo("UPWINDTYPE"),
             "No UPWINDTYPE defined in session.");

    // Do not forwards transform initial condition
    m_homoInitialFwd = false;

    // Set up locations of velocity vector.
    m_vecLocs    = Array<OneD, Array<OneD, NekDouble>>(1);
    m_vecLocs[0] = Array<OneD, NekDouble>(m_spacedim);
    for (int i = 0; i < m_spacedim; ++i)
    {
        m_vecLocs[0][i] = 1 + i;
    }

    // Loading parameters from session file
    InitialiseParameters();

    // Setting up advection and diffusion operators
    InitAdvection();

    // Create artificial diffusion with laplacian operator
    if (m_shockCaptureType == "NonSmooth")
    {
        m_artificialDiffusion = GetArtificialDiffusionFactory().CreateInstance(
            m_shockCaptureType, m_session, m_fields, m_spacedim);
    }

    // Forcing terms for the sponge region
    m_forcing = SolverUtils::Forcing::Load(m_session, shared_from_this(),
                                           m_fields, m_fields.size());

    // User-defined boundary conditions
    size_t cnt = 0;
    for (size_t n = 0; n < (size_t)m_fields[0]->GetBndConditions().size(); ++n)
    {
        std::string type = m_fields[0]->GetBndConditions()[n]->GetUserDefined();

        if (m_fields[0]->GetBndConditions()[n]->GetBoundaryConditionType() ==
            SpatialDomains::ePeriodic)
        {
            continue;
        }

        if (!type.empty())
        {
            m_bndConds.push_back(GetCFSBndCondFactory().CreateInstance(
                type, m_session, m_fields, m_traceNormals, m_gridVelocityTrace,
                m_spacedim, n, cnt));
        }
        cnt += m_fields[0]->GetBndCondExpansions()[n]->GetExpSize();
    }

    m_ode.DefineOdeRhs(&CompressibleFlowSystem::DoOdeRhs, this);
    m_ode.DefineProjection(&CompressibleFlowSystem::DoOdeProjection, this);

    SetBoundaryConditionsBwdWeight();
}

/**
 * @brief Load CFS parameters from the session file
 */
void CompressibleFlowSystem::InitialiseParameters()
{
    // Get gamma parameter from session file.
    m_session->LoadParameter("Gamma", m_gamma, 1.4);

    // Shock capture
    m_session->LoadSolverInfo("ShockCaptureType", m_shockCaptureType, "Off");

    // Check if the shock capture type is supported
    std::string err_msg = "Warning, ShockCaptureType = " + m_shockCaptureType +
                          " is not supported by this solver";
    ASSERTL0(v_SupportsShockCaptType(m_shockCaptureType), err_msg);

    // Load parameters for exponential filtering
    m_session->MatchSolverInfo("ExponentialFiltering", "True", m_useFiltering,
                               false);
    if (m_useFiltering)
    {
        m_session->LoadParameter("FilterAlpha", m_filterAlpha, 36);
        m_session->LoadParameter("FilterExponent", m_filterExponent, 16);
        m_session->LoadParameter("FilterCutoff", m_filterCutoff, 0);
    }

    // Load CFL for local time-stepping (for steady state)
    m_session->MatchSolverInfo("LocalTimeStep", "True", m_useLocalTimeStep,
                               false);
    if (m_useLocalTimeStep)
    {
        ASSERTL0(m_cflSafetyFactor != 0,
                 "Local time stepping requires CFL parameter.");
    }
}

/**
 * @brief Create advection and diffusion objects for CFS
 */
void CompressibleFlowSystem::InitAdvection()
{
    // Check if projection type is correct
    ASSERTL0(m_projectionType == MultiRegions::eDiscontinuous,
             "Unsupported projection type.");

    string advName, riemName;
    m_session->LoadSolverInfo("AdvectionType", advName, "WeakDG");

    m_advObject =
        SolverUtils::GetAdvectionFactory().CreateInstance(advName, advName);

    if (m_specHP_dealiasing)
    {
        m_advObject->SetFluxVector(
            &CompressibleFlowSystem::GetFluxVectorDeAlias, this);
    }
    else
    {
        m_advObject->SetFluxVector(&CompressibleFlowSystem::GetFluxVector,
                                   this);
    }

    // Setting up Riemann solver for advection operator
    m_session->LoadSolverInfo("UpwindType", riemName, "Average");

    SolverUtils::RiemannSolverSharedPtr riemannSolver;
    riemannSolver = SolverUtils::GetRiemannSolverFactory().CreateInstance(
        riemName, m_session);

    // Tell Riemann Solver if doing ALE and provide trace grid velocity
    riemannSolver->SetALEFlag(m_ALESolver);
    riemannSolver->SetVector("vgt", &ALEHelper::GetGridVelocityTrace, this);

    // Setting up parameters for advection operator Riemann solver
    riemannSolver->SetParam("gamma", &CompressibleFlowSystem::GetGamma, this);
    riemannSolver->SetAuxVec("vecLocs", &CompressibleFlowSystem::GetVecLocs,
                             this);
    riemannSolver->SetVector("N", &CompressibleFlowSystem::GetNormals, this);

    // Concluding initialisation of advection / diffusion operators
    m_advObject->SetRiemannSolver(riemannSolver);
    m_advObject->InitObject(m_session, m_fields);
}

/**
 * @brief Compute the right-hand side.
 */
void CompressibleFlowSystem::DoOdeRhs(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    LibUtilities::Timer timer;

    size_t nvariables = inarray.size();
    size_t nTracePts  = GetTraceTotPoints();

    // This converts our Mu in coefficient space to u in physical space for ALE
    Array<OneD, Array<OneD, NekDouble>> tmpIn(nvariables);
    if (m_ALESolver)
    {
        ALEHelper::ALEDoElmtInvMassBwdTrans(inarray, tmpIn);
    }
    else
    {
        tmpIn = inarray;
    }

    m_bndEvaluateTime = time;

    // Store forwards/backwards space along trace space
    Array<OneD, Array<OneD, NekDouble>> Fwd(nvariables);
    Array<OneD, Array<OneD, NekDouble>> Bwd(nvariables);

    if (m_HomogeneousType == eHomogeneous1D)
    {
        Fwd = NullNekDoubleArrayOfArray;
        Bwd = NullNekDoubleArrayOfArray;
    }
    else
    {
        for (size_t i = 0; i < nvariables; ++i)
        {
            Fwd[i] = Array<OneD, NekDouble>(nTracePts, 0.0);
            Bwd[i] = Array<OneD, NekDouble>(nTracePts, 0.0);
            m_fields[i]->GetFwdBwdTracePhys(tmpIn[i], Fwd[i], Bwd[i]);
        }
    }

    // Calculate advection
    timer.Start();
    DoAdvection(tmpIn, outarray, time, Fwd, Bwd);
    timer.Stop();
    timer.AccumulateRegion("DoAdvection");

    // Negate results
    for (size_t i = 0; i < nvariables; ++i)
    {
        Vmath::Neg(outarray[i].size(), outarray[i], 1);
    }

    // Add diffusion terms
    timer.Start();
    DoDiffusion(tmpIn, outarray, Fwd, Bwd);
    timer.Stop();
    timer.AccumulateRegion("DoDiffusion");

    // Add forcing terms
    for (auto &x : m_forcing)
    {
        x->Apply(m_fields, tmpIn, outarray, time);
    }

    if (m_useLocalTimeStep)
    {
        size_t nElements = m_fields[0]->GetExpSize();
        int nq, offset;
        NekDouble fac;
        Array<OneD, NekDouble> tmp;

        Array<OneD, NekDouble> tstep(nElements, 0.0);
        GetElmtTimeStep(tmpIn, tstep);

        // Loop over elements
        for (size_t n = 0; n < nElements; ++n)
        {
            nq     = m_fields[0]->GetExp(n)->GetTotPoints();
            offset = m_fields[0]->GetPhys_Offset(n);
            fac    = tstep[n] / m_timestep;
            for (size_t i = 0; i < nvariables; ++i)
            {
                Vmath::Smul(nq, fac, outarray[i] + offset, 1,
                            tmp = outarray[i] + offset, 1);
            }
        }
    }
}

/**
 * @brief Compute the projection and call the method for imposing the
 * boundary conditions in case of discontinuous projection.
 */
void CompressibleFlowSystem::DoOdeProjection(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    size_t nvariables = inarray.size();

    // Perform ALE movement
    if (m_ALESolver)
    {
        MoveMesh(time, m_traceNormals);
    }

    switch (m_projectionType)
    {
        case MultiRegions::eDiscontinuous:
        {
            // Just copy over array
            for (size_t i = 0; i < nvariables; ++i)
            {
                Vmath::Vcopy(inarray[i].size(), inarray[i], 1, outarray[i], 1);

                if (m_useFiltering)
                {
                    m_fields[i]->ExponentialFilter(outarray[i], m_filterAlpha,
                                                   m_filterExponent,
                                                   m_filterCutoff);
                }
            }
            SetBoundaryConditions(outarray, time);
            break;
        }
        case MultiRegions::eGalerkin:
        case MultiRegions::eMixed_CG_Discontinuous:
        {
            NEKERROR(ErrorUtil::efatal, "No Continuous Galerkin for full "
                                        "compressible Navier-Stokes equations");
            break;
        }
        default:
            NEKERROR(ErrorUtil::efatal, "Unknown projection scheme");
            break;
    }
}

/**
 * @brief Compute the advection terms for the right-hand side
 */
void CompressibleFlowSystem::DoAdvection(
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
    const Array<OneD, Array<OneD, NekDouble>> &pFwd,
    const Array<OneD, Array<OneD, NekDouble>> &pBwd)
{
    int nvariables = inarray.size();
    Array<OneD, Array<OneD, NekDouble>> advVel(m_spacedim);

    if (m_ALESolver)
    {
        auto advWeakDGObject =
            std::dynamic_pointer_cast<SolverUtils::AdvectionWeakDG>(
                m_advObject);
        advWeakDGObject->AdvectCoeffs(nvariables, m_fields, advVel, inarray,
                                      outarray, time, pFwd, pBwd);
    }
    else
    {
        m_advObject->Advect(nvariables, m_fields, advVel, inarray, outarray,
                            time, pFwd, pBwd);
    }
}

/**
 * @brief Add the diffusions terms to the right-hand side
 */
void CompressibleFlowSystem::DoDiffusion(
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray,
    const Array<OneD, Array<OneD, NekDouble>> &pFwd,
    const Array<OneD, Array<OneD, NekDouble>> &pBwd)
{
    v_DoDiffusion(inarray, outarray, pFwd, pBwd);
}

void CompressibleFlowSystem::SetBoundaryConditions(
    Array<OneD, Array<OneD, NekDouble>> &physarray, NekDouble time)
{
    size_t nTracePts  = GetTraceTotPoints();
    size_t nvariables = physarray.size();

    // This converts our Mu in coefficient space to u in physical space for ALE
    Array<OneD, Array<OneD, NekDouble>> tmpIn(nvariables);

    if (m_ALESolver && !m_ImplicitALESolver)
    {
        ALEHelper::ALEDoElmtInvMassBwdTrans(physarray, tmpIn);
    }
    else
    {
        tmpIn = physarray;
    }

    Array<OneD, Array<OneD, NekDouble>> Fwd(nvariables);
    for (size_t i = 0; i < nvariables; ++i)
    {
        Fwd[i] = Array<OneD, NekDouble>(nTracePts);
        m_fields[i]->ExtractTracePhys(tmpIn[i], Fwd[i]);
    }

    if (!m_bndConds.empty())
    {
        // Loop over user-defined boundary conditions
        for (auto &x : m_bndConds)
        {
            x->Apply(Fwd, tmpIn, time);
        }
    }
}
/**
 * @brief Set up a weight on physical boundaries for boundary condition
 * applications
 */
void CompressibleFlowSystem::SetBoundaryConditionsBwdWeight()
{
    if (m_bndConds.size())
    {
        // Loop over user-defined boundary conditions
        for (auto &x : m_bndConds)
        {
            x->ApplyBwdWeight();
        }
    }
}

/**
 * @brief Return the flux vector for the compressible Euler equations.
 *
 * @param physfield   Fields.
 * @param flux        Resulting flux.
 */
void CompressibleFlowSystem::GetFluxVector(
    const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    TensorOfArray3D<NekDouble> &flux)
{
    size_t nVariables = physfield.size();
    size_t nPts       = physfield[0].size();

    constexpr unsigned short maxVel = 3;
    constexpr unsigned short maxFld = 5;

    // hardcoding done for performance reasons
    ASSERTL1(nVariables <= maxFld, "GetFluxVector, hard coded max fields");

    for (size_t p = 0; p < nPts; ++p)
    {
        // local storage
        std::array<NekDouble, maxFld> fieldTmp;
        std::array<NekDouble, maxVel> velocity;

        // rearrenge and load data
        for (size_t f = 0; f < nVariables; ++f)
        {
            fieldTmp[f] = physfield[f][p]; // load
        }

        // 1 / rho
        NekDouble oneOrho = 1.0 / fieldTmp[0];

        for (size_t d = 0; d < m_spacedim; ++d)
        {
            // Flux vector for the rho equation
            flux[0][d][p] = fieldTmp[d + 1]; // store
            // compute velocity
            velocity[d] = fieldTmp[d + 1] * oneOrho;
        }

        NekDouble pressure = m_varConv->GetPressure(fieldTmp.data());
        NekDouble ePlusP   = fieldTmp[m_spacedim + 1] + pressure;
        for (size_t f = 0; f < m_spacedim; ++f)
        {
            // Flux vector for the velocity fields
            for (size_t d = 0; d < m_spacedim; ++d)
            {
                flux[f + 1][d][p] = velocity[d] * fieldTmp[f + 1]; // store
            }

            // Add pressure to appropriate field
            flux[f + 1][f][p] += pressure;

            // Flux vector for energy
            flux[m_spacedim + 1][f][p] = ePlusP * velocity[f]; // store
        }
    }

    // @TODO : for each row (3 columns) negative grid velocity component (d * d
    // + 2 (4 rows)) rho, rhou, rhov, rhow, E,
    // @TODO : top row is flux for rho etc... each row subtract v_g * conserved
    // variable for that row... For grid velocity subtract v_g * conserved
    // variable
    if (m_ALESolver)
    {
        for (int i = 0; i < m_spacedim + 2; ++i)
        {
            for (int j = 0; j < m_spacedim; ++j)
            {
                for (int k = 0; k < nPts; ++k)
                {
                    flux[i][j][k] -= physfield[i][k] * m_gridVelocity[j][k];
                }
            }
        }
    }
}

/**
 * @brief Return the flux vector for the compressible Euler equations
 * by using the de-aliasing technique.
 *
 * @param physfield   Fields.
 * @param flux        Resulting flux.
 */
void CompressibleFlowSystem::GetFluxVectorDeAlias(
    const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    TensorOfArray3D<NekDouble> &flux)
{
    int i, j;
    size_t nq         = physfield[0].size();
    size_t nVariables = m_fields.size();

    // Factor to rescale 1d points in dealiasing
    NekDouble OneDptscale = 2;
    nq                    = m_fields[0]->Get1DScaledTotPoints(OneDptscale);

    Array<OneD, NekDouble> pressure(nq);
    Array<OneD, Array<OneD, NekDouble>> velocity(m_spacedim);

    Array<OneD, Array<OneD, NekDouble>> physfield_interp(nVariables);
    TensorOfArray3D<NekDouble> flux_interp(nVariables);

    for (size_t i = 0; i < nVariables; ++i)
    {
        physfield_interp[i] = Array<OneD, NekDouble>(nq);
        flux_interp[i]      = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
        m_fields[0]->PhysInterp1DScaled(OneDptscale, physfield[i],
                                        physfield_interp[i]);

        for (j = 0; j < m_spacedim; ++j)
        {
            flux_interp[i][j] = Array<OneD, NekDouble>(nq);
        }
    }

    // Flux vector for the rho equation
    for (i = 0; i < m_spacedim; ++i)
    {
        velocity[i] = Array<OneD, NekDouble>(nq);

        // Galerkin project solution back to original space
        m_fields[0]->PhysGalerkinProjection1DScaled(
            OneDptscale, physfield_interp[i + 1], flux[0][i]);
    }

    m_varConv->GetVelocityVector(physfield_interp, velocity);
    m_varConv->GetPressure(physfield_interp, pressure);

    // Evaluation of flux vector for the velocity fields
    for (i = 0; i < m_spacedim; ++i)
    {
        for (j = 0; j < m_spacedim; ++j)
        {
            Vmath::Vmul(nq, velocity[j], 1, physfield_interp[i + 1], 1,
                        flux_interp[i + 1][j], 1);
        }

        // Add pressure to appropriate field
        Vmath::Vadd(nq, flux_interp[i + 1][i], 1, pressure, 1,
                    flux_interp[i + 1][i], 1);
    }

    // Galerkin project solution back to original space
    for (i = 0; i < m_spacedim; ++i)
    {
        for (j = 0; j < m_spacedim; ++j)
        {
            m_fields[0]->PhysGalerkinProjection1DScaled(
                OneDptscale, flux_interp[i + 1][j], flux[i + 1][j]);
        }
    }

    // Evaluation of flux vector for energy
    Vmath::Vadd(nq, physfield_interp[m_spacedim + 1], 1, pressure, 1, pressure,
                1);

    for (j = 0; j < m_spacedim; ++j)
    {
        Vmath::Vmul(nq, velocity[j], 1, pressure, 1,
                    flux_interp[m_spacedim + 1][j], 1);

        // Galerkin project solution back to original space
        m_fields[0]->PhysGalerkinProjection1DScaled(
            OneDptscale, flux_interp[m_spacedim + 1][j],
            flux[m_spacedim + 1][j]);
    }
}

/**
 * @brief Calculate the maximum timestep on each element
 *        subject to CFL restrictions.
 */
void CompressibleFlowSystem::GetElmtTimeStep(
    [[maybe_unused]] const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, NekDouble> &tstep)
{
    size_t nElements = m_fields[0]->GetExpSize();

    // Change value of m_timestep (in case it is set to zero)
    NekDouble tmp = m_timestep;
    m_timestep    = 1.0;

    Array<OneD, NekDouble> cfl(nElements);
    cfl = GetElmtCFLVals();

    // Factors to compute the time-step limit
    NekDouble alpha = MaxTimeStepEstimator();

    // Loop over elements to compute the time-step limit for each element
    for (size_t n = 0; n < nElements; ++n)
    {
        tstep[n] = m_cflSafetyFactor * alpha / cfl[n];
    }

    // Restore value of m_timestep
    m_timestep = tmp;
}

/**
 * @brief Calculate the maximum timestep subject to CFL restrictions.
 */
NekDouble CompressibleFlowSystem::v_GetTimeStep(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray)
{
    int nElements = m_fields[0]->GetExpSize();
    Array<OneD, NekDouble> tstep(nElements, 0.0);

    GetElmtTimeStep(inarray, tstep);

    // Get the minimum time-step limit and return the time-step
    NekDouble TimeStep = Vmath::Vmin(nElements, tstep, 1);
    m_comm->GetSpaceComm()->AllReduce(TimeStep, LibUtilities::ReduceMin);

    return TimeStep;
}

void CompressibleFlowSystem::v_GenerateSummary(SolverUtils::SummaryList &s)
{
    UnsteadySystem::v_GenerateSummary(s);
    if (m_session->DefinesSolverInfo("ICTYPE"))
    {
        if (boost::iequals(m_session->GetSolverInfo("ICTYPE"),
                           "IsentropicVortex"))
        {
            SolverUtils::AddSummaryItem(s, "Problem Type", "IsentropicVortex");
        }
        else if (boost::iequals(m_session->GetSolverInfo("ICTYPE"),
                                "RinglebFlow"))
        {
            SolverUtils::AddSummaryItem(s, "Problem Type", "RinglebFlow");
        }
        else
        {
            NEKERROR(ErrorUtil::efatal, "unknow initial condition");
        }
    }
}

/**
 * @brief Set up logic for residual calculation.
 */
void CompressibleFlowSystem::v_SetInitialConditions(
    NekDouble initialtime, bool dumpInitialConditions,
    [[maybe_unused]] const int domain)
{
    if (m_session->DefinesSolverInfo("ICTYPE") &&
        boost::iequals(m_session->GetSolverInfo("ICTYPE"), "IsentropicVortex"))
    {
        // Forward transform to fill the coefficient space
        for (int i = 0; i < m_fields.size(); ++i)
        {
            EvaluateIsentropicVortex(i, m_fields[i]->UpdatePhys(), initialtime);
            m_fields[i]->SetPhysState(true);
            m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),
                                  m_fields[i]->UpdateCoeffs());
        }
    }
    else
    {
        EquationSystem::v_SetInitialConditions(initialtime, false);
        m_nchk--; // Note: m_nchk has been incremented in EquationSystem.

        // insert white noise in initial condition
        NekDouble Noise;
        int phystot = m_fields[0]->GetTotPoints();
        Array<OneD, NekDouble> noise(phystot);

        m_session->LoadParameter("Noise", Noise, 0.0);
        int m_nConvectiveFields = m_fields.size();

        if (Noise > 0.0)
        {
            int seed = -m_comm->GetSpaceComm()->GetRank() * m_nConvectiveFields;
            for (int i = 0; i < m_nConvectiveFields; i++)
            {
                Vmath::FillWhiteNoise(phystot, Noise, noise, 1, seed);
                --seed;
                Vmath::Vadd(phystot, m_fields[i]->GetPhys(), 1, noise, 1,
                            m_fields[i]->UpdatePhys(), 1);
                m_fields[i]->FwdTransLocalElmt(m_fields[i]->GetPhys(),
                                               m_fields[i]->UpdateCoeffs());
            }
        }
    }

    if (dumpInitialConditions && m_nchk == 0 && m_checksteps &&
        !m_comm->IsParallelInTime())
    {
        Checkpoint_Output(0);
    }
    else if (dumpInitialConditions && m_nchk == 0 && m_comm->IsParallelInTime())
    {
        std::string newdir = m_sessionName + ".pit";
        if (!fs::is_directory(newdir))
        {
            fs::create_directory(newdir);
        }
        if (m_comm->GetTimeComm()->GetRank() == 0)
        {
            WriteFld(newdir + "/" + m_sessionName + "_0.fld");
        }
    }
    m_nchk++;
}

void CompressibleFlowSystem::v_EvaluateExactSolution(
    unsigned int field, Array<OneD, NekDouble> &outfield, const NekDouble time)
{

    if (!m_session->DefinesFunction("ExactSolution") &&
        m_session->DefinesSolverInfo("ICTYPE"))
    {
        if (boost::iequals(m_session->GetSolverInfo("ICTYPE"),
                           "IsentropicVortex"))
        {
            EvaluateIsentropicVortex(field, outfield, time);
        }
        else if (boost::iequals(m_session->GetSolverInfo("ICTYPE"),
                                "RinglebFlow"))
        {
            GetExactRinglebFlow(field, outfield);
        }
    }
    else
    {
        EquationSystem::v_EvaluateExactSolution(field, outfield, time);
    }
}

/**
 * @brief Compute the advection velocity in the standard space
 * for each element of the expansion.
 */
Array<OneD, NekDouble> CompressibleFlowSystem::v_GetMaxStdVelocity(
    const NekDouble SpeedSoundFactor)
{
    size_t nTotQuadPoints = GetTotPoints();
    size_t n_element      = m_fields[0]->GetExpSize();
    size_t expdim         = m_fields[0]->GetGraph()->GetMeshDimension();
    size_t nfields        = m_fields.size();
    int offset;
    Array<OneD, NekDouble> tmp;

    Array<OneD, Array<OneD, NekDouble>> physfields(nfields);
    for (size_t i = 0; i < nfields; ++i)
    {
        physfields[i] = m_fields[i]->GetPhys();
    }

    Array<OneD, NekDouble> stdV(n_element, 0.0);

    // Getting the velocity vector on the 2D normal space
    Array<OneD, Array<OneD, NekDouble>> velocity(m_spacedim);
    Array<OneD, Array<OneD, NekDouble>> stdVelocity(m_spacedim);
    Array<OneD, Array<OneD, NekDouble>> stdSoundSpeed(m_spacedim);
    Array<OneD, NekDouble> soundspeed(nTotQuadPoints);
    LibUtilities::PointsKeyVector ptsKeys;

    for (int i = 0; i < m_spacedim; ++i)
    {
        velocity[i]      = Array<OneD, NekDouble>(nTotQuadPoints);
        stdVelocity[i]   = Array<OneD, NekDouble>(nTotQuadPoints, 0.0);
        stdSoundSpeed[i] = Array<OneD, NekDouble>(nTotQuadPoints, 0.0);
    }

    m_varConv->GetVelocityVector(physfields, velocity);
    m_varConv->GetSoundSpeed(physfields, soundspeed);

    // Subtract Ug from the velocity for the ALE formulation
    for (size_t i = 0; i < m_spacedim; ++i)
    {
        Vmath::Vsub(nTotQuadPoints, velocity[i], 1, m_gridVelocity[i], 1,
                    velocity[i], 1);
    }

    for (int el = 0; el < n_element; ++el)
    {
        ptsKeys = m_fields[0]->GetExp(el)->GetPointsKeys();
        offset  = m_fields[0]->GetPhys_Offset(el);
        int nq  = m_fields[0]->GetExp(el)->GetTotPoints();

        const SpatialDomains::GeomFactorsSharedPtr metricInfo =
            m_fields[0]->GetExp(el)->GetGeom()->GetMetricInfo();
        const Array<TwoD, const NekDouble> &gmat =
            m_fields[0]
                ->GetExp(el)
                ->GetGeom()
                ->GetMetricInfo()
                ->GetDerivFactors(ptsKeys);

        // Convert to standard element
        //    consider soundspeed in all directions
        //    (this might overestimate the cfl)
        if (metricInfo->GetGtype() == SpatialDomains::eDeformed)
        {
            // d xi/ dx = gmat = 1/J * d x/d xi
            for (size_t i = 0; i < expdim; ++i)
            {
                Vmath::Vmul(nq, gmat[i], 1, velocity[0] + offset, 1,
                            tmp = stdVelocity[i] + offset, 1);
                Vmath::Vmul(nq, gmat[i], 1, soundspeed + offset, 1,
                            tmp = stdSoundSpeed[i] + offset, 1);
                for (size_t j = 1; j < expdim; ++j)
                {
                    Vmath::Vvtvp(nq, gmat[expdim * j + i], 1,
                                 velocity[j] + offset, 1,
                                 stdVelocity[i] + offset, 1,
                                 tmp = stdVelocity[i] + offset, 1);
                    Vmath::Vvtvp(nq, gmat[expdim * j + i], 1,
                                 soundspeed + offset, 1,
                                 stdSoundSpeed[i] + offset, 1,
                                 tmp = stdSoundSpeed[i] + offset, 1);
                }
            }
        }
        else
        {
            for (size_t i = 0; i < expdim; ++i)
            {
                Vmath::Smul(nq, gmat[i][0], velocity[0] + offset, 1,
                            tmp = stdVelocity[i] + offset, 1);
                Vmath::Smul(nq, gmat[i][0], soundspeed + offset, 1,
                            tmp = stdSoundSpeed[i] + offset, 1);
                for (size_t j = 1; j < expdim; ++j)
                {
                    Vmath::Svtvp(nq, gmat[expdim * j + i][0],
                                 velocity[j] + offset, 1,
                                 stdVelocity[i] + offset, 1,
                                 tmp = stdVelocity[i] + offset, 1);
                    Vmath::Svtvp(nq, gmat[expdim * j + i][0],
                                 soundspeed + offset, 1,
                                 stdSoundSpeed[i] + offset, 1,
                                 tmp = stdSoundSpeed[i] + offset, 1);
                }
            }
        }

        NekDouble vel;
        for (size_t i = 0; i < nq; ++i)
        {
            NekDouble pntVelocity = 0.0;
            for (size_t j = 0; j < expdim; ++j)
            {
                // Add sound speed
                vel = std::abs(stdVelocity[j][offset + i]) +
                      SpeedSoundFactor * std::abs(stdSoundSpeed[j][offset + i]);
                pntVelocity += vel * vel;
            }
            pntVelocity = sqrt(pntVelocity);
            if (pntVelocity > stdV[el])
            {
                stdV[el] = pntVelocity;
            }
        }
    }

    return stdV;
}

/**
 * @brief Set the denominator to compute the time step when a cfl
 * control is employed. This function is no longer used but is still
 * here for being utilised in the future.
 *
 * @param n   Order of expansion element by element.
 */
NekDouble CompressibleFlowSystem::GetStabilityLimit(int n)
{
    ASSERTL0(n <= 20, "Illegal modes dimension for CFL calculation "
                      "(P has to be less then 20)");

    NekDouble CFLDG[21] = {2.0000,   6.0000,   11.8424,  19.1569,  27.8419,
                           37.8247,  49.0518,  61.4815,  75.0797,  89.8181,
                           105.6700, 122.6200, 140.6400, 159.7300, 179.8500,
                           201.0100, 223.1800, 246.3600, 270.5300, 295.6900,
                           321.8300}; // CFLDG 1D [0-20]
    NekDouble CFL       = 0.0;

    if (m_projectionType == MultiRegions::eDiscontinuous)
    {
        CFL = CFLDG[n];
    }
    else
    {
        NEKERROR(ErrorUtil::efatal, "Continuous Galerkin stability "
                                    "coefficients not introduced yet.");
    }

    return CFL;
}

/**
 * @brief Compute the vector of denominators to compute the time step
 * when a cfl control is employed. This function is no longer used but
 * is still here for being utilised in the future.
 *
 * @param ExpOrder   Order of expansion element by element.
 */
Array<OneD, NekDouble> CompressibleFlowSystem::GetStabilityLimitVector(
    const Array<OneD, int> &ExpOrder)
{
    int i;
    Array<OneD, NekDouble> returnval(m_fields[0]->GetExpSize(), 0.0);
    for (i = 0; i < m_fields[0]->GetExpSize(); i++)
    {
        returnval[i] = GetStabilityLimit(ExpOrder[i]);
    }
    return returnval;
}

void CompressibleFlowSystem::v_ExtraFldOutput(
    std::vector<Array<OneD, NekDouble>> &fieldcoeffs,
    std::vector<std::string> &variables)
{
    bool extraFields;
    m_session->MatchSolverInfo("OutputExtraFields", "True", extraFields, true);
    if (extraFields)
    {
        const int nPhys   = m_fields[0]->GetNpoints();
        const int nCoeffs = m_fields[0]->GetNcoeffs();
        Array<OneD, Array<OneD, NekDouble>> tmp(m_fields.size());

        for (int i = 0; i < m_fields.size(); ++i)
        {
            tmp[i] = m_fields[i]->GetPhys();
        }

        Array<OneD, Array<OneD, NekDouble>> velocity(m_spacedim);
        Array<OneD, Array<OneD, NekDouble>> velFwd(m_spacedim);
        for (int i = 0; i < m_spacedim; ++i)
        {
            velocity[i] = Array<OneD, NekDouble>(nPhys);
            velFwd[i]   = Array<OneD, NekDouble>(nCoeffs);
        }

        Array<OneD, NekDouble> pressure(nPhys), temperature(nPhys);
        Array<OneD, NekDouble> entropy(nPhys);
        Array<OneD, NekDouble> soundspeed(nPhys), mach(nPhys);
        Array<OneD, NekDouble> sensor(nPhys), SensorKappa(nPhys);

        m_varConv->GetVelocityVector(tmp, velocity);
        m_varConv->GetPressure(tmp, pressure);
        m_varConv->GetTemperature(tmp, temperature);
        m_varConv->GetEntropy(tmp, entropy);
        m_varConv->GetSoundSpeed(tmp, soundspeed);
        m_varConv->GetMach(tmp, soundspeed, mach);

        int sensorOffset;
        m_session->LoadParameter("SensorOffset", sensorOffset, 1);
        m_varConv->GetSensor(m_fields[0], tmp, sensor, SensorKappa,
                             sensorOffset);

        Array<OneD, NekDouble> pFwd(nCoeffs), TFwd(nCoeffs);
        Array<OneD, NekDouble> sFwd(nCoeffs);
        Array<OneD, NekDouble> aFwd(nCoeffs), mFwd(nCoeffs);
        Array<OneD, NekDouble> sensFwd(nCoeffs);

        string velNames[3] = {"u", "v", "w"};
        for (int i = 0; i < m_spacedim; ++i)
        {
            m_fields[0]->FwdTransLocalElmt(velocity[i], velFwd[i]);
            variables.push_back(velNames[i]);
            fieldcoeffs.push_back(velFwd[i]);
        }

        m_fields[0]->FwdTransLocalElmt(pressure, pFwd);
        m_fields[0]->FwdTransLocalElmt(temperature, TFwd);
        m_fields[0]->FwdTransLocalElmt(entropy, sFwd);
        m_fields[0]->FwdTransLocalElmt(soundspeed, aFwd);
        m_fields[0]->FwdTransLocalElmt(mach, mFwd);
        m_fields[0]->FwdTransLocalElmt(sensor, sensFwd);

        variables.push_back("p");
        variables.push_back("T");
        variables.push_back("s");
        variables.push_back("a");
        variables.push_back("Mach");
        variables.push_back("Sensor");
        fieldcoeffs.push_back(pFwd);
        fieldcoeffs.push_back(TFwd);
        fieldcoeffs.push_back(sFwd);
        fieldcoeffs.push_back(aFwd);
        fieldcoeffs.push_back(mFwd);
        fieldcoeffs.push_back(sensFwd);

        if (m_artificialDiffusion)
        {
            // reuse pressure
            Array<OneD, NekDouble> sensorFwd(nCoeffs);
            m_artificialDiffusion->GetArtificialViscosity(tmp, pressure);
            m_fields[0]->FwdTransLocalElmt(pressure, sensorFwd);

            variables.push_back("ArtificialVisc");
            fieldcoeffs.push_back(sensorFwd);
        }

        if (m_ALESolver)
        {
            ExtraFldOutputGridVelocity(fieldcoeffs, variables);
        }
    }
}

/**
 *
 */
void CompressibleFlowSystem::v_GetPressure(
    const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, NekDouble> &pressure)
{
    m_varConv->GetPressure(physfield, pressure);
}

/**
 *
 */
void CompressibleFlowSystem::v_GetDensity(
    const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, NekDouble> &density)
{
    density = physfield[0];
}

/**
 *
 */
void CompressibleFlowSystem::v_GetVelocity(
    const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, Array<OneD, NekDouble>> &velocity)
{
    m_varConv->GetVelocityVector(physfield, velocity);
}

void CompressibleFlowSystem::v_SteadyStateResidual([[maybe_unused]] int step,
                                                   Array<OneD, NekDouble> &L2)
{
    const size_t nPoints = GetTotPoints();
    const size_t nFields = m_fields.size();
    Array<OneD, Array<OneD, NekDouble>> rhs(nFields);
    Array<OneD, Array<OneD, NekDouble>> inarray(nFields);
    for (size_t i = 0; i < nFields; ++i)
    {
        rhs[i]     = Array<OneD, NekDouble>(nPoints, 0.0);
        inarray[i] = m_fields[i]->UpdatePhys();
    }

    DoOdeRhs(inarray, rhs, m_time);

    // Holds L2 errors.
    Array<OneD, NekDouble> tmp;
    Array<OneD, NekDouble> RHSL2(nFields);
    Array<OneD, NekDouble> residual(nFields);

    for (size_t i = 0; i < nFields; ++i)
    {
        tmp = rhs[i];

        Vmath::Vmul(nPoints, tmp, 1, tmp, 1, tmp, 1);
        residual[i] = Vmath::Vsum(nPoints, tmp, 1);
    }

    m_comm->GetSpaceComm()->AllReduce(residual, LibUtilities::ReduceSum);

    NekDouble onPoints = 1.0 / NekDouble(nPoints);
    for (size_t i = 0; i < nFields; ++i)
    {
        L2[i] = sqrt(residual[i] * onPoints);
    }
}

/**
 *
 */
void CompressibleFlowSystem::EvaluateIsentropicVortex(
    unsigned int field, Array<OneD, NekDouble> &outfield, NekDouble time,
    const int o)
{
    NekDouble beta, u0, v0, x0, y0;

    int nTotQuadPoints = GetTotPoints();
    Array<OneD, NekDouble> x(nTotQuadPoints);
    Array<OneD, NekDouble> y(nTotQuadPoints);
    Array<OneD, NekDouble> z(nTotQuadPoints);
    Array<OneD, Array<OneD, NekDouble>> u(m_spacedim + 2);

    m_fields[0]->GetCoords(x, y, z);

    for (int i = 0; i < m_spacedim + 2; ++i)
    {
        u[i] = Array<OneD, NekDouble>(nTotQuadPoints);
    }
    m_session->LoadParameter("IsentropicBeta", beta, 5.0);
    m_session->LoadParameter("IsentropicU0", u0, 1.0);
    m_session->LoadParameter("IsentropicV0", v0, 0.5);
    m_session->LoadParameter("IsentropicX0", x0, 5.0);
    m_session->LoadParameter("IsentropicY0", y0, 0.0);

    // Flow parameters
    NekDouble r, xbar, ybar, tmp;
    NekDouble fac = 1.0 / (16.0 * m_gamma * M_PI * M_PI);

    // In 3D zero rhow field.
    if (m_spacedim == 3)
    {
        Vmath::Zero(nTotQuadPoints, &u[3][o], 1);
    }

    // Fill storage
    for (int i = 0; i < nTotQuadPoints; ++i)
    {
        xbar = x[i] - u0 * time - x0;
        ybar = y[i] - v0 * time - y0;
        r    = sqrt(xbar * xbar + ybar * ybar);
        tmp  = beta * exp(1 - r * r);
        u[0][i + o] =
            pow(1.0 - (m_gamma - 1.0) * tmp * tmp * fac, 1.0 / (m_gamma - 1.0));
        u[1][i + o] = u[0][i + o] * (u0 - tmp * ybar / (2 * M_PI));
        u[2][i + o] = u[0][i + o] * (v0 + tmp * xbar / (2 * M_PI));
        u[m_spacedim + 1][i + o] =
            pow(u[0][i + o], m_gamma) / (m_gamma - 1.0) +
            0.5 * (u[1][i + o] * u[1][i + o] + u[2][i + o] * u[2][i + o]) /
                u[0][i + o];
    }
    Vmath::Vcopy(nTotQuadPoints, u[field].data(), 1, outfield.data(), 1);
}

/**
 * @brief Compute the exact solution for the Ringleb flow problem.
 */
void CompressibleFlowSystem::GetExactRinglebFlow(
    int field, Array<OneD, NekDouble> &outarray)
{
    int nTotQuadPoints = GetTotPoints();

    Array<OneD, NekDouble> rho(nTotQuadPoints, 100.0);
    Array<OneD, NekDouble> rhou(nTotQuadPoints);
    Array<OneD, NekDouble> rhov(nTotQuadPoints);
    Array<OneD, NekDouble> E(nTotQuadPoints);
    Array<OneD, NekDouble> x(nTotQuadPoints);
    Array<OneD, NekDouble> y(nTotQuadPoints);
    Array<OneD, NekDouble> z(nTotQuadPoints);

    m_fields[0]->GetCoords(x, y, z);

    // Flow parameters
    NekDouble c, k, phi, r, J, VV, pp, sint, P, ss;
    NekDouble J11, J12, J21, J22, det;
    NekDouble Fx, Fy;
    NekDouble xi, yi;
    NekDouble dV;
    NekDouble dtheta;
    NekDouble par1;
    NekDouble theta     = M_PI / 4.0;
    NekDouble kExt      = 0.7;
    NekDouble V         = kExt * sin(theta);
    NekDouble toll      = 1.0e-8;
    NekDouble errV      = 1.0;
    NekDouble errTheta  = 1.0;
    NekDouble gamma     = m_gamma;
    NekDouble gamma_1_2 = (gamma - 1.0) / 2.0;

    for (int i = 0; i < nTotQuadPoints; ++i)
    {
        while ((abs(errV) > toll) || (abs(errTheta) > toll))
        {
            VV   = V * V;
            sint = sin(theta);
            c    = sqrt(1.0 - gamma_1_2 * VV);
            k    = V / sint;
            phi  = 1.0 / k;
            pp   = phi * phi;
            J    = 1.0 / c + 1.0 / (3.0 * c * c * c) +
                1.0 / (5.0 * c * c * c * c * c) -
                0.5 * log((1.0 + c) / (1.0 - c));

            r    = pow(c, 1.0 / gamma_1_2);
            xi   = 1.0 / (2.0 * r) * (1.0 / VV - 2.0 * pp) + J / 2.0;
            yi   = phi / (r * V) * sqrt(1.0 - VV * pp);
            par1 = 25.0 - 5.0 * VV;
            ss   = sint * sint;

            Fx = xi - x[i];
            Fy = yi - y[i];

            J11 =
                39062.5 / pow(par1, 3.5) * (1.0 / VV - 2.0 / VV * ss) * V +
                1562.5 / pow(par1, 2.5) *
                    (-2.0 / (VV * V) + 4.0 / (VV * V) * ss) +
                12.5 / pow(par1, 1.5) * V + 312.5 / pow(par1, 2.5) * V +
                7812.5 / pow(par1, 3.5) * V -
                0.25 *
                    (-1.0 / pow(par1, 0.5) * V / (1.0 - 0.2 * pow(par1, 0.5)) -
                     (1.0 + 0.2 * pow(par1, 0.5)) /
                         pow((1.0 - 0.2 * pow(par1, 0.5)), 2.0) /
                         pow(par1, 0.5) * V) /
                    (1.0 + 0.2 * pow(par1, 0.5)) * (1.0 - 0.2 * pow(par1, 0.5));

            J12 = -6250.0 / pow(par1, 2.5) / VV * sint * cos(theta);
            J21 = -6250.0 / (VV * V) * sint / pow(par1, 2.5) *
                      pow((1.0 - ss), 0.5) +
                  78125.0 / V * sint / pow(par1, 3.5) * pow((1.0 - ss), 0.5);

            // the matrix is singular when theta = pi/2
            if (abs(y[i]) < toll && abs(cos(theta)) < toll)
            {
                J22 = -39062.5 / pow(par1, 3.5) / V +
                      3125 / pow(par1, 2.5) / (VV * V) +
                      12.5 / pow(par1, 1.5) * V + 312.5 / pow(par1, 2.5) * V +
                      7812.5 / pow(par1, 3.5) * V -
                      0.25 *
                          (-1.0 / pow(par1, 0.5) * V /
                               (1.0 - 0.2 * pow(par1, 0.5)) -
                           (1.0 + 0.2 * pow(par1, 0.5)) /
                               pow((1.0 - 0.2 * pow(par1, 0.5)), 2.0) /
                               pow(par1, 0.5) * V) /
                          (1.0 + 0.2 * pow(par1, 0.5)) *
                          (1.0 - 0.2 * pow(par1, 0.5));

                // dV = -dV/dx * Fx
                dV     = -1.0 / J22 * Fx;
                dtheta = 0.0;
                theta  = M_PI / 2.0;
            }
            else
            {
                J22 = 3125.0 / VV * cos(theta) / pow(par1, 2.5) *
                          pow((1.0 - ss), 0.5) -
                      3125.0 / VV * ss / pow(par1, 2.5) / pow((1.0 - ss), 0.5) *
                          cos(theta);

                det = -1.0 / (J11 * J22 - J12 * J21);

                // [dV dtheta]' = -[invJ]*[Fx Fy]'
                dV     = det * (J22 * Fx - J12 * Fy);
                dtheta = det * (-J21 * Fx + J11 * Fy);
            }

            V     = V + dV;
            theta = theta + dtheta;

            errV     = abs(dV);
            errTheta = abs(dtheta);
        }

        c = sqrt(1.0 - gamma_1_2 * VV);
        r = pow(c, 1.0 / gamma_1_2);

        rho[i]  = r;
        rhou[i] = rho[i] * V * cos(theta);
        rhov[i] = rho[i] * V * sin(theta);
        P       = (c * c) * rho[i] / gamma;
        E[i]    = P / (gamma - 1.0) +
               0.5 * (rhou[i] * rhou[i] / rho[i] + rhov[i] * rhov[i] / rho[i]);

        // Resetting the guess value
        errV     = 1.0;
        errTheta = 1.0;
        theta    = M_PI / 4.0;
        V        = kExt * sin(theta);
    }

    switch (field)
    {
        case 0:
            outarray = rho;
            break;
        case 1:
            outarray = rhou;
            break;
        case 2:
            outarray = rhov;
            break;
        case 3:
            outarray = E;
            break;
        default:
            ASSERTL0(false, "Error in variable number!");
            break;
    }
}

void CompressibleFlowSystem::v_ALEInitObject(
    int spaceDim, Array<OneD, MultiRegions::ExpListSharedPtr> &fields)
{
    m_spaceDim  = spaceDim;
    m_fieldsALE = fields;

    // Initialise grid velocities as 0s
    m_gridVelocity      = Array<OneD, Array<OneD, NekDouble>>(m_spaceDim);
    m_gridVelocityTrace = Array<OneD, Array<OneD, NekDouble>>(m_spaceDim);
    for (int i = 0; i < spaceDim; ++i)
    {
        m_gridVelocity[i] =
            Array<OneD, NekDouble>(fields[0]->GetTotPoints(), 0.0);
        m_gridVelocityTrace[i] =
            Array<OneD, NekDouble>(fields[0]->GetTrace()->GetTotPoints(), 0.0);
    }

    ALEHelper::InitObject(spaceDim, fields);
}

} // namespace Nektar
