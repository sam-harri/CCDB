///////////////////////////////////////////////////////////////////////////////
//
// File: UnsteadySystem.cpp
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
// Description: Generic timestepping for Unsteady solvers
//
///////////////////////////////////////////////////////////////////////////////

#include <iomanip>
#include <iostream>
using namespace std;

#include <boost/format.hpp>

#include <LibUtilities/BasicUtils/Filesystem.hpp>
#include <LibUtilities/BasicUtils/Timer.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <SolverUtils/UnsteadySystem.h>

namespace Nektar::SolverUtils
{
std::string UnsteadySystem::cmdSetStartTime =
    LibUtilities::SessionReader::RegisterCmdLineArgument(
        "set-start-time", "", "Set the starting time of the simulation.");

std::string UnsteadySystem::cmdSetStartChkNum =
    LibUtilities::SessionReader::RegisterCmdLineArgument(
        "set-start-chknumber", "",
        "Set the starting number of the checkpoint file.");

/**
 * @class UnsteadySystem
 *
 * Provides the underlying timestepping framework for unsteady solvers
 * including the general timestepping routines. This class is not
 * intended to be directly instantiated, but rather is a base class
 * on which to define unsteady solvers.
 *
 * For details on implementing unsteady solvers see
 * \ref sectionADRSolverModuleImplementation here
 */

/**
 * Processes SolverInfo parameters from the session file and sets up
 * timestepping-specific code.
 * @param   pSession        Session object to read parameters from.
 */
UnsteadySystem::UnsteadySystem(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : EquationSystem(pSession, pGraph), SolverUtils::ALEHelper()

{
}

/**
 * Initialization object for UnsteadySystem class.
 */
void UnsteadySystem::v_InitObject(bool DeclareField)
{
    EquationSystem::v_InitObject(DeclareField);
    v_ALEInitObject(m_spacedim, m_fields);
    m_initialStep = 0;

    // Load SolverInfo parameters.
    m_session->MatchSolverInfo("DIFFUSIONADVANCEMENT", "Explicit",
                               m_explicitDiffusion, true);
    m_session->MatchSolverInfo("ADVECTIONADVANCEMENT", "Explicit",
                               m_explicitAdvection, true);
    m_session->MatchSolverInfo("REACTIONADVANCEMENT", "Explicit",
                               m_explicitReaction, true);
    m_session->LoadParameter("CheckAbortSteps", m_abortSteps, 1);

    // Steady state tolerance.
    m_session->LoadParameter("SteadyStateTol", m_steadyStateTol, 0.0);

    // Frequency for checking steady state.
    m_session->LoadParameter("SteadyStateSteps", m_steadyStateSteps, 1);

    // For steady problems, we do not initialise the time integration.
    if (m_session->DefinesSolverInfo("TimeIntegrationMethod") ||
        m_session->DefinesTimeIntScheme())
    {
        LibUtilities::TimeIntScheme timeInt;
        if (m_session->DefinesTimeIntScheme())
        {
            timeInt = m_session->GetTimeIntScheme();
        }
        else
        {
            timeInt.method = m_session->GetSolverInfo("TimeIntegrationMethod");
        }

        m_intScheme =
            LibUtilities::GetTimeIntegrationSchemeFactory().CreateInstance(
                timeInt.method, timeInt.variant, timeInt.order,
                timeInt.freeParams);

        // Load generic input parameters.
        m_session->LoadParameter("IO_InfoSteps", m_infosteps, 0);
        m_session->LoadParameter("IO_FiltersInfoSteps", m_filtersInfosteps,
                                 10 * m_infosteps);
        m_session->LoadParameter("CFL", m_cflSafetyFactor, 0.0);
        m_session->LoadParameter("CFLEnd", m_CFLEnd, 0.0);
        m_session->LoadParameter("CFLGrowth", m_CFLGrowth, 1.0);

        // Ensure that there is no conflict of parameters.
        if (m_cflSafetyFactor > 0.0)
        {
            // Check final condition.
            ASSERTL0(m_fintime == 0.0 || m_steps == 0,
                     "Final condition not unique: "
                     "fintime > 0.0 and Nsteps > 0");
            // Check timestep condition.
            ASSERTL0(m_timestep == 0.0,
                     "Timestep not unique: timestep > 0.0 & CFL > 0.0");
        }
        else
        {
            ASSERTL0(m_timestep != 0.0, "Need to set either TimeStep or CFL");
        }

        // Ensure that there is no conflict of parameters.
        if (m_CFLGrowth > 1.0)
        {
            // Check final condition.
            ASSERTL0(m_CFLEnd >= m_cflSafetyFactor,
                     "m_CFLEnd >= m_cflSafetyFactor required");
        }

        // Set up time to be dumped in field information.
        m_fieldMetaDataMap["Time"] = boost::lexical_cast<std::string>(m_time);
    }

    // By default attempt to forward transform initial condition.
    m_homoInitialFwd = true;

    // Set up filters.
    for (auto &x : m_session->GetFilters())
    {
        m_filters.push_back(make_pair(
            x.first, GetFilterFactory().CreateInstance(
                         x.first, m_session, shared_from_this(), x.second)));
    }
}

/**
 * Destructor for the class UnsteadyAdvection.
 */
UnsteadySystem::~UnsteadySystem()
{
}

/**
 * @brief Returns the maximum time estimator for CFL control.
 */
NekDouble UnsteadySystem::MaxTimeStepEstimator()
{
    return m_intScheme->GetTimeStability();
}

/**
 * @brief Returns the time integration scheme.
 */
LibUtilities::TimeIntegrationSchemeSharedPtr &UnsteadySystem::
    GetTimeIntegrationScheme()
{
    return m_intScheme;
}

/**
 * @brief Returns the time integration scheme operators.
 */
LibUtilities::TimeIntegrationSchemeOperators &UnsteadySystem::
    GetTimeIntegrationSchemeOperators()
{
    return m_ode;
}

/**
 * @brief Initialises the time integration scheme (as specified in the
 * session file), and perform the time integration.
 */
void UnsteadySystem::v_DoSolve()
{
    ASSERTL0(m_intScheme != nullptr, "No time integration scheme.");

    int i          = 1;
    int nvariables = 0;
    int nfields    = m_fields.size();
    if (m_intVariables.empty())
    {
        for (i = 0; i < nfields; ++i)
        {
            m_intVariables.push_back(i);
        }
        nvariables = nfields;
    }
    else
    {
        nvariables = m_intVariables.size();
    }

    // Integrate in wave-space if using homogeneous1D.
    if (m_HomogeneousType != eNotHomogeneous && m_homoInitialFwd)
    {
        for (i = 0; i < nfields; ++i)
        {
            m_fields[i]->HomogeneousFwdTrans(m_fields[i]->GetTotPoints(),
                                             m_fields[i]->GetPhys(),
                                             m_fields[i]->UpdatePhys());
            m_fields[i]->SetWaveSpace(true);
            m_fields[i]->SetPhysState(false);
        }
    }

    // Set up wrapper to fields data storage.
    Array<OneD, Array<OneD, NekDouble>> fields(nvariables);

    // Order storage to list time-integrated fields first.
    // @TODO: Refactor to take coeffs (FwdTrans) if boolean flag (in constructor
    // function) says to.
    for (i = 0; i < nvariables; ++i)
    {
        fields[i] = m_fields[m_intVariables[i]]->UpdatePhys();
        m_fields[m_intVariables[i]]->SetPhysState(false);
    }

    // @TODO: Virtual function that allows to transform the field space, embed
    // the MultiplyMassMatrix in here.
    // @TODO: Specify what the fields variables are physical or coefficient,
    // boolean in UnsteadySystem class...

    v_ALEPreMultiplyMass(fields);

    // Initialise time integration scheme.
    m_intScheme->InitializeScheme(m_timestep, fields, m_time, m_ode);

    // Initialise filters.
    for (auto &x : m_filters)
    {
        x.second->Initialise(m_fields, m_time);
    }

    LibUtilities::Timer timer;
    bool doCheckTime        = false;
    int step                = m_initialStep;
    int stepCounter         = 0;
    NekDouble intTime       = 0.0;
    NekDouble cpuTime       = 0.0;
    NekDouble cpuPrevious   = 0.0;
    NekDouble elapsed       = 0.0;
    NekDouble totFilterTime = 0.0;

    m_lastCheckTime = 0.0;

    Array<OneD, int> abortFlags(2, 0);
    string abortFile = "abort";
    if (m_session->DefinesSolverInfo("CheckAbortFile"))
    {
        abortFile = m_session->GetSolverInfo("CheckAbortFile");
    }

    NekDouble tmp_cflSafetyFactor = m_cflSafetyFactor;
    while ((step < m_steps || m_time < m_fintime - NekConstants::kNekZeroTol) &&
           abortFlags[1] == 0)
    {
        if (m_CFLGrowth > 1.0 && m_cflSafetyFactor < m_CFLEnd)
        {
            tmp_cflSafetyFactor =
                min(m_CFLEnd, m_CFLGrowth * tmp_cflSafetyFactor);
        }

        // Frozen preconditioner checks.
        if (!m_comm->IsParallelInTime())
        {
            if (v_UpdateTimeStepCheck())
            {
                m_cflSafetyFactor = tmp_cflSafetyFactor;

                if (m_cflSafetyFactor)
                {
                    m_timestep = GetTimeStep(fields);
                }

                // Ensure that the final timestep finishes at the final
                // time, or at a prescribed IO_CheckTime.
                if (m_time + m_timestep > m_fintime && m_fintime > 0.0)
                {
                    m_timestep = m_fintime - m_time;
                }
                else if (m_checktime &&
                         m_time + m_timestep - m_lastCheckTime >= m_checktime)
                {
                    m_lastCheckTime += m_checktime;
                    m_timestep  = m_lastCheckTime - m_time;
                    doCheckTime = true;
                }
            }
        }

        if (m_TimeIncrementFactor > 1.0)
        {
            NekDouble timeincrementFactor = m_TimeIncrementFactor;
            m_timestep *= timeincrementFactor;

            if (m_time + m_timestep > m_fintime && m_fintime > 0.0)
            {
                m_timestep = m_fintime - m_time;
            }
        }

        // Perform any solver-specific pre-integration steps.
        timer.Start();
        if (v_PreIntegrate(
                step)) // Could be possible to put a preintegrate step in the
                       // ALEHelper class, put in the Unsteady Advection class
        {
            break;
        }

        ASSERTL0(m_timestep > 0, "m_timestep < 0");

        fields = m_intScheme->TimeIntegrate(stepCounter, m_timestep);
        timer.Stop();

        m_time += m_timestep;
        elapsed = timer.TimePerTest(1);
        intTime += elapsed;
        cpuTime += elapsed;

        // Write out status information.
        v_PrintStatusInformation(step, cpuTime);
        if (m_infosteps &&
            m_session->GetComm()->GetSpaceComm()->GetRank() == 0 &&
            !((step + 1) % m_infosteps))
        {
            cpuPrevious = cpuTime;
            cpuTime     = 0.0;
        }

        // @TODO: Another virtual function with this in it based on if ALE or
        // not.
        if (m_ALESolver) // Change to advect coeffs, change flag to physical vs
                         // coefficent space
        {
            SetBoundaryConditions(m_time);
            ALEHelper::ALEDoElmtInvMass(m_traceNormals, fields, m_time);
        }
        else
        {
            // Transform data into coefficient space
            for (i = 0; i < nvariables; ++i)
            {
                // copy fields into ExpList::m_phys and assign the new
                // array to fields
                m_fields[m_intVariables[i]]->SetPhys(fields[i]);
                fields[i] = m_fields[m_intVariables[i]]->UpdatePhys();
                if (v_RequireFwdTrans())
                {
                    if (m_comm->IsParallelInTime())
                    {
                        m_fields[m_intVariables[i]]->FwdTrans(
                            m_fields[m_intVariables[i]]->GetPhys(),
                            m_fields[m_intVariables[i]]->UpdateCoeffs());
                    }
                    else
                    {
                        m_fields[m_intVariables[i]]->FwdTransLocalElmt(
                            m_fields[m_intVariables[i]]->GetPhys(),
                            m_fields[m_intVariables[i]]->UpdateCoeffs());
                    }
                }
                m_fields[m_intVariables[i]]->SetPhysState(false);
            }
        }

        // Perform any solver-specific post-integration steps.
        if (v_PostIntegrate(step))
        {
            break;
        }

        // Check for steady-state.
        if (m_steadyStateSteps && m_steadyStateTol > 0.0 &&
            (!((step + 1) % m_steadyStateSteps)))
        {
            if (CheckSteadyState(step, intTime))
            {
                if (m_comm->GetRank() == 0)
                {
                    cout << "Reached Steady State to tolerance "
                         << m_steadyStateTol << endl;
                }
                break;
            }
        }

        // Test for abort conditions (nan, or abort file).
        if (m_abortSteps && !((step + 1) % m_abortSteps))
        {
            abortFlags[0] = 0;
            for (i = 0; i < nvariables; ++i)
            {
                if (Vmath::Nnan(m_fields[m_intVariables[i]]->GetPhys().size(),
                                m_fields[m_intVariables[i]]->GetPhys(), 1) > 0)
                {
                    abortFlags[0] = 1;
                }
            }

            // Rank zero looks for abort file and deletes it
            // if it exists. Communicates the abort.
            if (m_session->GetComm()->GetSpaceComm()->GetRank() == 0)
            {
                if (fs::exists(abortFile))
                {
                    fs::remove(abortFile);
                    abortFlags[1] = 1;
                }
            }

            m_session->GetComm()->GetSpaceComm()->AllReduce(
                abortFlags, LibUtilities::ReduceMax);

            ASSERTL0(!abortFlags[0], "NaN found during time integration.");
        }

        // Update filters.
        for (auto &x : m_filters)
        {
            timer.Start();
            x.second->Update(m_fields, m_time);
            timer.Stop();
            elapsed = timer.TimePerTest(1);
            totFilterTime += elapsed;

            // Write out individual filter status information.
            if (m_filtersInfosteps && m_session->GetComm()->GetRank() == 0 &&
                !((step + 1) % m_filtersInfosteps) && !m_filters.empty() &&
                m_session->DefinesCmdLineArgument("verbose"))
            {
                stringstream s0;
                s0 << x.first << ":";
                stringstream s1;
                s1 << elapsed << "s";
                stringstream s2;
                s2 << elapsed / cpuPrevious * 100 << "%";
                cout << "CPU time for filter " << setw(25) << left << s0.str()
                     << setw(12) << left << s1.str() << endl
                     << "\t Percentage of time integration:     " << setw(10)
                     << left << s2.str() << endl;
            }
        }

        // Write out overall filter status information.
        if (m_filtersInfosteps && m_session->GetComm()->GetRank() == 0 &&
            !((step + 1) % m_filtersInfosteps) && !m_filters.empty())
        {
            stringstream ss;
            ss << totFilterTime << "s";
            cout << "Total filters CPU Time:\t\t\t     " << setw(10) << left
                 << ss.str() << endl;
        }
        totFilterTime = 0.0;

        // Write out checkpoint files.
        if ((m_checksteps && !((step + 1) % m_checksteps)) || doCheckTime)
        {
            if (m_HomogeneousType != eNotHomogeneous && !m_ALESolver)
            {
                // Transform to physical space for output.
                vector<bool> transformed(nfields, false);
                for (i = 0; i < nfields; i++)
                {
                    if (m_fields[i]->GetWaveSpace())
                    {
                        m_fields[i]->SetWaveSpace(false);
                        m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(),
                                              m_fields[i]->UpdatePhys());
                        transformed[i] = true;
                    }
                }
                Checkpoint_Output(m_nchk);
                m_nchk++;

                // Transform back to wave-space after output.
                for (i = 0; i < nfields; i++)
                {
                    if (transformed[i])
                    {
                        m_fields[i]->SetWaveSpace(true);
                        m_fields[i]->HomogeneousFwdTrans(
                            m_fields[i]->GetTotPoints(), m_fields[i]->GetPhys(),
                            m_fields[i]->UpdatePhys());
                        m_fields[i]->SetPhysState(false);
                    }
                }
            }
            else
            {
                Checkpoint_Output(m_nchk);
                m_nchk++;
            }
            doCheckTime = false;
        }

        // Step advance.
        ++step;
        ++stepCounter;
    }

    // Print out summary statistics.
    v_PrintSummaryStatistics(intTime);

    // If homogeneous, transform back into physical space if necessary.
    if (!m_ALESolver)
    {
        if (m_HomogeneousType != eNotHomogeneous)
        {
            for (i = 0; i < nfields; i++)
            {
                if (m_fields[i]->GetWaveSpace())
                {
                    m_fields[i]->SetWaveSpace(false);
                    m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(),
                                          m_fields[i]->UpdatePhys());
                }
            }
        }
        else
        {
            for (i = 0; i < nvariables; ++i)
            {
                m_fields[m_intVariables[i]]->SetPhysState(true);
            }
        }
    }
    // Finalise filters.
    for (auto &x : m_filters)
    {
        x.second->Finalise(m_fields, m_time);
    }

    // Print for 1D problems.
    if (m_spacedim == 1)
    {
        AppendOutput1D();
    }
}

void UnsteadySystem::v_PrintStatusInformation(const int step,
                                              const NekDouble cpuTime)
{
    if (m_infosteps && m_session->GetComm()->GetSpaceComm()->GetRank() == 0 &&
        !((step + 1) % m_infosteps))
    {
        if (m_comm->IsParallelInTime())
        {
            cout << "RANK " << m_session->GetComm()->GetTimeComm()->GetRank()
                 << " Steps: " << setw(8) << left << step + 1 << " "
                 << "Time: " << setw(12) << left << m_time;
        }
        else
        {
            cout << "Steps: " << setw(8) << left << step + 1 << " "
                 << "Time: " << setw(12) << left << m_time;
        }

        if (m_cflSafetyFactor)
        {
            cout << " Time-step: " << setw(12) << left << m_timestep;
        }

        stringstream ss;
        ss << cpuTime << "s";
        cout << " CPU Time: " << setw(8) << left << ss.str() << endl;
    }
}

void UnsteadySystem::v_PrintSummaryStatistics(const NekDouble intTime)
{
    if (m_session->GetComm()->GetRank() == 0)
    {
        if (m_cflSafetyFactor > 0.0)
        {
            cout << "CFL safety factor : " << m_cflSafetyFactor << endl
                 << "CFL time-step     : " << m_timestep << endl;
        }

        if (m_session->GetSolverInfo("Driver") != "SteadyState" &&
            m_session->GetSolverInfo("Driver") != "Parareal" &&
            m_session->GetSolverInfo("Driver") != "PFASST")
        {
            cout << "Time-integration  : " << intTime << "s" << endl;
        }
    }
}

/**
 * @brief Sets the initial conditions.
 */
void UnsteadySystem::v_DoInitialise(bool dumpInitialConditions)
{
    CheckForRestartTime(m_time, m_nchk);
    SetBoundaryConditions(m_time);
    SetInitialConditions(m_time, dumpInitialConditions);

    v_UpdateGridVelocity(m_time);

    InitializeSteadyState();
}

/**
 * @brief Prints a summary with some information regards the
 * time-stepping.
 */
void UnsteadySystem::v_GenerateSummary(SummaryList &s)
{
    EquationSystem::v_GenerateSummary(s);
    AddSummaryItem(s, "Advect. advancement",
                   m_explicitAdvection ? "explicit" : "implicit");

    AddSummaryItem(s, "Diffuse. advancement",
                   m_explicitDiffusion ? "explicit" : "implicit");

    if (m_session->GetSolverInfo("EQTYPE") ==
        "SteadyAdvectionDiffusionReaction")
    {
        AddSummaryItem(s, "React. advancement",
                       m_explicitReaction ? "explicit" : "implicit");
    }

    AddSummaryItem(s, "Time Step", m_timestep);
    AddSummaryItem(s, "No. of Steps", m_steps);
    AddSummaryItem(s, "Checkpoints (steps)", m_checksteps);
    if (m_intScheme)
    {
        AddSummaryItem(s, "Integration Type", m_intScheme->GetName());
    }
}

/**
 * Stores the solution in a file for 1D problems only. This method has
 * been implemented to facilitate the post-processing for 1D problems.
 */
void UnsteadySystem::AppendOutput1D(void)
{
    // Coordinates of the quadrature points in the real physical space.
    Array<OneD, NekDouble> x(GetNpoints());
    Array<OneD, NekDouble> y(GetNpoints());
    Array<OneD, NekDouble> z(GetNpoints());
    m_fields[0]->GetCoords(x, y, z);

    // Print out the solution in a txt file.
    ofstream outfile;
    outfile.open("solution1D.txt");
    for (int i = 0; i < GetNpoints(); i++)
    {
        outfile << scientific << setw(17) << setprecision(16) << x[i] << "  "
                << m_fields[m_intVariables[0]]->GetPhys()[i] << endl;
    }
    outfile << endl << endl;
    outfile.close();
}

/**
 *
 */
void UnsteadySystem::CheckForRestartTime(NekDouble &time, int &nchk)
{
    if (m_session->DefinesFunction("InitialConditions"))
    {
        for (int i = 0; i < m_fields.size(); ++i)
        {
            LibUtilities::FunctionType vType;

            vType = m_session->GetFunctionType("InitialConditions",
                                               m_session->GetVariable(i));

            if (vType == LibUtilities::eFunctionTypeFile)
            {
                std::string filename = m_session->GetFunctionFilename(
                    "InitialConditions", m_session->GetVariable(i));

                fs::path pfilename(filename);

                // Redefine path for parallel file which is in directory.
                if (fs::is_directory(pfilename))
                {
                    fs::path metafile("Info.xml");
                    fs::path fullpath = pfilename / metafile;
                    filename          = LibUtilities::PortablePath(fullpath);
                }
                LibUtilities::FieldIOSharedPtr fld =
                    LibUtilities::FieldIO::CreateForFile(m_session, filename);
                fld->ImportFieldMetaData(filename, m_fieldMetaDataMap);

                // Check to see if time defined.
                if (m_fieldMetaDataMap != LibUtilities::NullFieldMetaDataMap)
                {
                    auto iter = m_fieldMetaDataMap.find("Time");
                    if (iter != m_fieldMetaDataMap.end())
                    {
                        time = std::stod(iter->second);
                    }

                    iter = m_fieldMetaDataMap.find("ChkFileNum");
                    if (iter != m_fieldMetaDataMap.end())
                    {
                        nchk = std::stod(iter->second);
                    }
                }

                break;
            }
        }
    }
    if (m_session->DefinesCmdLineArgument("set-start-time"))
    {
        time = std::stod(
            m_session->GetCmdLineArgument<std::string>("set-start-time")
                .c_str());
    }
    if (m_session->DefinesCmdLineArgument("set-start-chknumber"))
    {
        nchk = boost::lexical_cast<int>(
            m_session->GetCmdLineArgument<std::string>("set-start-chknumber"));
    }
    ASSERTL0(time >= 0 && nchk >= 0,
             "Starting time and checkpoint number should be >= 0");
}

/**
 * @brief Return the timestep to be used for the next step in the
 * time-marching loop.
 *
 * @see UnsteadySystem::GetTimeStep
 */
NekDouble UnsteadySystem::v_GetTimeStep(
    [[maybe_unused]] const Array<OneD, const Array<OneD, NekDouble>> &inarray)
{
    NEKERROR(ErrorUtil::efatal, "Not defined for this class");
    return 0.0;
}

/**
 *
 */
bool UnsteadySystem::v_PreIntegrate([[maybe_unused]] int step)
{
    return false;
}

/**
 *
 */
bool UnsteadySystem::v_PostIntegrate([[maybe_unused]] int step)
{
    return false;
}

/**
 *
 */
bool UnsteadySystem::v_RequireFwdTrans(void)
{
    return true;
}

/**
 *
 */
bool UnsteadySystem::v_UpdateTimeStepCheck(void)
{
    return true;
}

/**
 *
 */
void UnsteadySystem::SVVVarDiffCoeff(
    const Array<OneD, Array<OneD, NekDouble>> vel,
    StdRegions::VarCoeffMap &varCoeffMap)
{
    int phystot = m_fields[0]->GetTotPoints();
    int nvel    = vel.size();

    Array<OneD, NekDouble> varcoeff(phystot), tmp;

    // Calculate magnitude of v.
    Vmath::Vmul(phystot, vel[0], 1, vel[0], 1, varcoeff, 1);
    for (int n = 1; n < nvel; ++n)
    {
        Vmath::Vvtvp(phystot, vel[n], 1, vel[n], 1, varcoeff, 1, varcoeff, 1);
    }
    Vmath::Vsqrt(phystot, varcoeff, 1, varcoeff, 1);

    for (int i = 0; i < m_fields[0]->GetNumElmts(); ++i)
    {
        int offset = m_fields[0]->GetPhys_Offset(i);
        int nq     = m_fields[0]->GetExp(i)->GetTotPoints();
        Array<OneD, NekDouble> unit(nq, 1.0);

        int nmodes = 0;

        for (int n = 0; n < m_fields[0]->GetExp(i)->GetNumBases(); ++n)
        {
            nmodes = max(nmodes, m_fields[0]->GetExp(i)->GetBasisNumModes(n));
        }

        NekDouble h = m_fields[0]->GetExp(i)->Integral(unit);
        h           = pow(h, (NekDouble)(1.0 / nvel)) / ((NekDouble)nmodes);

        Vmath::Smul(nq, h, varcoeff + offset, 1, tmp = varcoeff + offset, 1);
    }

    // Set up map with eVarCoffLaplacian key.
    varCoeffMap[StdRegions::eVarCoeffLaplacian] = varcoeff;
}

/**
 *
 */
void UnsteadySystem::InitializeSteadyState()
{
    if (m_steadyStateTol > 0.0)
    {
        const int nPoints = m_fields[0]->GetTotPoints();
        m_previousSolution =
            Array<OneD, Array<OneD, NekDouble>>(m_fields.size());

        for (int i = 0; i < m_fields.size(); ++i)
        {
            m_previousSolution[i] = Array<OneD, NekDouble>(nPoints);
            Vmath::Vcopy(nPoints, m_fields[i]->GetPhys(), 1,
                         m_previousSolution[i], 1);
        }

        if (m_comm->GetRank() == 0)
        {
            std::string fName =
                m_session->GetSessionName() + std::string(".resd");
            m_errFile.open(fName.c_str());
            m_errFile << setw(26) << left << "# Time";

            m_errFile << setw(26) << left << "CPU_Time";

            m_errFile << setw(26) << left << "Step";

            for (int i = 0; i < m_fields.size(); ++i)
            {
                m_errFile << setw(26) << m_session->GetVariables()[i];
            }

            m_errFile << endl;
        }
    }
}

/**
 * @brief Calculate whether the system has reached a steady state by
 * observing residuals to a user-defined tolerance.
 */
bool UnsteadySystem::CheckSteadyState(int step, const NekDouble &totCPUTime)
{
    const int nPoints = GetTotPoints();
    const int nFields = m_fields.size();

    // Holds L2 errors.
    Array<OneD, NekDouble> L2(nFields);

    SteadyStateResidual(step, L2);

    if (m_infosteps && m_comm->GetRank() == 0 &&
        (((step + 1) % m_infosteps == 0) || ((step == m_initialStep))))
    {
        // Output time.
        m_errFile << boost::format("%25.19e") % m_time;

        m_errFile << " " << boost::format("%25.19e") % totCPUTime;

        int stepp = step + 1;

        m_errFile << " " << boost::format("%25.19e") % stepp;

        // Output residuals.
        for (int i = 0; i < nFields; ++i)
        {
            m_errFile << " " << boost::format("%25.19e") % L2[i];
        }

        m_errFile << endl;
    }

    // Calculate maximum L2 error.
    NekDouble maxL2 = Vmath::Vmax(nFields, L2, 1);

    if (m_infosteps && m_session->DefinesCmdLineArgument("verbose") &&
        m_comm->GetRank() == 0 && ((step + 1) % m_infosteps == 0))
    {
        cout << "-- Maximum L^2 residual: " << maxL2 << endl;
    }

    if (maxL2 <= m_steadyStateTol)
    {
        return true;
    }

    for (int i = 0; i < m_fields.size(); ++i)
    {
        Vmath::Vcopy(nPoints, m_fields[i]->GetPhys(), 1, m_previousSolution[i],
                     1);
    }

    return false;
}

/**
 *
 */
void UnsteadySystem::v_SteadyStateResidual([[maybe_unused]] int step,
                                           Array<OneD, NekDouble> &L2)
{
    const int nPoints = GetTotPoints();
    const int nFields = m_fields.size();

    // Holds L2 errors.
    Array<OneD, NekDouble> RHSL2(nFields);
    Array<OneD, NekDouble> residual(nFields);
    Array<OneD, NekDouble> reference(nFields);

    for (int i = 0; i < nFields; ++i)
    {
        Array<OneD, NekDouble> tmp(nPoints);

        Vmath::Vsub(nPoints, m_fields[i]->GetPhys(), 1, m_previousSolution[i],
                    1, tmp, 1);
        Vmath::Vmul(nPoints, tmp, 1, tmp, 1, tmp, 1);
        residual[i] = Vmath::Vsum(nPoints, tmp, 1);

        Vmath::Vmul(nPoints, m_previousSolution[i], 1, m_previousSolution[i], 1,
                    tmp, 1);
        reference[i] = Vmath::Vsum(nPoints, tmp, 1);
    }

    m_comm->GetSpaceComm()->AllReduce(residual, LibUtilities::ReduceSum);
    m_comm->GetSpaceComm()->AllReduce(reference, LibUtilities::ReduceSum);

    // L2 error.
    for (int i = 0; i < nFields; ++i)
    {
        reference[i] = (reference[i] == 0) ? 1 : reference[i];
        L2[i]        = sqrt(residual[i] / reference[i]);
    }
}

/**
 *
 */
void UnsteadySystem::DoDummyProjection(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray,
    [[maybe_unused]] const NekDouble time)
{

    if (&inarray != &outarray)
    {
        for (int i = 0; i < inarray.size(); ++i)
        {
            Vmath::Vcopy(GetNpoints(), inarray[i], 1, outarray[i], 1);
        }
    }
}

} // namespace Nektar::SolverUtils
