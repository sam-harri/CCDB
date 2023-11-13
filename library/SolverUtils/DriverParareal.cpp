///////////////////////////////////////////////////////////////////////////////
//
// File DriverParareal.cpp
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
// Description: Driver class for the parareal solver
//
///////////////////////////////////////////////////////////////////////////////

#include <iomanip>

#include <LibUtilities/BasicUtils/Timer.h>
#include <SolverUtils/DriverParareal.h>
#include <boost/format.hpp>

namespace Nektar::SolverUtils
{
std::string DriverParareal::className =
    GetDriverFactory().RegisterCreatorFunction("Parareal",
                                               DriverParareal::create);
std::string DriverParareal::driverLookupId =
    LibUtilities::SessionReader::RegisterEnumValue("Driver", "Parareal", 0);

/**
 *
 */
DriverParareal::DriverParareal(
    const LibUtilities::SessionReaderSharedPtr pSession,
    const SpatialDomains::MeshGraphSharedPtr pGraph)
    : DriverParallelInTime(pSession, pGraph)
{
}

/**
 *
 */
void DriverParareal::v_InitObject(std::ostream &out)
{
    DriverParallelInTime::v_InitObject(out);
    DriverParallelInTime::GetParametersFromSession();
    DriverParallelInTime::PrintSolverInfo(out);
    DriverParallelInTime::InitialiseEqSystem(false);
    AllocateMemory();
}

/**
 *
 */
void DriverParareal::v_Execute(std::ostream &out)
{
    boost::ignore_unused(out);

    // Get and assert parameters from session file.
    AssertParameters();

    // Print solver info.
    m_totalTime = m_timestep[m_fineLevel] * m_nsteps[m_fineLevel];
    m_chunkTime = m_totalTime / m_numChunks / m_numWindowsPIT;
    m_nsteps[m_fineLevel] /= m_numChunks * m_numWindowsPIT;
    m_nsteps[m_coarseLevel] /= m_numChunks * m_numWindowsPIT;

    // Start iteration windows.
    m_comm->GetTimeComm()->Block();
    m_CPUtime = 0.0;
    m_timer.Start();
    UpdateInitialConditionFromSolver(m_fineLevel);
    for (size_t w = 0; w < m_numWindowsPIT; w++)
    {
        // Initialize time for the current window.
        m_time = (w * m_numChunks) * m_chunkTime;

        // Print window number.
        PrintHeader((boost::format("WINDOWS #%1%") % (w + 1)).str(), '*');

        // Run predictor
        UpdateSolverInitialCondition(m_coarseLevel);
        if (m_chunkRank > 0)
        {
            // Calculate coarse solution
            UpdateSolution(m_coarseLevel, m_time,
                           m_chunkRank * m_nsteps[m_coarseLevel]);
            UpdateInitialConditionFromSolver(m_coarseLevel);
            m_time += m_chunkRank * m_chunkTime;
        }
        UpdateSolution(m_coarseLevel, m_time, m_nsteps[m_coarseLevel]);

        // Interpolate coarse solution.
        InterpolateCoarseSolution();

        // Compute exact solution, if necessary.
        if (m_exactSolution)
        {
            EvaluateExactSolution(m_fineLevel, m_time + m_chunkTime);
        }

        // Solution convergence monitoring.
        CopyToPhysField(m_fineLevel, m_coarseSolution);
        SolutionConvergenceMonitoring(m_fineLevel, 0);

        // Start Parareal iteration.
        size_t iter         = 1;
        int convergenceCurr = false;
        int convergencePrev = (m_chunkRank == 0);
        while (iter <= m_iterMaxPIT && !convergenceCurr)
        {
            // Use previous parareal solution as "exact solution", if necessary.
            if (!m_exactSolution)
            {
                CopyFromPhysField(m_fineLevel, m_exactsoln);
            }

            // Calculate fine solution (parallel-in-time).
            UpdateSolverInitialCondition(m_fineLevel);
            UpdateSolution(m_fineLevel, m_time, m_nsteps[m_fineLevel], w, iter);

            // Compute F -> F - Gold
            CorrectionWithOldCoarseSolution();

            // Receive coarse solution from previous processor.
            RecvFromPreviousProc(m_initialCondition, convergencePrev);

            // Calculate coarse solution (serial-in-time).
            UpdateSolverInitialCondition(m_coarseLevel);
            UpdateSolution(m_coarseLevel, m_time, m_nsteps[m_coarseLevel]);

            // Compute F -> F + Gnew
            CorrectionWithNewCoarseSolution();

            // Solution convergence monitoring.
            SolutionConvergenceMonitoring(m_fineLevel, iter);

            // Check convergence of L2 error for each time chunk.
            convergenceCurr = (vL2ErrorMax() < m_tolerPIT && convergencePrev) ||
                              (m_chunkRank + 1 == iter);

            // Send solution to next processor.
            SendToNextProc(m_fineSolution, convergenceCurr);

            // Increment iteration index.
            iter++;
        }

        // Copy converged check points.
        CopyConvergedCheckPoints(w, iter);

        // Write time chunk solution to files.
        WriteTimeChunkOuput();

        // Apply windowing.
        ApplyWindowing(w);
    }
    m_timer.Stop();

    m_comm->GetTimeComm()->Block();
    PrintHeader("SUMMARY", '*');
    EvaluateExactSolution(m_fineLevel, m_time + m_chunkTime);
    SolutionConvergenceSummary(m_fineLevel);
    SpeedUpAnalysis();
}

/**
 *
 */
void DriverParareal::AllocateMemory(void)
{
    // Allocate storage for Parareal solver.
    m_initialCondition = Array<OneD, Array<OneD, NekDouble>>(m_nVar);
    m_fineSolution     = Array<OneD, Array<OneD, NekDouble>>(m_nVar);
    m_coarseSolution   = Array<OneD, Array<OneD, NekDouble>>(m_nVar);
    for (size_t i = 0; i < m_nVar; ++i)
    {
        m_initialCondition[i] =
            Array<OneD, NekDouble>(m_npts[m_fineLevel], 0.0);
        m_fineSolution[i] = m_EqSys[m_fineLevel]->UpdatePhysField(i);
        m_coarseSolution[i] =
            (m_npts[m_fineLevel] == m_npts[m_coarseLevel])
                ? m_EqSys[m_coarseLevel]->UpdatePhysField(i)
                : Array<OneD, NekDouble>(m_npts[m_fineLevel], 0.0);
    }
}

/**
 *
 */
void DriverParareal::AssertParameters(void)
{
    // Assert time-stepping parameters
    ASSERTL0(
        m_nsteps[m_fineLevel] % m_numChunks == 0,
        "Total number of fine step should be divisible by number of chunks.");

    ASSERTL0(
        m_nsteps[m_coarseLevel] % m_numChunks == 0,
        "Total number of coarse step should be divisible by number of chunks.");

    ASSERTL0(m_nsteps[m_fineLevel] % (m_numChunks * m_numWindowsPIT) == 0,
             "Total number of fine step should be divisible by number of "
             "windows times number of chunks.");

    ASSERTL0(m_nsteps[m_coarseLevel] % (m_numChunks * m_numWindowsPIT) == 0,
             "Total number of coarse step should be divisible by number of "
             "windows times number of chunks.");

    ASSERTL0(fabs(m_timestep[m_coarseLevel] * m_nsteps[m_coarseLevel] -
                  m_timestep[m_fineLevel] * m_nsteps[m_fineLevel]) < 1e-12,
             "Fine and coarse total computational times do not match");

    ASSERTL0(m_EqSys[m_fineLevel]
                     ->GetTimeIntegrationScheme()
                     ->GetNumIntegrationPhases() == 1,
             "Only single step time-integration schemes currently supported "
             "for Parareal");

    ASSERTL0(m_EqSys[m_coarseLevel]
                     ->GetTimeIntegrationScheme()
                     ->GetNumIntegrationPhases() == 1,
             "Only single step time-integration schemes currently supported "
             "for Parareal");

    // Assert I/O parameters
    if (m_infoSteps)
    {
        ASSERTL0(m_nsteps[m_fineLevel] %
                         (m_infoSteps * m_numChunks * m_numWindowsPIT) ==
                     0,
                 "number of IO_InfoSteps should divide number of fine steps "
                 "per time chunk");
    }

    if (m_checkSteps)
    {
        ASSERTL0(m_nsteps[m_fineLevel] %
                         (m_checkSteps * m_numChunks * m_numWindowsPIT) ==
                     0,
                 "number of IO_CheckSteps should divide number of fine steps "
                 "per time chunk");
    }
}

/**
 *
 */
void DriverParareal::UpdateInitialConditionFromSolver(const size_t timeLevel)
{
    // Interpolate solution to fine field.
    InterpExp1ToExp2(m_EqSys[timeLevel]->UpdateFields(),
                     m_EqSys[m_fineLevel]->UpdateFields(),
                     NullNekDoubleArrayOfArray, m_initialCondition);
}

/**
 *
 */
void DriverParareal::UpdateSolverInitialCondition(const size_t timeLevel)
{
    // Restrict fine field to coarse solution.
    InterpExp1ToExp2(m_EqSys[m_fineLevel]->UpdateFields(),
                     m_EqSys[timeLevel]->UpdateFields(), m_initialCondition,
                     NullNekDoubleArrayOfArray);
}

/**
 *
 */
void DriverParareal::UpdateSolution(const size_t timeLevel,
                                    const NekDouble time, const size_t nstep,
                                    const size_t wd, const size_t iter)
{
    if (timeLevel == m_fineLevel)
    {
        // Number of checkpoint by chunk.
        size_t nChkPts =
            m_checkSteps ? m_nsteps[m_fineLevel] / m_checkSteps : 1;

        // Checkpoint index.
        size_t iChkPts = (m_chunkRank + wd * m_numChunks) * nChkPts + 1;

        // Reinitialize check point number for each parallel-in-time
        // iteration.
        m_EqSys[m_fineLevel]->SetCheckpointNumber(iChkPts);

        // Update parallel-in-time iteration number.
        m_EqSys[m_fineLevel]->SetIterationNumberPIT(iter);

        // Update parallel-in-time window number.
        m_EqSys[m_fineLevel]->SetWindowNumberPIT(wd);
    }

    m_EqSys[timeLevel]->SetTime(time);
    m_EqSys[timeLevel]->SetSteps(nstep);
    m_EqSys[timeLevel]->DoSolve();
}

/**
 *
 */
void DriverParareal::CorrectionWithOldCoarseSolution(void)
{
    // Correct solution F -> F - Gold.
    for (size_t i = 0; i < m_nVar; ++i)
    {
        Vmath::Vsub(m_npts[m_fineLevel], m_fineSolution[i], 1,
                    m_coarseSolution[i], 1, m_fineSolution[i], 1);
    }
}

/**
 *
 */
void DriverParareal::CorrectionWithNewCoarseSolution(void)
{
    // Interpolate coarse solution.
    InterpolateCoarseSolution();

    // Correct solution F -> F + Gnew.
    for (size_t i = 0; i < m_nVar; ++i)
    {
        Vmath::Vadd(m_npts[m_fineLevel], m_fineSolution[i], 1,
                    m_coarseSolution[i], 1, m_fineSolution[i], 1);
    }
}

/**
 *
 */
void DriverParareal::InterpolateCoarseSolution(void)
{
    if (m_npts[m_fineLevel] != m_npts[m_coarseLevel])
    {
        // Interpolate coarse solution to fine field.
        InterpExp1ToExp2(m_EqSys[m_coarseLevel]->UpdateFields(),
                         m_EqSys[m_fineLevel]->UpdateFields(),
                         NullNekDoubleArrayOfArray, m_coarseSolution);
    }
}

/**
 *
 */
void DriverParareal::ApplyWindowing(const size_t w)
{
    if (w == m_numWindowsPIT - 1)
    {
        // No windowing required for the last window.
        return;
    }

    // Use last chunk solution as initial condition for the next
    // window.
    if (m_chunkRank == m_numChunks - 1)
    {
        for (size_t i = 0; i < m_nVar; ++i)
        {
            Vmath::Vcopy(m_npts[m_fineLevel],
                         m_EqSys[m_fineLevel]->UpdatePhysField(i), 1,
                         m_initialCondition[i], 1);
        }
    }

    // Broadcast I.C. for windowing.
    for (size_t i = 0; i < m_nVar; ++i)
    {
        m_comm->GetTimeComm()->Bcast(m_initialCondition[i], m_numChunks - 1);
    }
}

/**
 *
 */
void DriverParareal::CopyConvergedCheckPoints(const size_t w, const size_t k)
{
    // Determine max number of iteration.
    size_t kmax = k;
    m_comm->GetTimeComm()->AllReduce(kmax, Nektar::LibUtilities::ReduceMax);

    if (m_comm->GetSpaceComm()->GetRank() == 0 && m_checkSteps)
    {
        for (size_t j = k; j < kmax; j++)
        {
            // Copy converged solution files from directory corresponding to
            // iteration j - 1 to the directory corresponding to iteration j.

            // Input directory name.
            std::string indir = m_EqSys[m_fineLevel]->GetSessionName() + "_" +
                                boost::lexical_cast<std::string>(j - 1) +
                                ".pit";

            /// Output directory name.
            std::string outdir = m_EqSys[m_fineLevel]->GetSessionName() + "_" +
                                 boost::lexical_cast<std::string>(j) + ".pit";

            // Number of checkpoint by chunk.
            size_t nChkPts = m_nsteps[m_fineLevel] / m_checkSteps;

            // Checkpoint index.
            size_t iChkPts = (m_chunkRank + w * m_numChunks) * nChkPts + 1;

            for (size_t i = 0; i < nChkPts; i++)
            {
                // Filename corresponding to checkpoint iChkPts.
                std::string filename =
                    m_EqSys[m_fineLevel]->GetSessionName() + "_" +
                    boost::lexical_cast<std::string>(iChkPts + i) + ".chk";

                // Intput full file name.
                std::string infullname = indir + "/" + filename;

                // Output full file name.
                std::string outfullname = outdir + "/" + filename;

                // Remove output file if already existing.
                fs::remove_all(outfullname);

                // Copy converged solution files.
                fs::copy(infullname, outfullname);
            }
        }
    }
}

/**
 *
 */
void DriverParareal::WriteTimeChunkOuput(void)
{
    PrintHeader("PRINT SOLUTION FILES", '-');

    // Update field coefficients.
    UpdateFieldCoeffs(m_fineLevel);

    // Output solution files.
    m_EqSys[m_fineLevel]->Output();
}

/**
 *
 */
void DriverParareal::SpeedUpAnalysis(void)
{
    // Print header.
    PrintHeader("PARAREAL SPEED-UP ANALYSIS", '*');

    // Mean communication time.
    NekDouble commTime = EstimateCommunicationTime();
    PrintHeader("Mean Communication Time = " +
                    (boost::format("%1$.6e") % commTime).str() + "s",
                '-');

    // Mean restriction time.
    NekDouble restTime = EstimateRestrictionTime();
    PrintHeader("Mean Restriction Time = " +
                    (boost::format("%1$.6e") % restTime).str() + "s",
                '-');

    // Mean interpolation time.
    NekDouble interTime = EstimateInterpolationTime();
    PrintHeader("Mean Interpolation Time = " +
                    (boost::format("%1$.6e") % interTime).str() + "s",
                '-');

    // Mean coarse solver time.
    NekDouble coarseSolveTime = EstimateSolverTime(m_coarseLevel, 10);
    PrintHeader("Mean Coarse Solve Time = " +
                    (boost::format("%1$.6e") % coarseSolveTime).str() + "s",
                '-');

    // Mean fine solver time.
    NekDouble fineSolveTime = EstimateSolverTime(m_fineLevel, 100);
    PrintHeader("Mean Fine Solve Time = " +
                    (boost::format("%1$.6e") % fineSolveTime).str() + "s",
                '-');

    // Mean predictor time.
    NekDouble predictorTime = EstimatePredictorTime();
    PrintHeader("Mean Predictor Time = " +
                    (boost::format("%1$.6e") % predictorTime).str() + "s",
                '-');

    // Print speedup time.
    PrintSpeedUp(fineSolveTime, coarseSolveTime, restTime, interTime, commTime,
                 predictorTime);
}

/**
 *
 */
void DriverParareal::PrintSpeedUp(NekDouble fineSolveTime,
                                  NekDouble coarseSolveTime, NekDouble restTime,
                                  NekDouble interTime, NekDouble commTime,
                                  NekDouble predictTime)
{
    if (m_chunkRank == m_numChunks - 1 &&
        m_comm->GetSpaceComm()->GetRank() == 0)
    {
        // Print maximum theoretical speed-up.
        PrintHeader("Maximum Speed-up", '-');
        for (size_t k = 1; k <= m_numChunks; k++)
        {
            NekDouble speedup = ComputeSpeedUp(
                k, fineSolveTime, coarseSolveTime, 0.0, 0.0, 0.0, predictTime);
            std::cout << "Speed-up (" << k << ") = " << speedup << std::endl
                      << std::flush;
        }

        // Print speed-up with interpolation and restriction.
        PrintHeader("Speed-up with interp.", '-');
        for (size_t k = 1; k <= m_numChunks; k++)
        {
            NekDouble speedup =
                ComputeSpeedUp(k, fineSolveTime, coarseSolveTime, restTime,
                               interTime, 0.0, predictTime);
            std::cout << "Speed-up (" << k << ") = " << speedup << std::endl
                      << std::flush;
        }

        // Print speed-up with interpolation, restriction and communication.
        PrintHeader("Speed-up with comm. and interp.", '-');
        for (size_t k = 1; k <= m_numChunks; k++)
        {
            NekDouble speedup =
                ComputeSpeedUp(k, fineSolveTime, coarseSolveTime, restTime,
                               interTime, commTime, predictTime);
            std::cout << "Speed-up (" << k << ") = " << speedup << std::endl
                      << std::flush;
        }
        std::cout << "-------------------------------------------" << std::endl
                  << std::flush;
    }
}

/**
 *
 */
NekDouble DriverParareal::ComputeSpeedUp(
    const size_t iter, NekDouble fineSolveTime, NekDouble coarseSolveTime,
    NekDouble restTime, NekDouble interTime, NekDouble commTime,
    NekDouble predictorTime)
{
    // The speed-up estimate is based on "Lunet, T., Bodart, J., Gratton, S., &
    // Vasseur, X. (2018). Time-parallel simulation of the decay of homogeneous
    // turbulence using parareal with spatial coarsening. Computing and
    // Visualization in Science, 19, 31-44".

    size_t nComm             = (iter * (2 * m_numChunks - iter - 1)) / 2;
    NekDouble ratio          = double(iter) / m_numChunks;
    NekDouble ratioPredictor = predictorTime / fineSolveTime;
    NekDouble ratioSolve     = coarseSolveTime / fineSolveTime;
    NekDouble ratioInterp    = (restTime + interTime) / fineSolveTime;
    NekDouble ratioComm      = commTime / fineSolveTime;

    return 1.0 / (ratioPredictor + ratio * (1.0 + ratioSolve + ratioInterp) +
                  (ratioComm * nComm) / m_numChunks);
}

/**
 *
 */
NekDouble DriverParareal::EstimateCommunicationTime(void)
{
    // Allocate memory.
    Array<OneD, Array<OneD, NekDouble>> buffer1(m_nVar);
    Array<OneD, Array<OneD, NekDouble>> buffer2(m_nVar);
    for (size_t i = 0; i < m_nVar; ++i)
    {
        buffer1[i] = Array<OneD, NekDouble>(m_npts[m_fineLevel], 0.0);
        buffer2[i] = Array<OneD, NekDouble>(m_npts[m_fineLevel], 0.0);
    }

    // Estimate communication time.
    return DriverParallelInTime::EstimateCommunicationTime(buffer1, buffer2);
}

/**
 *
 */
NekDouble DriverParareal::EstimateRestrictionTime(void)
{
    if (m_npts[m_fineLevel] == m_npts[m_coarseLevel])
    {
        return 0.0;
    }

    // Average restriction time over niter iteration.
    size_t niter = 20;
    m_timer.Start();
    for (size_t n = 0; n < niter; n++)
    {
        UpdateSolverInitialCondition(m_coarseLevel);
    }
    m_timer.Stop();
    return m_timer.Elapsed().count() / niter;
}

/**
 *
 */
NekDouble DriverParareal::EstimateInterpolationTime(void)
{
    if (m_npts[m_fineLevel] == m_npts[m_coarseLevel])
    {
        return 0.0;
    }

    // Average interpolation time over niter iteration.
    size_t niter = 20;
    m_timer.Start();
    for (size_t n = 0; n < niter; n++)
    {
        InterpolateCoarseSolution();
    }
    m_timer.Stop();
    return m_timer.Elapsed().count() / niter;
}

/**
 *
 */
NekDouble DriverParareal::EstimateSolverTime(const size_t timeLevel,
                                             const size_t nstep)
{
    // Turnoff I/O.
    m_EqSys[timeLevel]->SetInfoSteps(0);
    m_EqSys[timeLevel]->SetCheckpointSteps(0);

    // Estimate solver time.
    m_timer.Start();
    UpdateSolution(timeLevel, m_time + m_chunkTime, nstep);
    m_timer.Stop();
    return m_timer.Elapsed().count() * m_nsteps[timeLevel] / nstep;
}

/**
 *
 */
NekDouble DriverParareal::EstimatePredictorTime(void)
{
    return EstimateSolverTime(m_coarseLevel, 10);
}

} // namespace Nektar::SolverUtils
