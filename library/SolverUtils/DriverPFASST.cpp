///////////////////////////////////////////////////////////////////////////////
//
// File DriverPFASST.cpp
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
// Description: Driver class for the PFASST solver
//
///////////////////////////////////////////////////////////////////////////////

#include <iomanip>

#include <LibUtilities/BasicUtils/Timer.h>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <SolverUtils/DriverPFASST.h>
#include <boost/format.hpp>

namespace Nektar::SolverUtils
{
std::string DriverPFASST::className =
    GetDriverFactory().RegisterCreatorFunction("PFASST", DriverPFASST::create);
std::string DriverPFASST::driverLookupId =
    LibUtilities::SessionReader::RegisterEnumValue("Driver", "PFASST", 0);

/**
 *
 */
DriverPFASST::DriverPFASST(const LibUtilities::SessionReaderSharedPtr pSession,
                           const SpatialDomains::MeshGraphSharedPtr pGraph)
    : DriverParallelInTime(pSession, pGraph)
{
}

/**
 *
 */
void DriverPFASST::v_InitObject(std::ostream &out)
{
    DriverParallelInTime::v_InitObject(out);
    DriverParallelInTime::GetParametersFromSession();
    DriverParallelInTime::PrintSolverInfo(out);
    DriverParallelInTime::InitialiseEqSystem(true);
    InitialiseSDCScheme();
    SetTimeInterpolator();
}

/**
 *
 */
void DriverPFASST::v_Execute([[maybe_unused]] std::ostream &out)
{
    size_t step     = m_chunkRank;
    m_totalTime     = m_timestep[0] * m_nsteps[0];
    m_numWindowsPIT = m_nsteps[0] / m_numChunks;
    m_chunkTime     = m_timestep[0];

    // Start iteration windows.
    m_comm->GetTimeComm()->Block();
    m_CPUtime = 0.0;
    for (size_t w = 0; w < m_numWindowsPIT; w++)
    {
        m_timer.Start();
        PrintHeader((boost::format("WINDOWS #%1%") % (w + 1)).str(), '*');

        // Compute initial guess for coarse solver.
        m_time = (w * m_numChunks) * m_chunkTime;
        ResidualEval(m_time, m_nTimeLevel - 1, 0);
        PropagateQuadratureSolutionAndResidual(m_nTimeLevel - 1, 0);
        for (size_t k = 0; k < m_chunkRank; k++)
        {
            RunSweep(m_time, m_nTimeLevel - 1);
            UpdateFirstQuadrature(m_nTimeLevel - 1);
            PropagateQuadratureSolutionAndResidual(m_nTimeLevel - 1, 0);
            m_time += m_chunkTime;
        }
        RunSweep(m_time, m_nTimeLevel - 1);

        // Interpolate coarse solution and residual to fine.
        for (size_t timeLevel = m_nTimeLevel - 1; timeLevel > 0; timeLevel--)
        {
            InterpolateSolution(timeLevel);
            InterpolateResidual(timeLevel);
        }

        // Start CorrectInitialPFASST iteration.
        size_t k            = 0;
        int convergenceCurr = 0;
        std::vector<int> convergencePrev(m_nTimeLevel, m_chunkRank == 0);
        while (k < m_iterMaxPIT && !convergenceCurr)
        {
            // The PFASST implementation follow "Bolten, M., Moser, D., & Speck,
            // R. (2017). A multigrid perspective on the parallel full
            // approximation scheme in space and time. Numerical Linear Algebra
            // with Applications, 24(6)".

            if (m_chunkRank == m_numChunks - 1 &&
                m_comm->GetSpaceComm()->GetRank() == 0)
            {
                std::cout << "Iteration " << k + 1 << std::endl << std::flush;
            }

            // Fine-to-coarse
            for (size_t timeLevel = 0; timeLevel < m_nTimeLevel - 1;
                 timeLevel++)
            {
                // Performe sweep (parallel-in-time).
                RunSweep(m_time, timeLevel, true);

                // Compute FAS correction (parallel-in-time).
                RestrictSolution(timeLevel);
                RestrictResidual(timeLevel);
                IntegratedResidualEval(timeLevel);
                IntegratedResidualEval(timeLevel + 1);
                ComputeFASCorrection(timeLevel + 1);

                // Check convergence.
                if (timeLevel == 0)
                {
                    // Evaluate SDC residual norm.
                    EvaluateSDCResidualNorm(timeLevel);

                    // Display L2norm.
                    PrintErrorNorm(timeLevel, true);
                }
            }

            // Perform coarse sweep (serial-in-time).
            RecvFromPreviousProc(m_SDCSolver[m_nTimeLevel - 1]
                                     ->UpdateFirstQuadratureSolutionVector(),
                                 convergencePrev[m_nTimeLevel - 1]);
            RunSweep(m_time, m_nTimeLevel - 1, true);
            convergenceCurr = (vL2ErrorMax() < m_tolerPIT &&
                               convergencePrev[m_nTimeLevel - 1]);
            SendToNextProc(m_SDCSolver[m_nTimeLevel - 1]
                               ->UpdateLastQuadratureSolutionVector(),
                           convergenceCurr);

            // Coarse-to-fine.
            for (size_t timeLevel = m_nTimeLevel - 1; timeLevel > 0;
                 timeLevel--)
            {
                // Correct solution and residual.
                /*SendToNextProc(m_SDCSolver[timeLevel - 1]
                                   ->UpdateLastQuadratureSolutionVector(),
                               convergenceCurr);*/
                CorrectSolution(timeLevel - 1);
                CorrectResidual(timeLevel - 1);
                /*RecvFromPreviousProc(
                    m_SDCSolver[timeLevel - 1]
                        ->UpdateFirstQuadratureSolutionVector(),
                    convergencePrev[timeLevel - 1]);
                CorrectInitialSolution(timeLevel - 1);*/
                if (timeLevel - 1 != 0)
                {
                    RunSweep(m_time, timeLevel - 1, true);
                }
            }
            k++;
        }

        // Apply windowing.
        if (w < m_numWindowsPIT - 1)
        {
            ApplyWindowing();
        }
        m_timer.Stop();

        // Update field and write output.
        WriteOutput(step, m_time);
        step += m_numChunks;
        m_CPUtime += m_timer.Elapsed().count();
    }

    m_comm->GetTimeComm()->Block();
    PrintHeader("SUMMARY", '*');
    EvaluateExactSolution(0, m_time + m_chunkTime);
    SolutionConvergenceSummary(0);
    SpeedUpAnalysis();
}

/**
 *
 */
void DriverPFASST::AssertParameters(void)
{
    // Assert time-stepping parameters.
    for (size_t timeLevel = 0; timeLevel < m_nTimeLevel; timeLevel++)
    {
        ASSERTL0(
            m_nsteps[timeLevel] % m_numChunks == 0,
            "Total number of steps should be divisible by number of chunks.");

        ASSERTL0(m_timestep[0] == m_timestep[timeLevel],
                 "All SDC levels should have the same timestep");

        ASSERTL0(m_nsteps[0] == m_nsteps[timeLevel],
                 "All SDC levels should have the same timestep");
    }

    // Assert I/O parameters.
    if (m_EqSys[0]->GetCheckpointSteps())
    {
        ASSERTL0(m_nsteps[0] % m_EqSys[0]->GetCheckpointSteps() == 0,
                 "number of IO_CheckSteps should divide number of steps "
                 "per time chunk");
    }

    if (m_EqSys[0]->GetInfoSteps())
    {
        ASSERTL0(m_nsteps[0] % m_EqSys[0]->GetInfoSteps() == 0,
                 "number of IO_InfoSteps should divide number of steps "
                 "per time chunk");
    }
}

/**
 *
 */
void DriverPFASST::InitialiseSDCScheme(void)
{
    m_SDCSolver =
        Array<OneD, std::shared_ptr<LibUtilities::TimeIntegrationSchemeSDC>>(
            m_nTimeLevel);
    for (size_t timeLevel = 0; timeLevel < m_nTimeLevel; timeLevel++)
    {
        // Cast pointer for TimeIntegrationSchemeSDC.
        m_SDCSolver[timeLevel] =
            std::dynamic_pointer_cast<LibUtilities::TimeIntegrationSchemeSDC>(
                m_EqSys[timeLevel]->GetTimeIntegrationScheme());

        // Assert if a SDC time-integration is used.
        ASSERTL0(m_SDCSolver[timeLevel] != nullptr,
                 "Should only be run with a SDC method");

        // Order storage to list time-integrated fields first.
        Array<OneD, Array<OneD, NekDouble>> fields(m_nVar);
        for (size_t i = 0; i < m_nVar; ++i)
        {
            fields[i] = m_EqSys[timeLevel]->UpdatePhysField(i);
        }

        // Initialize SDC scheme.
        m_SDCSolver[timeLevel]->SetPFASST(timeLevel != 0);
        m_SDCSolver[timeLevel]->InitializeScheme(
            m_timestep[timeLevel], fields, 0.0,
            m_EqSys[timeLevel]->GetTimeIntegrationSchemeOperators());
    }

    // Alocate memory.
    m_QuadPts = Array<OneD, size_t>(m_nTimeLevel);
    for (size_t timeLevel = 0; timeLevel < m_nTimeLevel; timeLevel++)
    {
        m_QuadPts[timeLevel] = m_SDCSolver[timeLevel]->GetQuadPtsNumber();
    }

    m_solutionRest = Array<OneD, SDCarray>(m_nTimeLevel - 1);
    m_residualRest = Array<OneD, SDCarray>(m_nTimeLevel - 1);
    m_integralRest = Array<OneD, SDCarray>(m_nTimeLevel - 1);
    m_correction   = Array<OneD, SDCarray>(m_nTimeLevel - 1);
    m_storage      = Array<OneD, SDCarray>(m_nTimeLevel - 1);
    for (size_t timeLevel = 0; timeLevel < m_nTimeLevel - 1; timeLevel++)
    {
        m_solutionRest[timeLevel] = SDCarray(m_QuadPts[timeLevel + 1]);
        m_residualRest[timeLevel] = SDCarray(m_QuadPts[timeLevel + 1]);
        m_integralRest[timeLevel] = SDCarray(m_QuadPts[timeLevel + 1]);
        m_correction[timeLevel]   = SDCarray(m_QuadPts[timeLevel + 1]);
        m_storage[timeLevel]      = SDCarray(m_QuadPts[timeLevel + 1]);
        for (size_t n = 0; n < m_QuadPts[timeLevel + 1]; ++n)
        {
            m_solutionRest[timeLevel][n] =
                Array<OneD, Array<OneD, NekDouble>>(m_nVar);
            m_residualRest[timeLevel][n] =
                Array<OneD, Array<OneD, NekDouble>>(m_nVar);
            m_integralRest[timeLevel][n] =
                Array<OneD, Array<OneD, NekDouble>>(m_nVar);
            m_correction[timeLevel][n] =
                Array<OneD, Array<OneD, NekDouble>>(m_nVar);
            m_storage[timeLevel][n] =
                Array<OneD, Array<OneD, NekDouble>>(m_nVar);
            for (size_t i = 0; i < m_nVar; ++i)
            {
                m_solutionRest[timeLevel][n][i] =
                    Array<OneD, NekDouble>(m_npts[timeLevel + 1], 0.0);
                m_residualRest[timeLevel][n][i] =
                    Array<OneD, NekDouble>(m_npts[timeLevel + 1], 0.0);
                m_integralRest[timeLevel][n][i] =
                    Array<OneD, NekDouble>(m_npts[timeLevel + 1], 0.0);
                m_correction[timeLevel][n][i] =
                    Array<OneD, NekDouble>(m_npts[timeLevel + 1], 0.0);
                m_storage[timeLevel][n][i] =
                    Array<OneD, NekDouble>(m_npts[timeLevel], 0.0);
            }
        }
    }
}

/**
 *
 */
void DriverPFASST::SetTimeInterpolator(void)
{
    // Initialize time interpolator.
    m_ImatFtoC = Array<OneD, Array<OneD, NekDouble>>(m_nTimeLevel - 1);
    m_ImatCtoF = Array<OneD, Array<OneD, NekDouble>>(m_nTimeLevel - 1);
    for (size_t timeLevel = 0; timeLevel < m_nTimeLevel - 1; timeLevel++)
    {
        LibUtilities::PointsKey fpoints =
            m_SDCSolver[timeLevel]->GetPointsKey();
        LibUtilities::PointsKey cpoints =
            m_SDCSolver[timeLevel + 1]->GetPointsKey();
        DNekMatSharedPtr ImatFtoC =
            LibUtilities::PointsManager()[fpoints]->GetI(cpoints);
        DNekMatSharedPtr ImatCtoF =
            LibUtilities::PointsManager()[cpoints]->GetI(fpoints);
        m_ImatFtoC[timeLevel] = Array<OneD, NekDouble>(
            m_QuadPts[timeLevel] * m_QuadPts[timeLevel + 1], 0.0);
        m_ImatCtoF[timeLevel] = Array<OneD, NekDouble>(
            m_QuadPts[timeLevel] * m_QuadPts[timeLevel + 1], 0.0);

        // Determine if Radau quadrature are used.
        size_t i0 = m_SDCSolver[timeLevel]->HasFirstQuadrature() ? 0 : 1;
        size_t j0 = m_SDCSolver[timeLevel + 1]->HasFirstQuadrature() ? 0 : 1;

        // Adapt fine to coarse time interpolator.
        for (size_t i = i0; i < m_QuadPts[timeLevel]; ++i)
        {
            for (size_t j = j0; j < m_QuadPts[timeLevel + 1]; ++j)
            {
                m_ImatFtoC[timeLevel][i * m_QuadPts[timeLevel + 1] + j] =
                    (ImatFtoC->GetPtr())[(i - i0) *
                                             (m_QuadPts[timeLevel + 1] - j0) +
                                         (j - j0)];
            }
        }
        if (j0 == 1)
        {
            m_ImatFtoC[timeLevel][0] = 1.0;
        }

        // Adapt coarse to fine time interpolator.
        for (size_t j = j0; j < m_QuadPts[timeLevel + 1]; ++j)
        {
            for (size_t i = i0; i < m_QuadPts[timeLevel]; ++i)
            {
                m_ImatCtoF[timeLevel][j * m_QuadPts[timeLevel] + i] =
                    (ImatCtoF
                         ->GetPtr())[(j - j0) * (m_QuadPts[timeLevel] - i0) +
                                     (i - i0)];
            }
        }
        if (i0 == 1)
        {
            m_ImatCtoF[timeLevel][0] = 1.0;
        }
    }
}

/**
 *
 */
bool DriverPFASST::IsNotInitialCondition(const size_t n)
{
    return !(n == 0 && m_chunkRank == 0);
}

/**
 *
 */
void DriverPFASST::PropagateQuadratureSolutionAndResidual(
    const size_t timeLevel, const size_t index)
{
    for (size_t n = 0; n < m_QuadPts[timeLevel]; ++n)
    {
        if (n != index)
        {
            for (size_t i = 0; i < m_nVar; ++i)
            {
                Vmath::Vcopy(
                    m_npts[timeLevel],
                    m_SDCSolver[timeLevel]->GetSolutionVector()[index][i], 1,
                    m_SDCSolver[timeLevel]->UpdateSolutionVector()[n][i], 1);
                Vmath::Vcopy(
                    m_npts[timeLevel],
                    m_SDCSolver[timeLevel]->GetResidualVector()[index][i], 1,
                    m_SDCSolver[timeLevel]->UpdateResidualVector()[n][i], 1);
            }
        }
    }
}

/**
 *
 */
void DriverPFASST::UpdateFirstQuadrature(const size_t timeLevel)
{
    m_SDCSolver[timeLevel]->UpdateFirstQuadrature();
}

/**
 *
 */
void DriverPFASST::RunSweep(const NekDouble time, const size_t timeLevel,
                            const bool update)
{
    size_t niter = m_SDCSolver[timeLevel]->GetOrder();

    if (update == true)
    {
        ResidualEval(m_time, timeLevel, 0);
    }

    // Start SDC iteration loop.
    m_SDCSolver[timeLevel]->SetTime(time);
    for (size_t k = 0; k < niter; k++)
    {
        m_SDCSolver[timeLevel]->SDCIterationLoop(m_chunkTime);
    }

    // Update last quadrature point.
    m_SDCSolver[timeLevel]->UpdateLastQuadrature();
}

/**
 *
 */
void DriverPFASST::ResidualEval(const NekDouble time, const size_t timeLevel,
                                const size_t n)
{
    m_SDCSolver[timeLevel]->SetTime(time);
    m_SDCSolver[timeLevel]->ResidualEval(m_chunkTime, n);
}

/**
 *
 */
void DriverPFASST::ResidualEval(const NekDouble time, const size_t timeLevel)
{
    m_SDCSolver[timeLevel]->SetTime(time);
    m_SDCSolver[timeLevel]->ResidualEval(m_chunkTime);
}

/**
 *
 */
void DriverPFASST::IntegratedResidualEval(const size_t timeLevel)
{
    m_SDCSolver[timeLevel]->UpdateIntegratedResidualQFint(m_chunkTime);
}

/**
 *
 */
void DriverPFASST::Interpolate(const size_t coarseLevel, const SDCarray &in,
                               const size_t fineLevel, SDCarray &out,
                               bool forced)
{
    // Interpolate solution in space.
    for (size_t n = 0; n < m_QuadPts[coarseLevel]; ++n)
    {
        Interpolate(m_EqSys[coarseLevel]->UpdateFields(),
                    m_EqSys[fineLevel]->UpdateFields(), in[n],
                    m_storage[fineLevel][n]);
    }

    // Interpolate solution in time.
    for (size_t n = 0; n < m_QuadPts[fineLevel]; ++n)
    {
        if (forced || IsNotInitialCondition(n))
        {
            for (size_t i = 0; i < m_nVar; ++i)
            {
                Vmath::Smul(m_npts[fineLevel], m_ImatCtoF[fineLevel][n],
                            m_storage[fineLevel][0][i], 1, out[n][i], 1);
                for (size_t k = 1; k < m_QuadPts[coarseLevel]; ++k)
                {
                    size_t index = k * m_QuadPts[fineLevel] + n;
                    Vmath::Svtvp(m_npts[fineLevel],
                                 m_ImatCtoF[fineLevel][index],
                                 m_storage[fineLevel][k][i], 1, out[n][i], 1,
                                 out[n][i], 1);
                }
            }
        }
    }
}

/**
 *
 */
void DriverPFASST::InterpolateSolution(const size_t timeLevel)
{
    size_t coarseLevel = timeLevel;
    size_t fineLevel   = timeLevel - 1;

    Interpolate(coarseLevel, m_SDCSolver[coarseLevel]->GetSolutionVector(),
                fineLevel, m_SDCSolver[fineLevel]->UpdateSolutionVector(),
                false);
}

/**
 *
 */
void DriverPFASST::InterpolateResidual(const size_t timeLevel)
{
    size_t coarseLevel = timeLevel;
    size_t fineLevel   = timeLevel - 1;

    Interpolate(coarseLevel, m_SDCSolver[coarseLevel]->GetResidualVector(),
                fineLevel, m_SDCSolver[fineLevel]->UpdateResidualVector(),
                true);
}

/**
 *
 */
void DriverPFASST::Restrict(const size_t fineLevel, const SDCarray &in,
                            const size_t coarseLevel, SDCarray &out)
{
    // Restrict fine solution in time.
    for (size_t n = 0; n < m_QuadPts[coarseLevel]; ++n)
    {
        for (size_t i = 0; i < m_nVar; ++i)
        {
            Vmath::Smul(m_npts[fineLevel], m_ImatFtoC[fineLevel][n], in[0][i],
                        1, m_storage[fineLevel][n][i], 1);
            for (size_t k = 1; k < m_QuadPts[fineLevel]; ++k)
            {
                size_t index = k * m_QuadPts[coarseLevel] + n;
                Vmath::Svtvp(m_npts[fineLevel], m_ImatFtoC[fineLevel][index],
                             in[k][i], 1, m_storage[fineLevel][n][i], 1,
                             m_storage[fineLevel][n][i], 1);
            }
        }
    }

    // Restrict fine solution in space.
    for (size_t n = 0; n < m_QuadPts[coarseLevel]; ++n)
    {
        Interpolate(m_EqSys[fineLevel]->UpdateFields(),
                    m_EqSys[coarseLevel]->UpdateFields(),
                    m_storage[fineLevel][n], out[n]);
    }
}

/**
 *
 */
void DriverPFASST::RestrictSolution(const size_t timeLevel)
{
    size_t fineLevel   = timeLevel;
    size_t coarseLevel = timeLevel + 1;

    Restrict(fineLevel, m_SDCSolver[fineLevel]->GetSolutionVector(),
             coarseLevel, m_SDCSolver[coarseLevel]->UpdateSolutionVector());

    for (size_t n = 0; n < m_QuadPts[coarseLevel]; ++n)
    {
        CopySolutionVector(m_SDCSolver[coarseLevel]->GetSolutionVector()[n],
                           m_solutionRest[fineLevel][n]);
    }
}

/**
 *
 */
void DriverPFASST::RestrictResidual(const size_t timeLevel)
{
    size_t fineLevel   = timeLevel;
    size_t coarseLevel = timeLevel + 1;

    ResidualEval(m_time, coarseLevel);

    if (!m_updateResidual)
    {
        for (size_t n = 0; n < m_QuadPts[coarseLevel]; ++n)
        {
            CopySolutionVector(m_SDCSolver[coarseLevel]->GetResidualVector()[n],
                               m_residualRest[fineLevel][n]);
        }
    }
}

/**
 *
 */
void DriverPFASST::ComputeFASCorrection(const size_t timeLevel)
{
    size_t fineLevel   = timeLevel - 1;
    size_t coarseLevel = timeLevel;

    if (fineLevel != 0)
    {
        // Restrict fine FAS correction term
        Restrict(fineLevel, m_SDCSolver[fineLevel]->UpdateFAScorrectionVector(),
                 coarseLevel,
                 m_SDCSolver[coarseLevel]->UpdateFAScorrectionVector());

        // Restrict fine integrated residual.
        Restrict(fineLevel,
                 m_SDCSolver[fineLevel]->GetIntegratedResidualVector(),
                 coarseLevel, m_integralRest[fineLevel]);
    }
    else
    {
        // Restrict fine integrated residual.
        Restrict(
            fineLevel, m_SDCSolver[fineLevel]->GetIntegratedResidualVector(),
            coarseLevel, m_SDCSolver[coarseLevel]->UpdateFAScorrectionVector());
    }

    // Compute coarse FAS correction terms.
    for (size_t n = 0; n < m_QuadPts[coarseLevel]; ++n)
    {
        for (size_t i = 0; i < m_nVar; ++i)
        {
            if (fineLevel != 0)
            {
                Vmath::Vadd(
                    m_npts[coarseLevel],
                    m_SDCSolver[coarseLevel]->GetFAScorrectionVector()[n][i], 1,
                    m_integralRest[fineLevel][n][i], 1,
                    m_SDCSolver[coarseLevel]->UpdateFAScorrectionVector()[n][i],
                    1);
            }

            Vmath::Vsub(
                m_npts[coarseLevel],
                m_SDCSolver[coarseLevel]->GetFAScorrectionVector()[n][i], 1,
                m_SDCSolver[coarseLevel]->GetIntegratedResidualVector()[n][i],
                1, m_SDCSolver[coarseLevel]->UpdateFAScorrectionVector()[n][i],
                1);
        }
    }
}

/**
 *
 */
void DriverPFASST::Correct(const size_t coarseLevel,
                           const Array<OneD, Array<OneD, NekDouble>> &in,
                           const size_t fineLevel,
                           Array<OneD, Array<OneD, NekDouble>> &out,
                           bool forced)
{
    if (forced || IsNotInitialCondition(0))
    {
        // Compute difference between coarse solution and restricted
        // solution.
        Interpolate(m_EqSys[fineLevel]->UpdateFields(),
                    m_EqSys[coarseLevel]->UpdateFields(), out,
                    m_correction[fineLevel][0]);
        for (size_t i = 0; i < m_nVar; ++i)
        {
            Vmath::Vsub(m_npts[coarseLevel], in[i], 1,
                        m_correction[fineLevel][0][i], 1,
                        m_correction[fineLevel][0][i], 1);
        }

        // Add correction to fine solution.
        Interpolate(m_EqSys[coarseLevel]->UpdateFields(),
                    m_EqSys[fineLevel]->UpdateFields(),
                    m_correction[fineLevel][0], m_storage[fineLevel][0]);
        for (size_t i = 0; i < m_nVar; ++i)
        {
            Vmath::Vadd(m_npts[fineLevel], m_storage[fineLevel][0][i], 1,
                        out[i], 1, out[i], 1);
        }
    }
}

/**
 *
 */
void DriverPFASST::CorrectInitialSolution(const size_t timeLevel)
{
    size_t fineLevel   = timeLevel;
    size_t coarseLevel = timeLevel + 1;

    Correct(coarseLevel, m_SDCSolver[coarseLevel]->GetSolutionVector()[0],
            fineLevel, m_SDCSolver[fineLevel]->UpdateSolutionVector()[0],
            false);
}

/**
 *
 */
void DriverPFASST::CorrectInitialResidual(const size_t timeLevel)
{
    size_t fineLevel   = timeLevel;
    size_t coarseLevel = timeLevel + 1;

    Correct(coarseLevel, m_SDCSolver[coarseLevel]->GetResidualVector()[0],
            fineLevel, m_SDCSolver[fineLevel]->UpdateResidualVector()[0],
            false);
}

/**
 *
 */
void DriverPFASST::Correct(const size_t coarseLevel, const SDCarray &rest,
                           const SDCarray &in, const size_t fineLevel,
                           SDCarray &out, bool forced)
{
    // Compute difference between coarse solution and restricted
    // solution.
    for (size_t n = 0; n < m_QuadPts[coarseLevel]; ++n)
    {
        if (forced || IsNotInitialCondition(n))
        {
            for (size_t i = 0; i < m_nVar; ++i)
            {
                Vmath::Vsub(m_npts[coarseLevel], in[n][i], 1, rest[n][i], 1,
                            m_correction[fineLevel][n][i], 1);
            }
        }
        else
        {
            for (size_t i = 0; i < m_nVar; ++i)
            {
                Vmath::Zero(m_npts[coarseLevel], m_correction[fineLevel][n][i],
                            1);
            }
        }
    }

    // Interpolate coarse solution delta in space.
    for (size_t n = 0; n < m_QuadPts[coarseLevel]; ++n)
    {
        Interpolate(m_EqSys[coarseLevel]->UpdateFields(),
                    m_EqSys[fineLevel]->UpdateFields(),
                    m_correction[fineLevel][n], m_storage[fineLevel][n]);
    }

    // Interpolate coarse solution delta in time and correct fine solution.
    for (size_t n = 0; n < m_QuadPts[fineLevel]; ++n)
    {
        if (forced || IsNotInitialCondition(n))
        {
            for (size_t i = 0; i < m_nVar; ++i)
            {
                for (size_t k = 0; k < m_QuadPts[coarseLevel]; ++k)
                {
                    size_t index = k * m_QuadPts[fineLevel] + n;
                    Vmath::Svtvp(m_npts[fineLevel],
                                 m_ImatCtoF[fineLevel][index],
                                 m_storage[fineLevel][k][i], 1, out[n][i], 1,
                                 out[n][i], 1);
                }
            }
        }
    }
}

/**
 *
 */
void DriverPFASST::CorrectSolution(const size_t timeLevel)
{
    size_t coarseLevel = timeLevel + 1;
    size_t fineLevel   = timeLevel;

    Correct(coarseLevel, m_solutionRest[fineLevel],
            m_SDCSolver[coarseLevel]->GetSolutionVector(), fineLevel,
            m_SDCSolver[fineLevel]->UpdateSolutionVector(), false);
}

/**
 *
 */
void DriverPFASST::CorrectResidual(const size_t timeLevel)
{
    size_t coarseLevel = timeLevel + 1;
    size_t fineLevel   = timeLevel;

    // Evaluate fine residual.
    if (m_updateResidual)
    {
        for (size_t n = 1; n < m_QuadPts[fineLevel]; ++n)
        {
            ResidualEval(fineLevel, n);
        }
    }
    // Correct fine residual.
    else
    {
        Correct(coarseLevel, m_residualRest[fineLevel],
                m_SDCSolver[coarseLevel]->GetResidualVector(), fineLevel,
                m_SDCSolver[fineLevel]->UpdateResidualVector(), false);
    }
}

/**
 *
 */
void DriverPFASST::ApplyWindowing(void)
{
    // Use last chunk solution as initial condition for the next window.
    if (m_chunkRank == m_numChunks - 1)
    {
        UpdateFirstQuadrature(0);
        for (size_t timeLevel = 0; timeLevel < m_nTimeLevel - 1; timeLevel++)
        {
            Interpolate(m_EqSys[timeLevel]->UpdateFields(),
                        m_EqSys[timeLevel + 1]->UpdateFields(),
                        m_SDCSolver[timeLevel]->GetSolutionVector()[0],
                        m_SDCSolver[timeLevel + 1]->UpdateSolutionVector()[0]);
        }
    }

    // Broadcast I.C. for windowing.
    for (size_t timeLevel = 0; timeLevel < m_nTimeLevel; timeLevel++)
    {
        for (size_t i = 0; i < m_nVar; ++i)
        {
            m_comm->GetTimeComm()->Bcast(
                m_SDCSolver[timeLevel]->UpdateSolutionVector()[0][i],
                m_numChunks - 1);
        }
    }
}

/**
 *
 */
void DriverPFASST::EvaluateSDCResidualNorm(const size_t timeLevel)
{
    // Compute SDC residual norm
    for (size_t i = 0; i < m_nVar; ++i)
    {
        Vmath::Vadd(
            m_npts[timeLevel],
            m_SDCSolver[timeLevel]->GetSolutionVector()[0][i], 1,
            m_SDCSolver[timeLevel]
                ->GetIntegratedResidualVector()[m_QuadPts[timeLevel] - 1][i],
            1, m_exactsoln[i], 1);
        m_EqSys[timeLevel]->CopyToPhysField(
            i, m_SDCSolver[timeLevel]
                   ->GetSolutionVector()[m_QuadPts[timeLevel] - 1][i]);
        m_vL2Errors[i]   = m_EqSys[timeLevel]->L2Error(i, m_exactsoln[i], 1);
        m_vLinfErrors[i] = m_EqSys[timeLevel]->LinfError(i, m_exactsoln[i]);
    }
}

/**
 *
 */
void DriverPFASST::WriteOutput(const size_t step, const NekDouble time)
{
    size_t timeLevel           = 0;
    static size_t IOChkStep    = m_EqSys[0]->GetCheckpointSteps()
                                     ? m_EqSys[0]->GetCheckpointSteps()
                                     : m_nsteps[timeLevel];
    static std::string dirname = m_session->GetSessionName() + ".pit";

    if ((step + 1) % IOChkStep == 0)
    {
        // Compute checkpoint index.
        size_t nchk = (step + 1) / IOChkStep;

        // Create directory if does not exist.
        if (!fs::is_directory(dirname))
        {
            fs::create_directory(dirname);
        }

        // Update solution field.
        UpdateFieldCoeffs(
            timeLevel,
            m_SDCSolver[timeLevel]->UpdateLastQuadratureSolutionVector());

        // Set filename.
        std::string filename = dirname + "/" + m_session->GetSessionName() +
                               "_" + std::to_string(nchk) + ".fld";

        // Set time metadata.
        m_EqSys[timeLevel]->SetTime(time);

        // Write checkpoint.
        m_EqSys[timeLevel]->WriteFld(filename);
    }
}

/**
 *
 */
void DriverPFASST::SpeedUpAnalysis(void)
{
    // Print header.
    PrintHeader("PFASST SPEED-UP ANALYSIS", '*');

    // Mean communication time.
    NekDouble commTime = EstimateCommunicationTime();
    PrintHeader("Mean Communication Time = " +
                    (boost::format("%1$.6e") % commTime).str() + "s",
                '-');

    // Mean FAS correction time.
    NekDouble fasTime = EstimateFASCorrectionTime();
    PrintHeader("Mean FAS Correction Time = " +
                    (boost::format("%1$.6e") % fasTime).str() + "s",
                '-');

    // Mean coarse solver time.
    NekDouble coarseSolveTime = EstimateSolverTime(m_nTimeLevel - 1);
    PrintHeader("Mean Coarse Solve Time = " +
                    (boost::format("%1$.6e") % coarseSolveTime).str() + "s",
                '-');

    // Mean fine solver time.
    NekDouble fineSolveTime = EstimateSolverTime(0);
    PrintHeader("Mean Fine Solve Time = " +
                    (boost::format("%1$.6e") % fineSolveTime).str() + "s",
                '-');

    // Mean predictor time.
    NekDouble predictorTime = EstimatePredictorTime();
    PrintHeader("Mean Predictor Time = " +
                    (boost::format("%1$.6e") % predictorTime).str() + "s",
                '-');

    // Mean overhead time.
    NekDouble overheadTime = EstimateOverheadTime();
    PrintHeader("Mean Overhead Time = " +
                    (boost::format("%1$.6e") % overheadTime).str() + "s",
                '-');

    // Print speedup time.
    PrintSpeedUp(fineSolveTime, coarseSolveTime, fasTime, commTime,
                 predictorTime, overheadTime);
}

/**
 *
 */
void DriverPFASST::PrintSpeedUp(NekDouble fineSolveTime,
                                NekDouble coarseSolveTime, NekDouble fasTime,
                                NekDouble commTime, NekDouble predictTime,
                                NekDouble overheadTime)
{
    if (m_chunkRank == m_numChunks - 1 &&
        m_comm->GetSpaceComm()->GetRank() == 0)
    {
        // Print maximum theoretical speed-up.
        PrintHeader("Maximum Speed-up", '-');
        for (size_t k = 1; k <= m_numChunks; k++)
        {
            NekDouble speedup =
                ComputeSpeedUp(k, fineSolveTime, coarseSolveTime, 0.0, 0.0,
                               predictTime, overheadTime);
            std::cout << "Speed-up (" << k << ") = " << speedup << std::endl
                      << std::flush;
        }

        // Print speed-up with fas correction.
        PrintHeader("Speed-up with fas", '-');
        for (size_t k = 1; k <= m_numChunks; k++)
        {
            NekDouble speedup =
                ComputeSpeedUp(k, fineSolveTime, coarseSolveTime, fasTime, 0.0,
                               predictTime, overheadTime);
            std::cout << "Speed-up (" << k << ") = " << speedup << std::endl
                      << std::flush;
        }

        // Print speed-up with fas correction and communication.
        PrintHeader("Speed-up with comm. and fas", '-');
        for (size_t k = 1; k <= m_numChunks; k++)
        {
            NekDouble speedup =
                ComputeSpeedUp(k, fineSolveTime, coarseSolveTime, fasTime,
                               commTime, predictTime, overheadTime);
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
NekDouble DriverPFASST::ComputeSpeedUp(const size_t iter,
                                       NekDouble fineSolveTime,
                                       NekDouble coarseSolveTime,
                                       NekDouble fasTime, NekDouble commTime,
                                       NekDouble predictorTime,
                                       NekDouble overheadTime)
{
    // The speed-up estimate based on "Emmett, M., & Minion, M. (2012). Toward
    // an efficient parallel in time method for partial differential equations.
    // Communications in Applied Mathematics and Computational Science, 7(1),
    // 105-132" and on "Lunet, T., Bodart, J., Gratton, S., &
    // Vasseur, X. (2018). Time-parallel simulation of the decay of homogeneous
    // turbulence using parareal with spatial coarsening. Computing and
    // Visualization in Science, 19, 31-44".

    size_t Kiter             = m_SDCSolver[0]->GetMaxOrder();
    size_t nComm             = (iter * (2 * m_numChunks - iter - 1)) / 2;
    NekDouble ratio          = double(iter) / m_numChunks;
    NekDouble ratioPredictor = predictorTime / fineSolveTime;
    NekDouble ratioSolve     = coarseSolveTime / fineSolveTime;
    NekDouble ratioFAS       = fasTime / fineSolveTime;
    NekDouble ratioComm      = commTime / fineSolveTime;
    NekDouble ratioOverhead  = overheadTime / fineSolveTime;

    // Speed-up relative to SDC.
    return Kiter / (ratioPredictor + ratio * (1.0 + ratioSolve + ratioFAS) +
                    (ratioComm * nComm + ratioOverhead) / m_numChunks);
    // Speed-up relative to MLSDC.
    /*return Kiter * (1.0 + ratioSolve + ratioFAS + ratioOverhead / m_numChunks)
       / (ratioPredictor + ratio * (1.0 + ratioSolve + ratioFAS) +
            (ratioComm * nComm + ratioOverhead) / m_numChunks);*/
}

/**
 *
 */
NekDouble DriverPFASST::EstimateCommunicationTime(void)
{
    // Allocate memory.
    NekDouble commTime = 0.0;
    Array<OneD, Array<OneD, NekDouble>> buffer1(m_nVar);
    Array<OneD, Array<OneD, NekDouble>> buffer2(m_nVar);
    for (size_t timeLevel = 0; timeLevel < m_nTimeLevel; timeLevel++)
    {
        for (size_t i = 0; i < m_nVar; ++i)
        {
            buffer1[i] = Array<OneD, NekDouble>(m_npts[timeLevel], 0.0);
            buffer2[i] = Array<OneD, NekDouble>(m_npts[timeLevel], 0.0);
        }
        commTime +=
            DriverParallelInTime ::EstimateCommunicationTime(buffer1, buffer2);
    }

    return commTime;
}

/**
 *
 */
NekDouble DriverPFASST::EstimatePredictorTime(void)
{
    // Estimate coarse overhead time.
    size_t niter = 20;
    Nektar::LibUtilities::Timer timer;
    timer.Start();
    for (size_t i = 0; i < niter; i++)
    {
        RunSweep(m_time, m_nTimeLevel - 1);
        UpdateFirstQuadrature(m_nTimeLevel - 1);
        PropagateQuadratureSolutionAndResidual(m_nTimeLevel - 1, 0);
    }
    timer.Stop();
    return timer.Elapsed().count() / niter;
}

/**
 *
 */
NekDouble DriverPFASST::EstimateFASCorrectionTime(void)
{
    // Average restriction time over niter iteration.
    size_t niter = 20;
    Nektar::LibUtilities::Timer timer;
    timer.Start();
    for (size_t n = 0; n < niter; n++)
    {
        for (size_t timeLevel = 0; timeLevel < m_nTimeLevel - 1; timeLevel++)
        {
            if (timeLevel != 0)
            {
                RunSweep(m_time, timeLevel, true);
            }
            RestrictSolution(timeLevel);
            RestrictResidual(timeLevel);
            IntegratedResidualEval(timeLevel);
            IntegratedResidualEval(timeLevel + 1);
            ComputeFASCorrection(timeLevel + 1);
        }

        for (size_t timeLevel = m_nTimeLevel - 1; timeLevel > 0; timeLevel--)
        {
            CorrectSolution(timeLevel - 1);
            CorrectResidual(timeLevel - 1);
            CorrectInitialSolution(timeLevel - 1);
            if (timeLevel - 1 != 0)
            {
                RunSweep(m_time, timeLevel - 1, true);
            }
        }
    }
    timer.Stop();
    return timer.Elapsed().count() / niter;
}

/**
 *
 */
NekDouble DriverPFASST::EstimateSolverTime(size_t timeLevel)
{
    // Estimate solver time.
    size_t niter = 20;
    Nektar::LibUtilities::Timer timer;
    timer.Start();
    for (size_t i = 0; i < niter; i++)
    {
        RunSweep(m_time, timeLevel, true);
    }
    timer.Stop();
    return timer.Elapsed().count() / niter;
}

/**
 *
 */
NekDouble DriverPFASST::EstimateOverheadTime(void)
{
    // Estimate overhead time.
    size_t niter = 20;
    Nektar::LibUtilities::Timer timer;
    timer.Start();
    for (size_t i = 0; i < niter; i++)
    {
        ResidualEval(m_time, m_nTimeLevel - 1, 0);
        PropagateQuadratureSolutionAndResidual(m_nTimeLevel - 1, 0);
        RunSweep(m_time, m_nTimeLevel - 1);
        for (size_t timeLevel = m_nTimeLevel - 1; timeLevel > 0; timeLevel--)
        {
            InterpolateSolution(timeLevel);
            InterpolateResidual(timeLevel);
        }
        ApplyWindowing();
    }
    timer.Stop();
    return timer.Elapsed().count() / niter;
}

} // namespace Nektar::SolverUtils
