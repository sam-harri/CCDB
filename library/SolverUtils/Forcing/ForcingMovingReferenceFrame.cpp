///////////////////////////////////////////////////////////////////////////////
//
// File: ForcingMovingReferenceFrame.cpp
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
// Description: Solving the absolute flow in a moving body frame,
// by adding (U0 + Omega X (x - x0)) . grad u - Omega X u
// as the body force.
// U0 is the translational velocity of the body frame.
// Omega is the angular velocity.
// x0 is the rotation pivot in the body frame.
// All vectors use the basis of the body frame.
// Translational motion is allowed for all dimensions.
// Rotation is not allowed for 1D, 2DH1D, 3DH2D.
// Rotation in z direction is allowed for 2D and 3DH1D.
// Rotation in 3 directions are allowed for 3D.
// TODO: add suport for 3D rotation using Quaternion
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <LibUtilities/BasicUtils/Vmath.hpp>
#include <MultiRegions/ExpList.h>
#include <SolverUtils/Filters/FilterInterfaces.hpp>
#include <SolverUtils/Forcing/ForcingMovingReferenceFrame.h>
#include <boost/format.hpp>
using namespace std;

namespace Nektar::SolverUtils
{

std::string ForcingMovingReferenceFrame::classNameBody =
    GetForcingFactory().RegisterCreatorFunction(
        "MovingReferenceFrame", ForcingMovingReferenceFrame::create,
        "Moving Frame");

/**
 * @brief
 * @param pSession
 * @param pEquation
 */
ForcingMovingReferenceFrame::ForcingMovingReferenceFrame(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const std::weak_ptr<EquationSystem> &pEquation)
    : Forcing(pSession, pEquation)
{
}

ForcingMovingReferenceFrame::~ForcingMovingReferenceFrame(void)
{
    if (m_rank == 0)
    {
        m_outputStream.close();
    }
}

/**
 * @brief Initialise the forcing module
 * @param pFields
 * @param pNumForcingFields
 * @param pForce
 */
void ForcingMovingReferenceFrame::v_InitObject(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    [[maybe_unused]] const unsigned int &pNumForcingFields,
    const TiXmlElement *pForce)
{
    m_session->MatchSolverInfo("Homogeneous", "1D", m_isH1d, false);
    m_session->MatchSolverInfo("Homogeneous", "2D", m_isH2d, false);
    bool singleMode, halfMode;
    m_session->MatchSolverInfo("ModeType", "SingleMode", singleMode, false);
    m_session->MatchSolverInfo("ModeType", "HalfMode", halfMode, false);
    if (singleMode || halfMode)
    {
        m_isH1d = false;
    }

    int npoints   = pFields[0]->GetNpoints();
    m_expdim      = m_isH2d ? 1 : pFields[0]->GetGraph()->GetMeshDimension();
    m_spacedim    = m_expdim + (m_isH1d ? 1 : 0) + (m_isH2d ? 2 : 0);
    m_NumVariable = m_spacedim;

    m_hasPlane0 = true;
    if (m_isH1d)
    {
        m_hasPlane0 = pFields[0]->GetZIDs()[0] == 0;
    }

    // initialize variables
    m_velXYZ      = Array<OneD, NekDouble>(3, 0.0);
    m_velxyz      = Array<OneD, NekDouble>(3, 0.0);
    m_omegaXYZ    = Array<OneD, NekDouble>(3, 0.0);
    m_omegaxyz    = Array<OneD, NekDouble>(3, 0.0);
    m_extForceXYZ = Array<OneD, NekDouble>(6, 0.0);
    m_hasVel      = Array<OneD, bool>(3, false);
    m_hasOmega    = Array<OneD, bool>(3, false);
    m_disp        = Array<OneD, NekDouble>(9, 0.0);
    m_pivotPoint  = Array<OneD, NekDouble>(3, 0.0);
    m_travelWave  = Array<OneD, NekDouble>(3, 0.);

    // initialize time
    NekDouble time = 0.;
    CheckForRestartTime(pFields, time);
    m_startTime = time;
    // load moving frame infomation
    LoadParameters(pForce, time);
    if (m_hasRotation)
    {
        m_coords = Array<OneD, Array<OneD, NekDouble>>(3);
        for (int j = 0; j < m_spacedim; ++j)
        {
            m_coords[j] = Array<OneD, NekDouble>(npoints);
        }
        pFields[0]->GetCoords(m_coords[0], m_coords[1], m_coords[2]);
        // move the origin to the pivot point
        for (int i = 0; i < m_spacedim; ++i)
        {
            Vmath::Sadd(npoints, -m_pivotPoint[i], m_coords[i], 1, m_coords[i],
                        1);
        }
    }

    // init body solver
    m_rank = pFields[0]->GetComm()->GetRank();
    InitBodySolver(pForce, m_spacedim, m_rank, time);

    // Initialize theta with the data from NS class
    // This ensure correct moving coordinate orientation with respect to the
    // stationary inertial frame when we restart the simulation
    auto equ = m_equ.lock();
    ASSERTL0(equ, "Weak pointer to the equation system is expired");
    auto FluidEq = std::dynamic_pointer_cast<FluidInterface>(equ);
    FluidEq->SetMovingFrameDisp(m_disp);

    // initiate the rotation matrix to zero,
    // Note that the rotation matrix is assumed for the rotation around z axis
    // TODO: Generalize this for the 3D case with possiblity of rotation around
    // each of axis. Probabley we only can support rotation around one axis. In
    // that case the generalization means the user can provide Omega for one of
    // x, y or z axis. Not sure how complicated to consider a general axis of
    // rotation
    //
    // Note that these rotation matrices should be extrinsic rotations
    m_ProjMatZ = bn::ublas::zero_matrix<NekDouble>(3, 3);
    // populate the rotation matrix R(z)
    {
        NekDouble sn, cs;
        sn = sin(m_disp[5]);
        cs = cos(m_disp[5]);

        m_ProjMatZ(0, 0) = cs;
        m_ProjMatZ(0, 1) = sn;
        m_ProjMatZ(1, 0) = -1. * sn;
        m_ProjMatZ(1, 1) = cs;
        m_ProjMatZ(2, 2) = 1.0;
    }
    // account for the effect of rotation
    // Omega_X results in having v and w even if not defined by user
    // Omega_Y results in having u and w even if not defined by user
    // Omega_Z results in having u and v even if not defined by user
    //
    for (int i = 0; i < 3; ++i)
    {
        int j = (i + 1) % 3;
        int k = (i + 2) % 3;

        if (m_hasOmega[i])
        {
            m_hasVel[j] = true;
            m_hasVel[k] = true;
        }
    }
    // initialize aeroforce filter
    if (m_hasFreeMotion)
    {
        InitialiseFilter(m_session, pFields, pForce);
    }
}

NekDouble ForcingMovingReferenceFrame::EvaluateExpression(
    std::string expression)
{
    NekDouble value = 0.;
    try
    {
        LibUtilities::Equation expession(m_session->GetInterpreter(),
                                         expression);
        value = expession.Evaluate();
    }
    catch (const std::runtime_error &)
    {
        NEKERROR(ErrorUtil::efatal, "Error evaluating expression" + expression);
    }
    return value;
}

void ForcingMovingReferenceFrame::LoadParameters(const TiXmlElement *pForce,
                                                 const NekDouble time)
{
    const TiXmlElement *funcNameElmt;
    // load frame velocity
    funcNameElmt = pForce->FirstChildElement("FRAMEVELOCITY");
    if (funcNameElmt)
    {
        std::string FuncName = funcNameElmt->GetText();
        ASSERTL0(m_session->DefinesFunction(FuncName),
                 "Function '" + FuncName + "' is not defined in the session.");
        // linear velocity
        for (int i = 0; i < m_spacedim; ++i)
        {
            std::string var = m_session->GetVariable(i);
            if (m_session->DefinesFunction(FuncName, var))
            {
                m_frameVelFunction[i] = m_session->GetFunction(FuncName, var);
                m_hasVel[i]           = true;
            }
        }
        // angular velocities
        m_hasRotation                       = false;
        std::vector<std::string> angularVar = {"Omega_x", "Omega_y", "Omega_z"};
        for (int i = 0; i < 3; ++i)
        {
            std::string var = angularVar[i];
            if (m_session->DefinesFunction(FuncName, var))
            {
                m_frameVelFunction[i + 3] =
                    m_session->GetFunction(FuncName, var);
                m_hasOmega[i] = true;
                m_hasRotation = true;
            }
        }
        // TODO: add the support for general rotation
        for (int i = 0; i < 2; ++i)
        {
            ASSERTL0(!m_hasOmega[i], "Currently only Omega_z is supported");
        }
    }
    if (m_expdim == 1)
    {
        m_hasRotation = false;
    }

    // load external force
    funcNameElmt = pForce->FirstChildElement("EXTERNALFORCE");
    if (funcNameElmt)
    {
        std::string FuncName = funcNameElmt->GetText();
        if (m_session->DefinesFunction(FuncName))
        {
            std::vector<std::string> forceVar = {"fx", "fy", "fz"};
            for (int i = 0; i < m_spacedim; ++i)
            {
                std::string var = forceVar[i];
                if (m_session->DefinesFunction(FuncName, var))
                {
                    m_extForceFunction[i] =
                        m_session->GetFunction(FuncName, var);
                }
            }
            std::vector<std::string> momentVar = {"Mx", "My", "Mz"};
            for (int i = 0; i < 3; ++i)
            {
                std::string var = momentVar[i];
                if (m_session->DefinesFunction(FuncName, var))
                {
                    m_extForceFunction[i + 3] =
                        m_session->GetFunction(FuncName, var);
                }
            }
        }
    }

    // load pitching pivot
    const TiXmlElement *pivotElmt = pForce->FirstChildElement("PIVOTPOINT");
    if (pivotElmt)
    {
        std::vector<std::string> values;
        std::string mssgStr = pivotElmt->GetText();
        ParseUtils::GenerateVector(mssgStr, values);
        for (int j = 0; j < m_spacedim; ++j)
        {
            m_pivotPoint[j] = EvaluateExpression(values[j]);
            m_disp[6 + j]   = m_pivotPoint[j];
        }
    }

    // load travelling wave speed
    const TiXmlElement *TWSElmt =
        pForce->FirstChildElement("TRAVELINGWAVESPEED");
    if (TWSElmt)
    {
        std::vector<std::string> values;
        std::string mssgStr = TWSElmt->GetText();
        ParseUtils::GenerateVector(mssgStr, values);
        for (int j = 0; j < m_spacedim; ++j)
        {
            NekDouble tmp = EvaluateExpression(values[j]);
            if (fabs(tmp) > NekConstants::kNekMachineEpsilon)
            {
                m_travelWave[j] = tmp;
                m_hasVel[j]     = true;
            }
        }
    }

    // load initial displacement
    std::map<int, LibUtilities::EquationSharedPtr> initDisplacement;
    funcNameElmt = pForce->FirstChildElement("INITIALDISPLACEMENT");
    if (funcNameElmt)
    {
        std::string FuncName = funcNameElmt->GetText();
        ASSERTL0(m_session->DefinesFunction(FuncName),
                 "Function '" + FuncName + "' is not defined in the session.");
        // linear displacement
        std::vector<std::string> dispVar = {"x", "y", "z"};
        for (int i = 0; i < m_spacedim; ++i)
        {
            std::string var = dispVar[i];
            if (m_session->DefinesFunction(FuncName, var))
            {
                m_disp[i] = m_session->GetFunction(FuncName, var)
                                ->Evaluate(0., 0., 0., time);
            }
        }
        // angular displacement
        std::vector<std::string> angleVar = {"theta_x", "theta_y", "theta_z"};
        for (int i = 0; i < 3; ++i)
        {
            std::string var = angleVar[i];
            if (m_session->DefinesFunction(FuncName, var))
            {
                m_disp[i + 3] =
                    m_session->GetFunction(FuncName, var)->Evaluate();
            }
        }
    }
}

void ForcingMovingReferenceFrame::InitBodySolver(const TiXmlElement *pForce,
                                                 const int dim, const int rank,
                                                 const NekDouble time)
{
    int NumDof = dim + 1;
    std::set<int> DirBCs;
    const TiXmlElement *mssgTag;
    std::string mssgStr;
    // read free motion DoFs
    m_DirDoFs.clear();
    for (int i = 0; i < NumDof; ++i)
    {
        m_DirDoFs.insert(i);
    }
    mssgTag = pForce->FirstChildElement("MOTIONPRESCRIBED");
    if (mssgTag)
    {
        std::vector<std::string> values;
        mssgStr = mssgTag->GetText();
        ParseUtils::GenerateVector(mssgStr, values);
        ASSERTL0(values.size() >= NumDof,
                 "MOTIONPRESCRIBED vector should be of size " +
                     std::to_string(NumDof));
        for (int i = 0; i < NumDof; ++i)
        {
            if (EvaluateExpression(values[i]) == 0)
            {
                m_DirDoFs.erase(i);
                if (i < dim)
                {
                    m_hasVel[i] = true;
                }
                else if (i == dim)
                {
                    m_hasOmega[2] = true;
                }
            }
        }
    }
    m_hasFreeMotion = m_DirDoFs.size() < NumDof;
    // read mass matrix
    DNekMatSharedPtr M =
        MemoryManager<DNekMat>::AllocateSharedPtr(NumDof, NumDof, 0.0, eFULL);
    mssgTag = pForce->FirstChildElement("MASS");
    ASSERTL0(m_DirDoFs.size() == NumDof || mssgTag, "Mass matrix is required.");
    if (mssgTag)
    {
        std::vector<std::string> values;
        mssgStr = mssgTag->GetText();
        ParseUtils::GenerateVector(mssgStr, values);
        ASSERTL0(values.size() >= NumDof * NumDof,
                 "Mass matrix should be of size " + std::to_string(NumDof) +
                     "X" + std::to_string(NumDof));
        for (int i = 0; i < NumDof; ++i)
        {
            for (int j = 0; j < NumDof; ++j)
            {
                M->SetValue(i, j, EvaluateExpression(values[i * NumDof + j]));
            }
        }
    }
    // read damping matrix
    DNekMatSharedPtr C =
        MemoryManager<DNekMat>::AllocateSharedPtr(NumDof, NumDof, 0.0, eFULL);
    mssgTag = pForce->FirstChildElement("DAMPING");
    if (mssgTag)
    {
        std::vector<std::string> values;
        mssgStr = mssgTag->GetText();
        ParseUtils::GenerateVector(mssgStr, values);
        ASSERTL0(values.size() >= NumDof * NumDof,
                 "Damping matrix should be of size " + std::to_string(NumDof) +
                     "X" + std::to_string(NumDof));
        for (int i = 0; i < NumDof; ++i)
        {
            for (int j = 0; j < NumDof; ++j)
            {
                C->SetValue(i, j, EvaluateExpression(values[i * NumDof + j]));
            }
        }
    }
    // read rigidity matrix
    DNekMatSharedPtr K =
        MemoryManager<DNekMat>::AllocateSharedPtr(NumDof, NumDof, 0.0, eFULL);
    mssgTag = pForce->FirstChildElement("RIGIDITY");
    if (mssgTag)
    {
        std::vector<std::string> values;
        mssgStr = mssgTag->GetText();
        ParseUtils::GenerateVector(mssgStr, values);
        ASSERTL0(values.size() >= NumDof * NumDof,
                 "Rigidity matrix should be of size " + std::to_string(NumDof) +
                     "X" + std::to_string(NumDof));
        for (int i = 0; i < NumDof; ++i)
        {
            for (int j = 0; j < NumDof; ++j)
            {
                K->SetValue(i, j, EvaluateExpression(values[i * NumDof + j]));
            }
        }
    }
    // read Newmark Beta paramters
    m_timestep      = m_session->GetParameter("TimeStep");
    NekDouble beta  = 0.25;
    NekDouble gamma = 0.75;
    if (m_session->DefinesParameter("NewmarkBeta"))
    {
        beta = m_session->GetParameter("NewmarkBeta");
    }
    if (m_session->DefinesParameter("NewmarkGamma"))
    {
        gamma = m_session->GetParameter("NewmarkGamma");
    }
    m_bodySolver.SetNewmarkBeta(beta, gamma, m_timestep, M, C, K, m_DirDoFs);
    // OutputFile
    string filename;
    mssgTag = pForce->FirstChildElement("OutputFile");
    if (mssgTag)
    {
        filename = mssgTag->GetText();
    }
    else
    {
        filename = m_session->GetSessionName();
    }
    if (!(filename.length() >= 4 &&
          filename.substr(filename.length() - 4) == ".mrf"))
    {
        filename += ".mrf";
    }
    if (rank == 0)
    {
        m_outputStream.open(filename.c_str());
        if (dim == 2)
        {
            m_outputStream
                << "Variables = t, x, ux, ax, y, uy, ay, theta, omega, domega"
                << endl;
        }
        else if (dim == 3)
        {
            m_outputStream << "Variables = t, x, ux, ax, y, uy, ay, z, uz, az, "
                              "theta, omega, domega"
                           << endl;
        }
    }
    // output frequency
    m_index           = 0;
    m_outputFrequency = 1;
    mssgTag           = pForce->FirstChildElement("OutputFrequency");
    if (mssgTag)
    {
        std::vector<std::string> values;
        mssgStr           = mssgTag->GetText();
        m_outputFrequency = round(EvaluateExpression(mssgStr));
    }
    ASSERTL0(m_outputFrequency > 0,
             "OutputFrequency should be greater than zero.");
    // initialize displacement velocity
    m_bodyVel = Array<OneD, Array<OneD, NekDouble>>(3);
    for (size_t i = 0; i < 3; ++i)
    {
        m_bodyVel[i] = Array<OneD, NekDouble>(NumDof, 0.);
    }
    for (int i = 0; i < m_spacedim; ++i)
    {
        m_bodyVel[0][i] = m_disp[i];
    }
    m_bodyVel[0][NumDof - 1] = m_disp[5];
    for (auto it : m_frameVelFunction)
    {
        if (it.first < 3)
        {
            m_bodyVel[1][it.first] = it.second->Evaluate(0., 0., 0., time);
        }
        else if (it.first == 5)
        {
            m_bodyVel[1][NumDof - 1] = it.second->Evaluate(0., 0., 0., time);
        }
    }
}

/**
 * @brief Updates the forcing array with the current required forcing.
 * @param pFields
 * @param time
 */
void ForcingMovingReferenceFrame::Update(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    // compute the velociites whos functions are provided in inertial frame
    for (auto it : m_frameVelFunction)
    {
        if (it.first < 3)
        {
            m_velXYZ[it.first] = it.second->Evaluate(0., 0., 0., time);
        }
        else
        {
            m_omegaXYZ[it.first - 3] = it.second->Evaluate(0., 0., 0., time);
        }
    }
    for (auto it : m_extForceFunction)
    {
        m_extForceXYZ[it.first] = it.second->Evaluate(0., 0., 0., time);
    }
    auto equ = m_equ.lock();
    ASSERTL0(equ, "Weak pointer to the equation system is expired");
    auto FluidEq = std::dynamic_pointer_cast<FluidInterface>(equ);
    if (time > m_startTime)
    {
        Array<OneD, NekDouble> forcebody(6, 0.); // fluid force
        std::map<int, NekDouble> Dirs;           // prescribed Motion
        for (auto it : m_DirDoFs)
        {
            if (it < m_spacedim)
            {
                Dirs[it] = m_velXYZ[it];
            }
            else if (it == m_spacedim)
            {
                Dirs[it] = m_omegaXYZ[2];
            }
        }
        if (m_hasFreeMotion)
        {
            m_aeroforceFilter->GetForces(pFields, NullNekDouble1DArray, time);
        }
        FluidEq->GetAeroForce(forcebody);
        SolveBodyMotion(m_bodyVel, forcebody, Dirs);
    }
    if (m_rank == 0 && m_index % m_outputFrequency == 0)
    {
        m_outputStream << boost::format("%25.19e") % time << " ";
        for (size_t i = 0; i < m_bodyVel[0].size(); ++i)
        {
            m_outputStream << boost::format("%25.19e") % m_bodyVel[0][i] << " "
                           << boost::format("%25.19e") % m_bodyVel[1][i] << " "
                           << boost::format("%25.19e") % m_bodyVel[2][i] << " ";
        }
        m_outputStream << endl;
    }
    // overwirte m_omegaXYZ with u, also update theta
    for (int i = 0; i < m_spacedim; ++i)
    {
        m_velXYZ[i] = m_bodyVel[1][i];
    }
    m_omegaXYZ[2] = m_bodyVel[1][m_spacedim];
    for (int i = 0; i < m_spacedim; ++i)
    {
        m_disp[i] = m_bodyVel[0][i] + m_travelWave[i] * time;
    }
    m_disp[5] = m_bodyVel[0][m_spacedim];

    // include the effect of rotation
    if (m_hasRotation)
    {
        UpdateRotMat();
        m_frame.SetAngle(m_disp + 3);
        m_frame.IneritalToBody(3, m_velxyz, m_velXYZ);
        m_frame.IneritalToBody(3, m_omegaxyz, m_omegaXYZ);
    }
    else
    {
        // for translation only,
        for (int i = 0; i < m_spacedim; ++i)
        {
            m_velxyz[i] = m_velXYZ[i];
        }
    }
    ++m_index;
}

void ForcingMovingReferenceFrame::SolveBodyMotion(
    Array<OneD, Array<OneD, NekDouble>> &bodyVel,
    const Array<OneD, NekDouble> &forcebody, std::map<int, NekDouble> &Dirs)
{
    if (!m_hasRotation || Dirs.find(m_spacedim) != Dirs.end())
    {
        Array<OneD, NekDouble> force(6, 0.), tmp;
        if (Dirs.find(m_spacedim) != Dirs.end())
        {
            Array<OneD, NekDouble> angle(3, 0.);
            angle[2] =
                bodyVel[0][m_spacedim] +
                0.5 * m_timestep * (bodyVel[1][m_spacedim] + Dirs[m_spacedim]);
            m_frame.SetAngle(angle);
            m_frame.BodyToInerital(3, forcebody, force);
            m_frame.BodyToInerital(3, forcebody + 3, tmp = force + 3);
        }
        else
        {
            Vmath::Vcopy(6, forcebody, 1, force, 1);
        }
        for (int i = 0; i < m_spacedim; ++i)
        {
            force[i] += m_extForceXYZ[i];
        }
        force[m_spacedim] = force[5] + m_extForceXYZ[5];
        m_bodySolver.Solve(bodyVel, force, Dirs);
    }
    else
    {
        Array<OneD, Array<OneD, NekDouble>> tmpbodyVel(bodyVel.size());
        for (size_t i = 0; i < bodyVel.size(); ++i)
        {
            tmpbodyVel[i] = Array<OneD, NekDouble>(bodyVel[i].size());
        }
        Array<OneD, NekDouble> angle(3, 0.);
        angle[2] = bodyVel[0][m_spacedim];
        for (int iter = 0; iter < 2; ++iter)
        {
            // copy initial condition
            for (size_t i = 0; i < bodyVel.size(); ++i)
            {
                Vmath::Vcopy(bodyVel[i].size(), bodyVel[i], 1, tmpbodyVel[i],
                             1);
            }
            Array<OneD, NekDouble> force(6, 0.), tmp;
            m_frame.SetAngle(angle);
            m_frame.BodyToInerital(3, forcebody, force);
            m_frame.BodyToInerital(3, forcebody + 3, tmp = force + 3);
            for (int i = 0; i < m_spacedim; ++i)
            {
                force[i] += m_extForceXYZ[i];
            }
            force[m_spacedim] = force[5] + m_extForceXYZ[5];
            m_bodySolver.Solve(tmpbodyVel, force, Dirs);
            angle[2] = tmpbodyVel[0][m_spacedim];
        }
        // copy final results
        for (size_t i = 0; i < bodyVel.size(); ++i)
        {
            Vmath::Vcopy(bodyVel[i].size(), tmpbodyVel[i], 1, bodyVel[i], 1);
        }
    }
}

void ForcingMovingReferenceFrame::UpdateRotMat()
{
    // update the rotation matrix
    NekDouble sn, cs;
    sn = sin(m_disp[5]);
    cs = cos(m_disp[5]);

    m_ProjMatZ(0, 0) = cs;
    m_ProjMatZ(0, 1) = sn;
    m_ProjMatZ(1, 0) = -sn;
    m_ProjMatZ(1, 1) = cs;
    m_ProjMatZ(2, 2) = 1.0;
}

/**
 * @brief Adds the body force, -Omega X u.
 * @param fields
 * @param inarray
 * @param outarray
 * @param time
 */
void ForcingMovingReferenceFrame::v_Apply(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble &time)
{
    // If there is no rotation, body force is zero,
    // nothing needs to be done here.
    if (m_hasRotation)
    {
        int npoints = fields[0]->GetNpoints();
        addRotation(npoints, outarray, -1., inarray, outarray);
    }
    UpdateBoundaryConditions(time);
}

/**
 * @brief Set velocity boundary condition at the next time step
 */
void ForcingMovingReferenceFrame::UpdateBoundaryConditions(NekDouble time)
{
    time += m_timestep;
    // compute the velocities whose functions are provided in inertial frame
    for (auto it : m_frameVelFunction)
    {
        if (it.first < 3)
        {
            m_velXYZ[it.first] = it.second->Evaluate(0., 0., 0., time);
        }
        else
        {
            m_omegaXYZ[it.first - 3] = it.second->Evaluate(0., 0., 0., time);
        }
    }
    for (auto it : m_extForceFunction)
    {
        m_extForceXYZ[it.first] = it.second->Evaluate(0., 0., 0., time);
    }
    std::map<int, NekDouble> Dirs; // prescribed Motion
    for (auto it : m_DirDoFs)
    {
        if (it < m_spacedim)
        {
            Dirs[it] = m_velXYZ[it];
        }
        else if (it == m_spacedim)
        {
            Dirs[it] = m_omegaXYZ[2];
        }
    }

    auto equ = m_equ.lock();
    ASSERTL0(equ, "Weak pointer to the equation system is expired");
    auto FluidEq = std::dynamic_pointer_cast<FluidInterface>(equ);
    Array<OneD, NekDouble> forcebody(6, 0.); // fluid force
    FluidEq->GetAeroForce(forcebody);

    Array<OneD, Array<OneD, NekDouble>> bodyVel(m_bodyVel.size());
    for (size_t i = 0; i < m_bodyVel.size(); ++i)
    {
        bodyVel[i] = Array<OneD, NekDouble>(m_bodyVel[i].size());
    }
    // copy initial condition
    for (size_t i = 0; i < m_bodyVel.size(); ++i)
    {
        Vmath::Vcopy(m_bodyVel[i].size(), m_bodyVel[i], 1, bodyVel[i], 1);
    }
    SolveBodyMotion(bodyVel, forcebody, Dirs);
    // to stablize the velocity bcs, use the old values if extrapolation is
    // involved
    for (size_t i = 0; i < bodyVel[0].size(); ++i)
    {
        if (Dirs.find(i) == Dirs.end())
        {
            for (size_t k = 0; k < bodyVel.size(); ++k)
            {
                bodyVel[k][i] = m_bodyVel[k][i];
            }
        }
    }
    // set displacements at the next time step
    for (int i = 0; i < m_spacedim; ++i)
    {
        m_disp[i] = bodyVel[0][i] + m_travelWave[i] * time;
    }
    m_disp[5] = bodyVel[0][m_spacedim];
    // copy the velocities and rotation matrix for enforcing bc
    // copy linear and angular velocities and accelerations
    Array<OneD, NekDouble> vel(12, 0.0), tmp;
    for (int i = 0; i < m_spacedim; ++i)
    {
        vel[i]     = bodyVel[1][i];
        vel[i + 6] = bodyVel[2][i];
    }
    vel[5]  = bodyVel[1][m_spacedim];
    vel[11] = bodyVel[2][m_spacedim];
    if (m_hasRotation)
    {
        UpdateRotMat();
        m_frame.SetAngle(m_disp + 3);
        m_frame.IneritalToBody(3, vel, vel);
        m_frame.IneritalToBody(3, tmp = vel + 6, vel + 6);
    }

    // to set the boundary condition of the next time step
    // update the frame velocities
    FluidEq->SetMovingFrameVelocities(vel);
    // update the projection matrix
    FluidEq->SetMovingFrameProjectionMat(m_ProjMatZ);
    // update the frame angle (with respect to the inertial frame)
    // this angle is used to update the meta data,
    // on the other hand, for boundary conditions the projection matrix is used
    FluidEq->SetMovingFrameDisp(m_disp);
}

/**
 * @brief outarray = inarray0 + angVelScale Omega x inarray1
 */
void ForcingMovingReferenceFrame::addRotation(
    int nPnts, // number of points
    const Array<OneD, Array<OneD, NekDouble>> &inarray0, NekDouble angVelScale,
    const Array<OneD, Array<OneD, NekDouble>> &inarray1,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    ASSERTL0(&inarray1 != &outarray, "inarray1 and outarray "
                                     "should not be the same.");

    // TODO: In case of having support for all three components of Omega,
    // they should be transformed into the rotating frame first!

    // In case that the inarray0 and outarry are different, to avoid using
    // un-initialized array, copy the array first
    if (&inarray0 != &outarray)
    {
        ASSERTL0(inarray0.size() == outarray.size(),
                 "inarray0 and outarray must have same dimentions");
        for (int i = 0; i < inarray0.size(); ++i)
        {
            Vmath::Vcopy(nPnts, inarray0[i], 1, outarray[i], 1);
        }
    }

    if (m_spacedim >= 2 && m_hasOmega[2])
    {
        NekDouble cp = m_omegaxyz[2] * angVelScale;
        NekDouble cm = -1. * cp;

        Vmath::Svtvp(nPnts, cm, inarray1[1], 1, outarray[0], 1, outarray[0], 1);
        Vmath::Svtvp(nPnts, cp, inarray1[0], 1, outarray[1], 1, outarray[1], 1);
    }

    if (m_spacedim == 3 && m_hasOmega[0])
    {
        NekDouble cp = m_omegaxyz[0] * angVelScale;
        NekDouble cm = -1. * cp;

        Vmath::Svtvp(nPnts, cp, inarray1[1], 1, outarray[2], 1, outarray[2], 1);
        Vmath::Svtvp(nPnts, cm, inarray1[2], 1, outarray[1], 1, outarray[1], 1);
    }

    if (m_spacedim == 3 && m_hasOmega[1])
    {
        NekDouble cp = m_omegaxyz[1] * angVelScale;
        NekDouble cm = -1. * cp;

        Vmath::Svtvp(nPnts, cp, inarray1[2], 1, outarray[0], 1, outarray[0], 1);
        Vmath::Svtvp(nPnts, cm, inarray1[0], 1, outarray[2], 1, outarray[2], 1);
    }
}

void ForcingMovingReferenceFrame::v_PreApply(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble &time)
{
    Update(fields, time);
    int npoints = fields[0]->GetNpoints();
    if (m_isH2d && fields[0]->GetWaveSpace())
    {
        for (int i = 0; i < m_spacedim; ++i)
        {
            if (m_hasVel[i])
            {
                Array<OneD, NekDouble> tmpphys(npoints,
                                               -m_velxyz[i] - m_travelWave[i]);
                Array<OneD, NekDouble> tmpcoef(npoints);
                fields[0]->HomogeneousFwdTrans(npoints, tmpphys, tmpcoef);
                Vmath::Vadd(npoints, tmpcoef, 1, inarray[i], 1, outarray[i], 1);
            }
            else if (&inarray != &outarray)
            {
                Vmath::Vcopy(npoints, inarray[i], 1, outarray[i], 1);
            }
        }
    }
    else
    {
        int npoints0 = npoints;
        if (m_isH1d && fields[0]->GetWaveSpace())
        {
            npoints0 = m_hasPlane0 ? fields[0]->GetPlane(0)->GetNpoints() : 0;
        }
        for (int i = 0; i < m_spacedim; ++i)
        {
            if (m_hasVel[i])
            {
                Vmath::Sadd(npoints0, -m_velxyz[i] - m_travelWave[i],
                            inarray[i], 1, outarray[i], 1);
                if (&inarray != &outarray && npoints != npoints0)
                {
                    Array<OneD, NekDouble> tmp = outarray[i] + npoints0;
                    Vmath::Vcopy(npoints - npoints0, inarray[i] + npoints0, 1,
                                 tmp, 1);
                }
            }
            else if (&inarray != &outarray)
            {
                Vmath::Vcopy(npoints, inarray[i], 1, outarray[i], 1);
            }
        }
        if (m_hasRotation)
        {
            addRotation(npoints0, outarray, -1., m_coords, outarray);
        }
    }
}

void ForcingMovingReferenceFrame::CheckForRestartTime(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields, NekDouble &time)
{
    std::map<std::string, std::string> fieldMetaDataMap;
    if (m_session->DefinesFunction("InitialConditions"))
    {
        for (int i = 0; i < pFields.size(); ++i)
        {
            LibUtilities::FunctionType vType;

            vType = m_session->GetFunctionType("InitialConditions",
                                               m_session->GetVariable(i));

            if (vType == LibUtilities::eFunctionTypeFile)
            {
                std::string filename = m_session->GetFunctionFilename(
                    "InitialConditions", m_session->GetVariable(i));

                fs::path pfilename(filename);

                // redefine path for parallel file which is in directory
                if (fs::is_directory(pfilename))
                {
                    fs::path metafile("Info.xml");
                    fs::path fullpath = pfilename / metafile;
                    filename          = LibUtilities::PortablePath(fullpath);
                }
                LibUtilities::FieldIOSharedPtr fld =
                    LibUtilities::FieldIO::CreateForFile(m_session, filename);
                fld->ImportFieldMetaData(filename, fieldMetaDataMap);

                // check to see if time is defined
                if (fieldMetaDataMap != LibUtilities::NullFieldMetaDataMap)
                {
                    std::string strt = "Time";
                    if (fieldMetaDataMap.find(strt) != fieldMetaDataMap.end())
                    {
                        time = boost::lexical_cast<NekDouble>(
                            fieldMetaDataMap[strt]);
                        break;
                    }
                }
            }
        }
    }
}

void ForcingMovingReferenceFrame::InitialiseFilter(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    const TiXmlElement *pForce)
{
    std::map<std::string, std::string> vParams;
    vParams["OutputFile"]     = ".dummyfilename";
    const TiXmlElement *param = pForce->FirstChildElement("BOUNDARY");
    ASSERTL0(param, "Body surface should be assigned");
    vParams["Boundary"] = param->GetText();
    m_aeroforceFilter   = MemoryManager<FilterAeroForces>::AllocateSharedPtr(
        pSession, m_equ.lock(), vParams);
    m_aeroforceFilter->Initialise(pFields, 0.0);
}

void Newmark_BetaSolver::SetNewmarkBeta(NekDouble beta, NekDouble gamma,
                                        NekDouble dt, DNekMatSharedPtr M,
                                        DNekMatSharedPtr C, DNekMatSharedPtr K,
                                        std::set<int> DirDoFs)
{
    m_coeffs    = Array<OneD, NekDouble>(5, 0.);
    m_coeffs[0] = 1. / (gamma * dt);
    m_coeffs[1] = 1. / gamma - 1.;
    m_coeffs[2] = beta * dt / gamma;
    m_coeffs[3] = dt * (1. - beta / gamma);
    m_coeffs[4] = (0.5 - beta / gamma) * dt * dt;

    m_rows = M->GetRows();
    m_coeffMatrix =
        MemoryManager<DNekMat>::AllocateSharedPtr(m_rows, m_rows, 0.0, eFULL);
    m_inverseMatrix =
        MemoryManager<DNekMat>::AllocateSharedPtr(m_rows, m_rows, 0.0, eFULL);
    for (int i = 0; i < m_rows; ++i)
    {
        for (int j = 0; j < m_rows; ++j)
        {
            NekDouble value = m_coeffs[0] * M->GetValue(i, j) +
                              C->GetValue(i, j) +
                              m_coeffs[2] * K->GetValue(i, j);
            m_coeffMatrix->SetValue(i, j, value);
            if (DirDoFs.find(i) != DirDoFs.end() ||
                DirDoFs.find(j) != DirDoFs.end())
            {
                value = (i == j) ? 1. : 0.;
            }
            m_inverseMatrix->SetValue(i, j, value);
        }
    }
    m_inverseMatrix->Invert();
    m_M = M;
    m_C = C;
    m_K = K;
}

void Newmark_BetaSolver::Solve(Array<OneD, Array<OneD, NekDouble>> u,
                               Array<OneD, NekDouble> force,
                               std::map<int, NekDouble> motionPrescribed)
{
    NekVector<NekDouble> bm(m_rows);
    for (int i = 0; i < m_rows; ++i)
    {
        bm[i] = m_coeffs[0] * u[1][i] + m_coeffs[1] * u[2][i];
    }
    NekVector<NekDouble> fbm(m_rows);
    Multiply(fbm, *m_M, bm);
    NekVector<NekDouble> bk(m_rows);
    for (int i = 0; i < m_rows; ++i)
    {
        bk[i] = u[0][i] + m_coeffs[3] * u[1][i] + m_coeffs[4] * u[2][i];
    }
    NekVector<NekDouble> fbk(m_rows);
    Multiply(fbk, *m_K, bk);
    NekVector<NekDouble> rhs(m_rows);
    for (int i = 0; i < m_rows; ++i)
    {
        rhs[i] = force[i] - fbk[i] + fbm[i];
    }
    // apply Dirichelet DoFs
    for (int i = 0; i < 3; ++i)
    {
        if (motionPrescribed.find(i) != motionPrescribed.end())
        {
            rhs[i] = motionPrescribed[i];
        }
        else
        {
            for (auto it : motionPrescribed)
            {
                rhs[i] -= it.second * m_coeffMatrix->GetValue(it.first, i);
            }
        }
    }
    // solve
    NekVector<NekDouble> b(m_rows);
    Multiply(b, *m_inverseMatrix, rhs);
    for (int i = 0; i < m_rows; ++i)
    {
        u[1][i] = b[i];
        u[0][i] = m_coeffs[2] * u[1][i] + bk[i];
        u[2][i] = m_coeffs[0] * u[1][i] - bm[i];
    }
}

FrameTransform::FrameTransform()
{
    m_matrix = Array<OneD, NekDouble>(2, 0.);
}

void FrameTransform::SetAngle(const Array<OneD, NekDouble> theta)
{
    m_matrix[0] = cos(theta[2]);
    m_matrix[1] = sin(theta[2]);
}
void FrameTransform::BodyToInerital(const int dim,
                                    const Array<OneD, NekDouble> &body,
                                    Array<OneD, NekDouble> &inertial)
{
    NekDouble xi = body[0], eta = body[1];
    inertial[0] = xi * m_matrix[0] - eta * m_matrix[1];
    inertial[1] = eta * m_matrix[0] + xi * m_matrix[1];
    if (body != inertial && dim > 2)
    {
        for (int i = 2; i < dim; ++i)
        {
            inertial[i] = body[i];
        }
    }
}
void FrameTransform::IneritalToBody(const int dim, Array<OneD, NekDouble> &body,
                                    const Array<OneD, NekDouble> &inertial)
{
    NekDouble x = inertial[0], y = inertial[1];
    body[0] = x * m_matrix[0] + y * m_matrix[1];
    body[1] = y * m_matrix[0] - x * m_matrix[1];
    if (body != inertial && dim > 2)
    {
        for (int i = 2; i < dim; ++i)
        {
            body[i] = inertial[i];
        }
    }
}

} // namespace Nektar::SolverUtils
