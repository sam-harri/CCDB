///////////////////////////////////////////////////////////////////////////////
//
// File: Driver.cpp
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
// Description: Base class for Drivers
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Driver.h>
#include <SpatialDomains/MeshGraphIO.h>

using namespace std;

namespace Nektar::SolverUtils
{

std::string Driver::evolutionOperatorLookupIds[6] = {
    LibUtilities::SessionReader::RegisterEnumValue("EvolutionOperator",
                                                   "Nonlinear", eNonlinear),
    LibUtilities::SessionReader::RegisterEnumValue("EvolutionOperator",
                                                   "Direct", eDirect),
    LibUtilities::SessionReader::RegisterEnumValue("EvolutionOperator",
                                                   "Adjoint", eAdjoint),
    LibUtilities::SessionReader::RegisterEnumValue(
        "EvolutionOperator", "TransientGrowth", eTransientGrowth),
    LibUtilities::SessionReader::RegisterEnumValue(
        "EvolutionOperator", "SkewSymmetric", eSkewSymmetric),
    LibUtilities::SessionReader::RegisterEnumValue(
        "EvolutionOperator", "AdaptiveSFD", eAdaptiveSFD)};
std::string Driver::evolutionOperatorDef =
    LibUtilities::SessionReader::RegisterDefaultSolverInfo("EvolutionOperator",
                                                           "Nonlinear");
std::string Driver::driverDefault =
    LibUtilities::SessionReader::RegisterDefaultSolverInfo("Driver",
                                                           "Standard");

/**
 *
 */
DriverFactory &GetDriverFactory()
{
    static DriverFactory instance;
    return instance;
}

/**
 *
 */
Driver::Driver(const LibUtilities::SessionReaderSharedPtr pSession,
               const SpatialDomains::MeshGraphSharedPtr pGraph)
    : m_comm(pSession->GetComm()), m_session(pSession), m_graph(pGraph)
{
}

/**
 *
 */
void Driver::v_InitObject(ostream &out)
{
    try
    {
        // Retrieve the equation system to solve.
        ASSERTL0(m_session->DefinesSolverInfo("EqType"),
                 "EqType SolverInfo tag must be defined.");
        std::string vEquation = m_session->GetSolverInfo("EqType");
        if (m_session->DefinesSolverInfo("SolverType"))
        {
            vEquation = m_session->GetSolverInfo("SolverType");
        }

        // Check such a module exists for this equation.
        ASSERTL0(
            GetEquationSystemFactory().ModuleExists(vEquation),
            "EquationSystem '" + vEquation +
                "' is not defined.\n"
                "Ensure equation name is correct and module is compiled.\n");

        // Retrieve the type of evolution operator to use
        /// @todo At the moment this is Navier-Stokes specific - generalise?
        m_EvolutionOperator =
            m_session->GetSolverInfoAsEnum<EvolutionOperatorType>(
                "EvolutionOperator");

        m_nequ = ((m_EvolutionOperator == eTransientGrowth ||
                   m_EvolutionOperator == eAdaptiveSFD)
                      ? 2
                      : 1);

        m_equ = Array<OneD, EquationSystemSharedPtr>(m_nequ);

        // Set the AdvectiveType tag and create EquationSystem objects.
        switch (m_EvolutionOperator)
        {
            case eNonlinear:
                m_session->SetTag("AdvectiveType", "Convective");
                m_equ[0] = GetEquationSystemFactory().CreateInstance(
                    vEquation, m_session, m_graph);
                break;
            case eDirect:
                m_session->SetTag("AdvectiveType", "Linearised");
                m_equ[0] = GetEquationSystemFactory().CreateInstance(
                    vEquation, m_session, m_graph);
                break;
            case eAdjoint:
                m_session->SetTag("AdvectiveType", "Adjoint");
                m_equ[0] = GetEquationSystemFactory().CreateInstance(
                    vEquation, m_session, m_graph);
                break;
            case eTransientGrowth:
                // forward timestepping
                m_session->SetTag("AdvectiveType", "Linearised");
                m_equ[0] = GetEquationSystemFactory().CreateInstance(
                    vEquation, m_session, m_graph);

                // backward timestepping
                m_session->SetTag("AdvectiveType", "Adjoint");
                m_equ[1] = GetEquationSystemFactory().CreateInstance(
                    vEquation, m_session, m_graph);
                break;
            case eSkewSymmetric:
                m_session->SetTag("AdvectiveType", "SkewSymmetric");
                m_equ[0] = GetEquationSystemFactory().CreateInstance(
                    vEquation, m_session, m_graph);
                break;
            case eAdaptiveSFD:
            {
                // use-opt-file
                bool useOptFile =
                    m_session->DefinesCmdLineArgument("use-opt-file");
                std::string optfilename =
                    useOptFile ? m_session->GetFilenames()[0] : "";

                char *argv[] = {const_cast<char *>("IncNavierStokesSolver"),
                                const_cast<char *>("--use-opt-file"),
                                const_cast<char *>(optfilename.c_str()),
                                nullptr};

                size_t argc = useOptFile ? 3 : 1;

                std::vector<std::string> LinNSFilename;
                if (m_comm->GetSize() == 1)
                {
                    // Use .xml.gz file for serial case
                    LinNSFilename.push_back(m_session->GetSessionName() +
                                            ".xml.gz");
                }
                else
                {
                    // Use previous partition for parallel case
                    LinNSFilename.push_back(m_session->GetSessionName() +
                                            "_xml");
                }
                LinNSFilename.push_back(m_session->GetSessionName() +
                                        "_LinNS.xml");

                session_LinNS = LibUtilities::SessionReader::CreateInstance(
                    argc, argv, LinNSFilename, m_comm);

                // Set graph for LinNs solver.
                SpatialDomains::MeshGraphSharedPtr graph_linns =
                    SpatialDomains::MeshGraphIO::Read(session_LinNS);

                // For running stability analysis
                session_LinNS->SetTag("AdvectiveType", "Linearised");
                m_equ[0] = GetEquationSystemFactory().CreateInstance(
                    vEquation, session_LinNS, graph_linns);

                // For running the SFD method on the nonlinear problem
                m_session->SetTag("AdvectiveType", "Convective");
                m_equ[1] = GetEquationSystemFactory().CreateInstance(
                    vEquation, m_session, m_graph);
            }
            break;
            default:
                ASSERTL0(false, "Unrecognised evolution operator.");
        }
    }
    catch (int e)
    {
        ASSERTL0(e == -1, "No such class class defined.");
        out << "An error occurred during driver initialisation." << endl;
    }
}

} // namespace Nektar::SolverUtils
