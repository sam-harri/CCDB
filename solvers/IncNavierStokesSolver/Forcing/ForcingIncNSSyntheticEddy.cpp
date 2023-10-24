///////////////////////////////////////////////////////////////////////////////
//
// File: ForcingIncNSSyntheticEddy.cpp
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
// Description: Derived base class - Synthetic turbulence forcing
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/Forcing/ForcingIncNSSyntheticEddy.h>
#include <fstream>
#include <iomanip>
#include <ctime>

using namespace std;

namespace Nektar
{
namespace SolverUtils
{

std::string ForcingIncNSSyntheticEddy::className =
    GetForcingFactory().RegisterCreatorFunction(
        "IncNSSyntheticTurbulence", ForcingIncNSSyntheticEddy::create,
        "Inc NS Synthetic Eddy Turbulence Forcing (Generation)");

ForcingIncNSSyntheticEddy::ForcingIncNSSyntheticEddy(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const std::weak_ptr<EquationSystem> &pEquation)
    : Forcing(pSession, pEquation)
{
}

void ForcingIncNSSyntheticEddy::v_InitObject(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    const unsigned int &pNumForcingFields, const TiXmlElement *pForce)
{
    boost::ignore_unused(pNumForcingFields);
    m_session->MatchSolverInfo("Homogeneous", "1D", m_isH1D, false);

    // Check if it is not 2D.
    if (pFields[0]->GetGraph()->GetMeshDimension() < 3)
    {
        if (!m_isH1D)
        {
            NEKERROR(Nektar::ErrorUtil::efatal, "Sythetic eddy method "
                "is only available for three-dimensional simulations");
        }
    }

    // Space dimension
    m_spacedim = pFields[0]->GetGraph()->GetMeshDimension() + (m_isH1D ? 1 : 0);

    // Get gamma parameter
    m_session->LoadParameter("Gamma", m_gamma, 1.4);

    const TiXmlElement *elmtInfTurb;

    // Reynolds Stresses
    elmtInfTurb = pForce->FirstChildElement("ReynoldsStresses");
    ASSERTL0(elmtInfTurb, "Requires ReynoldsStresses tag specifying function "
                          "name which prescribes the reynodls stresses.");

    std::string reyStressesFuncName = elmtInfTurb->GetText();
    ASSERTL0(m_session->DefinesFunction(reyStressesFuncName),
             "Function '" + reyStressesFuncName +
                 "' is not defined in the session.");

    // Reynolds stresses variables. Do not change the order of the variables.
    std::vector<std::string> reyStressesVars = {"r00", "r10", "r20",
                                                "r11", "r21", "r22"};

    for (size_t i = 0; i < reyStressesVars.size(); ++i)
    {
        std::string varStr = reyStressesVars[i];
        if (m_session->DefinesFunction(reyStressesFuncName, varStr))
        {
            m_R[i] = m_session->GetFunction(reyStressesFuncName, varStr);
        }
        else
        {
            NEKERROR(Nektar::ErrorUtil::efatal,
                     "Missing parameter '" + varStr +
                         "' in the FUNCTION NAME = " + reyStressesFuncName +
                         ".");
        }
    }

    // Characteristic length scales
    elmtInfTurb = pForce->FirstChildElement("CharLengthScales");
    ASSERTL0(elmtInfTurb, "Requires CharLengthScales tag specifying function "
                          "name which prescribes the characteristic "
                          "length scales.");

    std::string charLenScalesFuncName = elmtInfTurb->GetText();
    ASSERTL0(m_session->DefinesFunction(charLenScalesFuncName),
             "Function '" + charLenScalesFuncName +
                 "' is not defined in the session.");

    // Characteristic length scale variables
    // Do not change the order of the variables
    std::vector<std::string> lenScalesVars = {"l00", "l10", "l20", "l01", "l11",
                                              "l21", "l02", "l12", "l22"};
    LibUtilities::EquationSharedPtr clsAux;
    m_l = {m_spacedim * m_spacedim, 0.0};

    for (size_t i = 0; i < lenScalesVars.size(); ++i)
    {
        std::string varStr = lenScalesVars[i];
        if (m_session->DefinesFunction(charLenScalesFuncName, varStr))
        {
            clsAux = m_session->GetFunction(charLenScalesFuncName, varStr);
            m_l[i] = (NekDouble)clsAux->Evaluate();
        }
        else
        {
            NEKERROR(Nektar::ErrorUtil::efatal,
                     "Missing parameter '" + varStr +
                         "' in the FUNCTION NAME = " + charLenScalesFuncName +
                         ".");
        }
    }

    // Read box of eddies parameters
    m_rc = {m_spacedim, 0.0};
    // Array<OneD, NekDouble> m_lyz(m_spacedim - 1, 0.0);
    m_lyz       = {m_spacedim - 1, 0.0};
    elmtInfTurb = pForce->FirstChildElement("BoxOfEddies");
    ASSERTL0(elmtInfTurb,
             "Unable to find BoxOfEddies tag. in SyntheticTurbulence forcing");

    if (elmtInfTurb)
    {
        std::stringstream boxStream;
        std::string boxStr = elmtInfTurb->GetText();
        boxStream.str(boxStr);
        size_t countVar = 0;
        for (size_t i = 0; i < (2 * m_spacedim - 1); ++i)
        {
            boxStream >> boxStr;
            if (i < m_spacedim)
            {
                m_rc[i] = boost::lexical_cast<NekDouble>(boxStr);
            }
            else
            {
                m_lyz[i - m_spacedim] = boost::lexical_cast<NekDouble>(boxStr);
            }
            countVar += 1;
        }
        if (countVar != (2 * m_spacedim - 1))
        {
            NEKERROR(Nektar::ErrorUtil::efatal,
                     "Missing parameter in the BoxOfEddies tag");
        }
    }

    // Read sigma
    elmtInfTurb = pForce->FirstChildElement("Sigma");
    ASSERTL0(elmtInfTurb,
             "Unable to find Sigma tag. in SyntheticTurbulence forcing");
    std::string sigmaStr = elmtInfTurb->GetText();
    m_sigma              = boost::lexical_cast<NekDouble>(sigmaStr);

    // Read bulk velocity
    elmtInfTurb = pForce->FirstChildElement("BulkVelocity");
    ASSERTL0(elmtInfTurb,
             "Unable to find BulkVelocity tag. in SyntheticTurbulence forcing");
    std::string bVelStr = elmtInfTurb->GetText();
    m_Ub                = boost::lexical_cast<NekDouble>(bVelStr);

    // Set Cholesky decomposition of the Reynolds Stresses in the domain
    SetCholeskyReyStresses(pFields);
    // Compute reference lengths
    ComputeRefLenghts();
    // Compute Xi maximum
    ComputeXiMax();
    // Set Number of Eddies
    SetNumberOfEddies();
    // Set mask
    SetBoxOfEddiesMask(pFields);
    // Compute initial location of the eddies in the box
    ComputeInitialRandomLocationOfEddies();

    // Seed to generate random positions for the eddies
    srand(time(0));

    // Initialise member from the base class
    m_Forcing = Array<OneD, Array<OneD, NekDouble>>(pFields.size());
    for (int i = 0; i < pFields.size(); ++i)
    {
        m_Forcing[i] = Array<OneD, NekDouble>(pFields[0]->GetTotPoints(), 0.0);
    }

    // Initialise eddies ID vector
    for (size_t n = 0; n < m_N; ++n)
    {
        m_eddiesIDForcing.push_back(n);
    }
}

void ForcingIncNSSyntheticEddy::v_Apply(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble &time)
{
    // Number of Variables
    int nVars = fields.size();
    // Total number of coefficients
    unsigned int nqTot = fields[0]->GetTotPoints();

    // Only apply in the first time step and when an eddy leaves
    // the box
    if (m_calcForcing)
    {
        CalculateForcing(fields);

        for (size_t i = 0; i < (nVars-1); ++i) // Only velocity: nVars - 1
        {
            Vmath::Vadd(nqTot, m_Forcing[i], 1, outarray[i], 1, 
                outarray[i], 1);
        }
        m_calcForcing = false;
    }

    // Check for update
    UpdateEddiesPositions();
}

void ForcingIncNSSyntheticEddy::CalculateForcing(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields)
{
    // Total number of quadrature points
    int nqTot = fields[0]->GetTotPoints();
    // Number of Variables
    int nVars = fields.size();

    // Compute Stochastic Signal
    Array<OneD, Array<OneD, NekDouble>> stochasticSignal;
    stochasticSignal = ComputeStochasticSignal(fields);

    // Compute velocity flsuctuation
    Array<OneD, Array<OneD, NekDouble>> velFluc;
    velFluc = ComputeVelocityFluctuation(fields, stochasticSignal);

    // Compute characteristic convective turbulent time
    Array<OneD, Array<OneD, NekDouble>> convTurbTime;
    convTurbTime = ComputeCharConvTurbTime(fields);

    // Compute smoothing factor
    Array<OneD, Array<OneD, NekDouble>> smoothFac;
    smoothFac = ComputeSmoothingFactor(fields, convTurbTime);

    // Check if eddies left the box
    if(!m_eddiesIDForcing.empty())
    {
        // Clean the m_Forcing member
        for (size_t j = 0; j < nVars; ++j)
        {
            for (int i = 0; i < nqTot; ++i)
            {
                m_Forcing[j][i] = 0.0;
            } 
        }
        // Select the eddies to apply the forcing. Superposition.
        for (auto& n : m_eddiesIDForcing)
        {
            for (size_t i = 0; i < nqTot; ++i)
            {
                if (m_mask[i])
                {
                    //  velocity term
                    for (size_t j = 0; j < m_spacedim; ++j)
                    {
                        m_Forcing[j][i] += 
                            (velFluc[n][j * nqTot + i] * 
                            smoothFac[j][i]) / convTurbTime[j][i];
                    }
                }
            }
        }
        // delete eddies 
        m_eddiesIDForcing.erase(m_eddiesIDForcing.begin(), 
            m_eddiesIDForcing.end());
    }
    else 
    {
        NEKERROR(ErrorUtil::efatal, "Failed: Eddies ID vector is empty.");
    }
}

/**
 * @brief   Compute characteristic convective turbulent time.
 */
Array<OneD, Array<OneD, NekDouble>> ForcingIncNSSyntheticEddy::
    ComputeCharConvTurbTime(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields)
{
    // Total number of quadrature points
    int nqTot = pFields[0]->GetTotPoints();
    // Characteristic convective turbulent time
    Array<OneD, Array<OneD, NekDouble>> convTurbTime(m_spacedim);

    for (size_t k = 0; k < m_spacedim; ++k)
    {
        convTurbTime[k] = Array<OneD, NekDouble>(nqTot, 0.0);
        for (size_t i = 0; i < nqTot; ++i)
        {
            NekDouble convTurbLength = m_xiMaxMin * m_lref[0];
            if ((m_l[k] > m_xiMaxMin * m_lref[0]) && (m_mask[i]))
            {
                convTurbLength = m_l[k];
            }
            convTurbTime[k][i] = convTurbLength / m_Ub;
        }
    }

    return convTurbTime;
}

/**
 * @brief   Compute smoothing factor to avoid strong variations
 *          of the source term across the domain.
 */
Array<OneD, Array<OneD, NekDouble>> ForcingIncNSSyntheticEddy::
    ComputeSmoothingFactor(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        Array<OneD, Array<OneD, NekDouble>> convTurbTime)
{
    // Total number of quadrature points
    int nqTot = pFields[0]->GetTotPoints();
    // Number of elements
    int nElmts = pFields[0]->GetNumElmts();
    // Total number of quadrature points of each element
    int nqe;
    // Characteristic convective turbulent time
    Array<OneD, Array<OneD, NekDouble>> smoothFac(m_spacedim);
    // Counter
    int count = 0;
    // module
    NekDouble mod;
    // Create Array with zeros
    for (size_t i = 0; i < m_spacedim; ++i)
    {
        smoothFac[i] = Array<OneD, NekDouble>(nqTot, 0.0);
    }

    for (size_t e = 0; e < nElmts; ++e)
    {
        nqe = pFields[0]->GetExp(e)->GetTotPoints();

        Array<OneD, NekDouble> coords0(nqe, 0.0);
        Array<OneD, NekDouble> coords1(nqe, 0.0);
        Array<OneD, NekDouble> coords2(nqe, 0.0);

        pFields[0]->GetExp(e)->GetCoords(coords0, coords1, coords2);

        for (size_t i = 0; i < nqe; ++i)
        {
            if (m_mask[count + i])
            {
                // Calculate here all three directions at once.
                mod = (coords0[i] - m_rc[0]) * (coords0[i] - m_rc[0]) +
                      (coords1[i] - m_rc[1]) * (coords1[i] - m_rc[1]) +
                      (coords2[i] - m_rc[2]) * (coords2[i] - m_rc[2]);

                smoothFac[0][count + i] =
                    -0.5 * M_PI * mod * convTurbTime[0][count + i] * m_Ub;
                smoothFac[1][count + i] =
                    -0.5 * M_PI * mod * convTurbTime[1][count + i] * m_Ub;
                smoothFac[2][count + i] =
                    -0.5 * M_PI * mod * convTurbTime[2][count + i] * m_Ub;
            }
        }
        count += nqe;
    }

    return smoothFac;
}

/**
 * @brief Calculate velocity fluctuation for the source term
 *
 * @param stochasticSignal  Stochastic signal
 * @return velFluc         Velocity fluctuation
 */
Array<OneD, Array<OneD, NekDouble>> ForcingIncNSSyntheticEddy::
    ComputeVelocityFluctuation(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        Array<OneD, Array<OneD, NekDouble>> stochasticSignal)
{
    // Total number of quadrature points
    int nqTot = pFields[0]->GetTotPoints();
    // Velocity fluctuation vector
    Array<OneD, Array<OneD, NekDouble>> velFluc(m_N);
    // Control loop for the m_Cholesky
    int l;

    for (size_t n = 0; n < m_N; ++n)
    {   
        velFluc[n] = Array<OneD, NekDouble>(nqTot * m_spacedim, 0.0);

        for (size_t k = 0; k < m_spacedim; ++k)
        {
            for (size_t j = 0; j < k + 1; ++j)
            {
                for (size_t i = 0; i < nqTot; ++i)
                {
                    if (m_mask[i])
                    {
                        l = k + j * (2 * m_spacedim - j - 1) * 0.5;
                        velFluc[n][k * nqTot + i] += m_Cholesky[i][l] 
                            * stochasticSignal[n][j * nqTot + i];
                    }
                }
            }
        }
    }

    return velFluc;
}

/**
 * @brief Compute stochastic signal
*/
Array<OneD, Array<OneD, NekDouble>> ForcingIncNSSyntheticEddy::
    ComputeStochasticSignal(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields)
{
    // Total number of quadrature points
    int nqTot = pFields[0]->GetTotPoints();
    // Number of elements
    int nElmts = pFields[0]->GetNumElmts();
    // Total number of quadrature points of each element
    int nqe;
    // Stochastic Signal vector
    Array<OneD, Array<OneD, NekDouble>> stochasticSignal(m_N);
    // Random numbers: -1 and 1
    Array<OneD, Array<OneD, int>> epsilonSign;

    epsilonSign = GenerateRandomOneOrMinusOne();

    // Calculate the stochastic signal for all eddies.
    for (size_t n = 0; n < m_N; ++n)
    {
        stochasticSignal[n] = Array<OneD, NekDouble>(nqTot * m_spacedim, 0.0);  

        // Evaluate the function at interpolation points for each component
        for (size_t j = 0; j < m_spacedim; ++j)
        {
            // Count the number of quadrature points
            int nqeCount = 0;

            for (size_t e = 0; e < nElmts; ++e)
            {
                nqe = pFields[0]->GetExp(e)->GetTotPoints();

                Array<OneD, NekDouble> coords0(nqe, 0.0);
                Array<OneD, NekDouble> coords1(nqe, 0.0);
                Array<OneD, NekDouble> coords2(nqe, 0.0);

                pFields[0]->GetExp(e)->GetCoords(coords0, coords1, coords2);

                // i: degrees of freedom, j: direction, n: eddies                
                for (size_t i = 0; i < nqe; ++i)
                {
                    if (m_mask[nqeCount + i])
                    {   
                        stochasticSignal[n][j * nqTot + nqeCount + i] = 
                            epsilonSign[j][n] *
                            ComputeGaussian((coords0[i] - m_eddyPos[n][0]) / 
                                m_lref[0], m_xiMax[j * m_spacedim + 0], 
                                    ComputeConstantC(0, j)) *
                            ComputeGaussian((coords1[i] - m_eddyPos[n][1]) / 
                                m_lref[1], m_xiMax[j * m_spacedim + 1], 
                                    ComputeConstantC(1, j)) *
                            ComputeGaussian((coords2[i] - m_eddyPos[n][2]) / 
                                m_lref[2], m_xiMax[j * m_spacedim + 2], 
                                    ComputeConstantC(2, j));                     }
                }  
                nqeCount += nqe;
            }
        }
    }
    
    return stochasticSignal;
}

/**
 * @brief Update the position of the eddies for every time step.
 */
void ForcingIncNSSyntheticEddy::UpdateEddiesPositions()
{
    NekDouble dt = m_session->GetParameter("TimeStep");

    for (size_t n = 0; n < m_N; ++n)
    {
        // Check if any eddy went through the outlet plane (box)
        if ((m_eddyPos[n][0] + m_Ub * dt) < (m_rc[0] + m_lref[0]))
        {
            m_eddyPos[n][0] = m_eddyPos[n][0] + m_Ub * dt;
        }
        else
        {
            // Generate a new one in the inlet plane
            m_eddyPos[n][0] = m_rc[0] - m_lref[0];
            m_eddyPos[n][1] = (m_rc[1] - 0.5 * m_lyz[0]) + 
                    (NekSingle(std::rand()) / NekSingle(RAND_MAX)) * m_lyz[0];
            m_eddyPos[n][2] = (m_rc[2] - 0.5 * m_lyz[1]) + 
                    (NekSingle(std::rand()) / NekSingle(RAND_MAX)) * m_lyz[1];

            m_eddiesIDForcing.push_back(n);
            m_calcForcing = true;
        }
    }
}

/**
 * @brief   Calculate distribution of eddies in the box.
 */
void ForcingIncNSSyntheticEddy::ComputeInitialRandomLocationOfEddies()
{
    m_eddyPos = Array<OneD, Array<OneD, NekDouble>>(m_N);

    for (size_t n = 0; n < m_N; ++n)
    {
        m_eddyPos[n] = Array<OneD, NekDouble>(m_spacedim);
        // Generate randomly eddies inside the box
        m_eddyPos[n][0] = (m_rc[0] - m_lref[0]) + 
                (NekSingle(std::rand()) / NekSingle(RAND_MAX)) * 2 * m_lref[0];
        m_eddyPos[n][1] = (m_rc[1] - 0.5 * m_lyz[0]) + 
                (NekSingle(std::rand()) / NekSingle(RAND_MAX)) * m_lyz[0];
        m_eddyPos[n][2] = (m_rc[2] - 0.5 * m_lyz[1]) + 
                (NekSingle(std::rand()) / NekSingle(RAND_MAX)) * m_lyz[1];
    }
}

/**
 * @brief         Compute standard Gaussian with zero mean
 * @param coord   Coordianate
 * @return        Gaussian value for the coordinate
 */
NekDouble ForcingIncNSSyntheticEddy::ComputeGaussian(NekDouble coord,
                                                      NekDouble xiMaxVal,
                                                      NekDouble constC)
{

    NekDouble epsilon = 1e-6;
    if (abs(coord) <= (xiMaxVal + epsilon))
    {
        return ((1.0 / (m_sigma * sqrt(2.0 * M_PI * constC))) *
                exp(-0.5 * (coord / m_sigma) * (coord / m_sigma)));
    }
    else
    {
        return 0.0;
    }
}

/**
 * @brief   Compute constant C for the gaussian funcion
 */
NekDouble ForcingIncNSSyntheticEddy::ComputeConstantC(int row, int col)
{ 
    NekDouble sizeLenScale = m_xiMax[col * m_spacedim + row];

    // Integration
    NekDouble sum  = 0.0;
    NekDouble step = 0.01;
    NekDouble xi0 = - 1;
    NekDouble xif = 1;
    while (xi0 < xif)
    {
        sum += (ComputeGaussian(xi0 + step, sizeLenScale) *
                ComputeGaussian(xi0 + step, sizeLenScale) +
                ComputeGaussian(xi0, sizeLenScale) * 
                ComputeGaussian(xi0, sizeLenScale));
        xi0 += step;
    }

    return (0.5 * 0.5 * step * sum);
}

/**
 * @brief       Generate random 1 or -1 values to be use to compute
 *              the stochastic signal.
 * @return      ramdom 1 or -1 values
 */
Array<OneD, Array<OneD, int>> ForcingIncNSSyntheticEddy::
    GenerateRandomOneOrMinusOne()
{
    Array<OneD, Array<OneD, int>> epsilonSign(m_spacedim);

    // j: component of the stochastic signal
    // n: eddy
    for (size_t j = 0; j < m_spacedim; ++j)
    {
        epsilonSign[j] = Array<OneD, int>(m_N, 0.0);
        for (size_t n = 0; n < m_N; ++n)
        {
            // Convert to -1 or to 1
            epsilonSign[j][n] =
                ((NekSingle(std::rand()) / NekSingle(RAND_MAX)) <= 0.5) ? -1
                                                                        : 1;
        }
    }

    return epsilonSign;
}

/**
 * @brief   Set box of eddies mask to be use to seprate the
 *          degrees of freedom (coordinates) inside and outside
 *          the box of eddies.
 */
void ForcingIncNSSyntheticEddy::SetBoxOfEddiesMask(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields)
{
    // Total number of quadrature points
    int nqTot = pFields[0]->GetTotPoints();
    // Number of elements
    int nElmts = pFields[0]->GetNumElmts();
    // Total number of quadrature points of each element
    int nqe;
    // Mask
    m_mask = {nqTot, 0}; // 0 for ouside, 1 for inside
    // Counter
    int count = 0;

    for (size_t e = 0; e < nElmts; ++e)
    {
        nqe = pFields[0]->GetExp(e)->GetTotPoints();

        Array<OneD, NekDouble> coords0(nqe, 0.0);
        Array<OneD, NekDouble> coords1(nqe, 0.0);
        Array<OneD, NekDouble> coords2(nqe, 0.0);

        pFields[0]->GetExp(e)->GetCoords(coords0, coords1, coords2);

        for (size_t i = 0; i < nqe; ++i)
        {
            if (InsideBoxOfEddies(coords0[i], coords1[i], coords2[i]))
            {
                // 0 for outside, 1 for inside
                m_mask[count + i] = 1;
            }
        }
        count += nqe;
    }
}

/**
 * @brief Check it point is inside the box of eddies.
 *
 * @param coord0    coordinate in the x-direction
 * @param coord1    coordinate in the y-direction
 * @param coord2    coordinate in the z-direction
 * @return flag     true or false
 */
bool ForcingIncNSSyntheticEddy::InsideBoxOfEddies(NekDouble coord0,
                                                   NekDouble coord1,
                                                   NekDouble coord2)
{
    if ((coord0 > (m_rc[0] - m_lref[0])) && (coord0 < (m_rc[0] + m_lref[0])) &&
        (coord1 > (m_rc[1] - 0.5 * m_lyz[0])) &&
        (coord1 < (m_rc[1] + 0.5 * m_lyz[0])) &&
        (coord2 > (m_rc[2] - 0.5 * m_lyz[1])) &&
        (coord2 < (m_rc[2] + 0.5 * m_lyz[1])))
    {
        return true;
    }

    return false;
}

void ForcingIncNSSyntheticEddy::ComputeRefLenghts()
{
    m_lref    = {m_spacedim, 0.0};
    m_lref[0] = m_l[0];
    m_lref[1] = m_l[1];
    m_lref[2] = m_l[2];

    // The l_{x}^{ref}, l_{y}^{ref} and l_{z}^{ref}
    // are the maximum among the velocity components
    // over the domain.
    for (size_t i = 0; i < m_spacedim; ++i)
    {
        for (size_t j = 1; j < m_spacedim; ++j)
        {
            if (m_l[m_spacedim * j + i] > m_lref[i])
            {
                m_lref[i] = m_l[m_spacedim * j + i];
            }
        }
    }
}

void ForcingIncNSSyntheticEddy::ComputeXiMax()
{
    NekDouble value;
    m_xiMax = {m_spacedim * m_spacedim, 0.0};

    for (size_t i = 0; i < m_spacedim; i++)
    {
        for (size_t j = 0; j < m_spacedim; j++)
        {
            value = (m_l[m_spacedim * j + i] / m_lref[i]);
            if (value > m_xiMaxMin)
            {
                m_xiMax[m_spacedim * j + i] = value;
            }
            else
            {
                m_xiMax[m_spacedim * j + i] = m_xiMaxMin;
            }
        }
    }
}

/**
 * @brief Calculates the Cholesky decomposition of the Reynolds Stresses
 *        in each degree of freedom of the mesh.
 */
void ForcingIncNSSyntheticEddy::SetCholeskyReyStresses(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields)
{
    // Total number of quadrature points
    int nqTot = pFields[0]->GetTotPoints();
    // Total number of quadrature points of each element
    int nqe;
    // Number of elements
    int nElmts = pFields[0]->GetNumElmts();
    // Count the degrees of freedom
    int nqeCount = 0;
    // Size of Cholesky decomposition matrix - aux vector
    Array<OneD, NekDouble> A(m_spacedim * (m_spacedim + 1) * 0.5, 0.0);
    // Cholesky decomposition matrix for the box domain
    m_Cholesky = Array<OneD, Array<OneD, NekDouble>>(nqTot);

    for (size_t e = 0; e < nElmts; ++e)
    {
        nqe = pFields[0]->GetExp(e)->GetTotPoints();

        Array<OneD, NekDouble> coords0(nqe, 0.0);
        Array<OneD, NekDouble> coords1(nqe, 0.0);
        Array<OneD, NekDouble> coords2(nqe, 0.0);

        // Coordinates (for each degree of freedom) for the element k.
        pFields[0]->GetExp(e)->GetCoords(coords0, coords1, coords2);

        // Evaluate Cholesky decomposition
        for (size_t i = 0; i < nqe; ++i)
        {
            int l = 0;
            while (l < m_spacedim * (m_spacedim + 1) / 2)
            {
                // Variable to compute the Cholesky decomposition for each
                // degree of freedom
                A[l] = m_R[l]->Evaluate(coords0[i], coords1[i], coords2[i]);
                l++;
            }
            int info = 0;
            Lapack::Dpptrf('L', m_spacedim, A.get(), info);
            if (info < 0)
            {
                std::string message =
                    "ERROR: The " + std::to_string(-info) +
                    "th parameter had an illegal parameter for dpptrf";
                NEKERROR(ErrorUtil::efatal, message.c_str());
            }
            /*else if (info > 0)
            {
                std::string message = "ERROR: The leading minor of order " +
                                      std::to_string(info) +
                                      " is not positive definite from dpptrf";
                NEKERROR(ErrorUtil::efatal, message.c_str());
            }*/
            m_Cholesky[nqeCount + i] =
                Array<OneD, NekDouble>(m_spacedim * (m_spacedim + 1) / 2, A);
        }
        nqeCount += nqe;
    }
}

void ForcingIncNSSyntheticEddy::SetNumberOfEddies()
{
    m_N = int((m_lyz[0] * m_lyz[1]) /
              (4 * m_lref[m_spacedim - 2] * m_lref[m_spacedim - 1]));
}


} // namespace SolverUtils
} // namespace Nektar