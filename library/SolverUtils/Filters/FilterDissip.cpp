///////////////////////////////////////////////////////////////////////////////
//
// File: FilterDissip.cpp
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
// Description: Output dissipation rate (divided by viscosity).
//
///////////////////////////////////////////////////////////////////////////////

#include <iomanip>

#include <SolverUtils/Filters/FilterDissip.h>
#include <SolverUtils/Filters/FilterInterfaces.hpp>

using namespace std;

namespace Nektar::SolverUtils
{
std::string FilterDissip::className =
    SolverUtils::GetFilterFactory().RegisterCreatorFunction(
        "Dissipation", FilterDissip::create);

FilterDissip::FilterDissip(const LibUtilities::SessionReaderSharedPtr &pSession,
                           const std::shared_ptr<EquationSystem> &pEquation,
                           const ParamMap &pParams)
    : Filter(pSession, pEquation), m_index(-1), m_homogeneous(false), m_planes()
{
    std::string outName;

    // OutputFile
    auto it = pParams.find("OutputFile");
    if (it == pParams.end())
    {
        outName = m_session->GetSessionName();
    }
    else
    {
        ASSERTL0(it->second.length() > 0, "Missing parameter 'OutputFile'.");
        outName = it->second;
    }
    outName += ".csv";

    m_comm = pSession->GetComm();
    if (m_comm->GetRank() == 0)
    {
        m_outFile.open(outName.c_str());
        ASSERTL0(m_outFile.good(), "Unable to open: '" + outName + "'");
        m_outFile.setf(ios::scientific, ios::floatfield);
        m_outFile << "Time, Dissipation rate" << endl;
    }
    pSession->LoadParameter("LZ", m_homogeneousLength, 0.0);

    // OutputFrequency
    it = pParams.find("OutputFrequency");
    ASSERTL0(it != pParams.end(), "Missing parameter 'OutputFrequency'.");
    LibUtilities::Equation equ(m_session->GetInterpreter(), it->second);
    m_outputFrequency = round(equ.Evaluate());
}

FilterDissip::~FilterDissip()
{
}

void FilterDissip::v_Initialise(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    m_index = -1;
    MultiRegions::ExpListSharedPtr areaField;

    ASSERTL0(pFields[0]->GetExpType() != MultiRegions::e1D,
             "1D expansion not supported for dissipation filter");

    ASSERTL0(pFields[0]->GetExpType() != MultiRegions::e2D,
             "2D expansion not supported for dissipation filter");

    ASSERTL0(pFields[0]->GetExpType() != MultiRegions::e3DH2D,
             "Homogeneous 2D expansion not supported for dissipation filter");

    if (pFields[0]->GetExpType() == MultiRegions::e3DH1D)
    {
        m_homogeneous = true;
    }

    // Calculate area/volume of domain.
    if (m_homogeneous)
    {
        m_planes  = pFields[0]->GetZIDs();
        areaField = pFields[0]->GetPlane(0);
    }
    else
    {
        areaField = pFields[0];
    }

    Array<OneD, NekDouble> inarray(areaField->GetNpoints(), 1.0);
    m_area = areaField->Integral(inarray);

    if (m_homogeneous)
    {
        m_area *= m_homogeneousLength;
    }

    // Output values at initial time.
    v_Update(pFields, time);
}

void FilterDissip::v_Update(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    int i, nPoints = pFields[0]->GetNpoints();

    m_index++;

    if (m_index % m_outputFrequency > 0)
    {
        return;
    }

    // Lock equation system pointer
    auto equ = m_equ.lock();
    ASSERTL0(equ, "Weak pointer expired");

    auto fluidEqu = std::dynamic_pointer_cast<FluidInterface>(equ);
    ASSERTL0(fluidEqu, "Dissip filter is incompatible with this solver.");

    // Store physical values in an array
    Array<OneD, Array<OneD, NekDouble>> physfields(pFields.size());
    for (i = 0; i < pFields.size(); ++i)
    {
        physfields[i] = pFields[i]->GetPhys();
    }

    bool waveSpace[3] = {pFields[0]->GetWaveSpace(), pFields[1]->GetWaveSpace(),
                         pFields[2]->GetWaveSpace()};

    // First calculate rate of viscous dissipation field.
    NekDouble Ek = 0.0;
    Array<OneD, NekDouble> tmp(nPoints, 0.0);
    Array<OneD, NekDouble> density;
    Array<OneD, Array<OneD, NekDouble>> u(3);
    for (int i = 0; i < 3; ++i)
    {
        u[i] = Array<OneD, NekDouble>(nPoints);
    }
    fluidEqu->GetVelocity(physfields, u);
    
    Array<OneD, NekDouble> tmp2(nPoints), tmp3(nPoints);
    for (int i = 0; i < 3; ++i)
    {
        for (int j = i; j < 3; j++) 
        {
        	pFields[i]->PhysDeriv(j, u[i], tmp2);
        	pFields[j]->PhysDeriv(i, u[j], tmp3);
        	Vmath::Vadd(nPoints, tmp2, 1, tmp3, 1, tmp2, 1);
        	Vmath::Vmul(nPoints, tmp2, 1, tmp2, 1, tmp2, 1);
        	if (j != i)
        	{
        		Vmath::Smul(nPoints, 2.0, tmp2, 1, tmp2, 1); //off-diagonals get added twice (symmetric matrix)
        	}
        	Vmath::Vadd(nPoints, tmp2, 1, tmp, 1, tmp, 1);
        } 
    } 
    Vmath::Smul(nPoints, 0.5, tmp, 1, tmp, 1); // over 2 so multiply by viscosity for dissipation per volume

    if (!fluidEqu->HasConstantDensity())
    {
        density = Array<OneD, NekDouble>(nPoints);
        fluidEqu->GetDensity(physfields, density);
        Vmath::Vmul(nPoints, density, 1, tmp, 1, tmp, 1);
        pFields[0]->PhysDeriv(0, u[0], tmp3);
        pFields[1]->PhysDeriv(1, u[1], tmp2);
        Vmath::Vadd(nPoints, tmp2, 1, tmp3, 1, tmp3, 1);
        pFields[2]->PhysDeriv(2, u[2], tmp2);
        Vmath::Vadd(nPoints, tmp2, 1, tmp3, 1, tmp3, 1);
        Vmath::Vmul(nPoints, tmp3, 1, tmp3, 1, tmp3, 1);
        Vmath::Svtvp(nPoints, (-2.0/3.0), tmp3, 1, tmp, 1, tmp, 1);
    }

    if (m_homogeneous)
    {
        for (int i = 0; i < 3; ++i)
        {
            pFields[i]->SetWaveSpace(waveSpace[i]);
        }
        pFields[0]->HomogeneousFwdTrans(nPoints, tmp, tmp);
        Ek = pFields[0]->GetPlane(0)->Integral(tmp) * m_homogeneousLength;
    }
    else
    {
        Ek = pFields[0]->Integral(tmp);
    }

    // Ek /= 2.0 * m_area;

    if (m_comm->GetRank() == 0)
    {
        m_outFile << time << "," << Ek << endl;
    }
}

void FilterDissip::v_Finalise(
    [[maybe_unused]] const Array<OneD, const MultiRegions::ExpListSharedPtr>
        &pFields,
    [[maybe_unused]] const NekDouble &time)
{
    m_outFile.close();
}

bool FilterDissip::v_IsTimeDependent()
{
    return true;
}

} // namespace Nektar::SolverUtils
