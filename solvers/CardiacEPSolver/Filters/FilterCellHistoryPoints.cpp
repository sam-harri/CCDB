///////////////////////////////////////////////////////////////////////////////
//
// File: FilterCellHistoryPoints.cpp
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
// Description: Outputs values at specific points during time-stepping.
//
///////////////////////////////////////////////////////////////////////////////

#include <CardiacEPSolver/Filters/FilterCellHistoryPoints.h>
#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <boost/format.hpp>
#include <iomanip>

using namespace std;

namespace Nektar
{

std::string FilterCellHistoryPoints::className =
    SolverUtils::GetFilterFactory().RegisterCreatorFunction(
        "CellHistoryPoints", FilterCellHistoryPoints::create);

/**
 *
 */
FilterCellHistoryPoints::FilterCellHistoryPoints(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const std::shared_ptr<SolverUtils::EquationSystem> &pEquation,
    const ParamMap &pParams)
    : FilterHistoryPoints(pSession, pEquation, pParams)
{
}

/**
 *
 */
FilterCellHistoryPoints::~FilterCellHistoryPoints()
{
}

/**
 *
 */
void FilterCellHistoryPoints::v_Update(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    // Only output every m_outputFrequency.
    if ((m_index++) % m_outputFrequency)
    {
        return;
    }

    const size_t numPoints            = m_historyPoints.size();
    const size_t numFields            = m_cell->GetNumCellVariables();
    LibUtilities::CommSharedPtr vComm = pFields[0]->GetComm();
    Array<OneD, NekDouble> data(numPoints * numFields, 0.0);
    Array<OneD, NekDouble> gloCoord(3, 0.0);
    Array<OneD, NekDouble> physvals;
    Array<OneD, NekDouble> locCoord;

    // Pull out data values field by field
    for (size_t j = 0; j < numFields; ++j)
    {
        if (m_isHomogeneous1D)
        {
            // Number of points per plane
            const int nppp = pFields[0]->GetPlane(0)->GetTotPoints();

            for (auto &x : m_historyList)
            {
                locCoord        = std::get<1>(x);
                const int indx  = std::get<2>(x);
                const int expId = std::get<3>(x);

                physvals = m_cell->GetCellSolution(j) + m_outputPlane * nppp +
                           pFields[j]->GetPhys_Offset(expId);

                // interpolate point
                data[indx * numFields + j] =
                    pFields[0]->GetExp(expId)->StdPhysEvaluate(locCoord,
                                                               physvals);
            }
        }
        else
        {
            for (auto &x : m_historyList)
            {
                locCoord        = std::get<1>(x);
                const int indx  = std::get<2>(x);
                const int expId = std::get<3>(x);

                physvals = m_cell->GetCellSolution(j) +
                           pFields[0]->GetPhys_Offset(expId);

                // interpolate point
                data[indx * numFields + j] =
                    pFields[0]->GetExp(expId)->StdPhysEvaluate(locCoord,
                                                               physvals);
            }
        }
    }

    // Exchange history data
    // This could be improved to reduce communication but works for now
    vComm->AllReduce(data, LibUtilities::ReduceSum);
    v_WriteData(vComm->GetRank(), data, numFields, time);
}

void FilterCellHistoryPoints::v_WriteData(const int &rank,
                                          const Array<OneD, NekDouble> &data,
                                          const int &numFields,
                                          const NekDouble &time)
{
    // Only the root process writes out history data
    if (rank > 0)
    {
        return;
    }

    Array<OneD, NekDouble> gloCoord(3, 0.0);
    if (!m_outputOneFile || m_index == 1)
    {
        std::stringstream vTmpFilename;
        std::string vOutputFilename;
        // get the file extension
        std::string ext = fs::path(m_outputFile).extension().string();
        ext             = (ext == "") ? ".his" : ext;
        if (m_outputOneFile)
        {
            vTmpFilename << fs::path(m_outputFile).replace_extension(ext);
        }
        else
        {
            vTmpFilename
                << fs::path(m_outputFile).replace_extension("").string() << "_"
                << m_outputIndex << ext;
        }
        // back up the file if already exists and backup switch is turned on
        vOutputFilename = Filter::SetupOutput(ext, vTmpFilename.str());

        ++m_outputIndex;
        if (m_adaptive)
        {
            m_outputStream.open(vOutputFilename.c_str(), ofstream::app);
        }
        else
        {
            m_outputStream.open(vOutputFilename.c_str());
        }

        if (m_isHomogeneous1D)
        {
            m_outputStream << ") at points:" << endl;
        }
        else
        {
            m_outputStream << ") at points:" << endl;
        }

        for (int i = 0; i < m_historyPoints.size(); ++i)
        {
            gloCoord = m_historyPoints[i];

            m_outputStream << "# " << boost::format("%6.0f") % i;
            m_outputStream << " " << boost::format("%15.9e") % gloCoord[0];
            m_outputStream << " " << boost::format("%15.9e") % gloCoord[1];
            m_outputStream << " " << boost::format("%15.9e") % gloCoord[2];
            m_outputStream << endl;
        }

        if (m_isHomogeneous1D)
        {
            if (m_waveSpace)
            {
                m_outputStream << "# (in Wavespace)" << endl;
            }
        }
    }

    // Write data values point by point
    for (int k = 0; k < m_historyPoints.size(); ++k)
    {
        m_outputStream << boost::format("%15.9e") % time;
        for (int j = 0; j < numFields; ++j)
        {
            m_outputStream << " "
                           << boost::format("%15.9e") % data[k * numFields + j];
        }
        m_outputStream << endl;
    }

    if (!m_outputOneFile)
    {
        m_outputStream.close();
    }
}

} // namespace Nektar
