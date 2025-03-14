///////////////////////////////////////////////////////////////////////////////
//
// File: FilterElectrogram.h
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
// Description: Outputs virtual electrograms at specific locations.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_FILTERS_FILTERELECTROGRAM_H
#define NEKTAR_SOLVERUTILS_FILTERS_FILTERELECTROGRAM_H

#include <CardiacEPSolver/CellModels/CellModel.h>
#include <SolverUtils/Filters/Filter.h>

namespace Nektar
{

class FilterElectrogram : public SolverUtils::Filter
{
public:
    friend class MemoryManager<FilterElectrogram>;

    /// Creates an instance of this class
    static SolverUtils::FilterSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const std::shared_ptr<SolverUtils::EquationSystem> &pEquation,
        const ParamMap &pParams)
    {
        SolverUtils::FilterSharedPtr p =
            MemoryManager<FilterElectrogram>::AllocateSharedPtr(
                pSession, pEquation, pParams);
        return p;
    }

    /// Name of the class
    static std::string className;

    /// Electrogram filter constructor
    FilterElectrogram(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const std::shared_ptr<SolverUtils::EquationSystem> &pEquation,
        const ParamMap &pParams);

    /// Electrogram filter destructor
    ~FilterElectrogram() override;

protected:
    /// Initialises the electrogram filter and open output file.
    void v_Initialise(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time) override;
    /// Compute extracellular potential at egm points at current time.
    void v_Update(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time) override;
    /// Finalise the electrogram filter and close output file.
    void v_Finalise(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time) override;
    /// Filter is time-dependent and should be called at each time-step.
    bool v_IsTimeDependent() override;

private:
    /// Gradient of the radius from each electrogram point in x-direction
    Array<OneD, Array<OneD, NekDouble>> m_grad_R_x;
    /// Gradient of the radius from each electrogram point in y-direction
    Array<OneD, Array<OneD, NekDouble>> m_grad_R_y;
    /// Gradient of the radius from each electrogram point in z-direction
    Array<OneD, Array<OneD, NekDouble>> m_grad_R_z;
    /// List of electrogram points
    SpatialDomains::PointGeomVector m_electrogramPoints;
    /// Counts number of calls to update (number of timesteps)
    unsigned int m_index;
    /// Number of timesteps between outputs
    unsigned int m_outputFrequency;
    /// Filename to output electrogram data to
    std::string m_outputFile;
    /// Output file stream for electrogram data
    std::ofstream m_outputStream;
    /// Point coordinate input string
    std::stringstream m_electrogramStream;
};
} // namespace Nektar

#endif
