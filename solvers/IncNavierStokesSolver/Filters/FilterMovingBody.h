///////////////////////////////////////////////////////////////////////////////
//
// File: FilterMovingBody.h
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

#ifndef NEKTAR_INCNAVIERSTOKES_FILTERS_FILTERMOVINGBODY_H
#define NEKTAR_INCNAVIERSTOKES_FILTERS_FILTERMOVINGBODY_H

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <SolverUtils/Filters/Filter.h>

namespace Nektar
{
class FilterMovingBody;

typedef std::shared_ptr<FilterMovingBody> FilterMovingBodySharedPtr;
typedef std::map<std::string, std::string> FilterParams;
typedef std::pair<std::string, FilterParams> FilterMap;

class FilterMovingBody : public SolverUtils::Filter
{
public:
    friend class MemoryManager<FilterMovingBody>;

    /// Creates an instance of this class
    static SolverUtils::FilterSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const std::shared_ptr<SolverUtils::EquationSystem> &pEquation,
        const ParamMap &pParams)
    {
        SolverUtils::FilterSharedPtr p =
            MemoryManager<FilterMovingBody>::AllocateSharedPtr(
                pSession, pEquation, pParams);
        return p;
    }

    static std::string className;

    FilterMovingBody(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const std::shared_ptr<SolverUtils::EquationSystem> &pEquation,
        const ParamMap &pParams);
    ~FilterMovingBody() override;

    void v_Initialise(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time) override;

    void v_Update(
        [[maybe_unused]] const Array<OneD, const MultiRegions::ExpListSharedPtr>
            &pFields,
        [[maybe_unused]] const NekDouble &time) override
    {
    }

    void UpdateForce(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        Array<OneD, NekDouble> &Aeroforces, const NekDouble &time);

    void UpdateMotion(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        Array<OneD, NekDouble> &MotionVars, const NekDouble &time);

    void v_Finalise(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time) override;

    bool v_IsTimeDependent() override;

private:
    /// ID's of boundary regions where we want the forces
    std::vector<unsigned int> m_boundaryRegionsIdList;
    /// Determines if a given Boundary Region is in
    /// m_boundaryRegionsIdList
    std::vector<bool> m_boundaryRegionIsInList;
    size_t m_index_f;
    size_t m_index_m;
    size_t m_outputFrequency;
    /// plane to take history point from if using a homogeneous1D
    /// expansion
    size_t m_outputPlane;
    bool m_isHomogeneous1D;
    LibUtilities::BasisSharedPtr m_homogeneousBasis;
    std::string m_BoundaryString;
    /// number of planes for homogeneous1D expansion
    size_t m_planes;
    Array<OneD, std::ofstream> m_outputStream;
    std::string m_outputFile_fce;
    std::string m_outputFile_mot;
};

} // namespace Nektar

#endif /* NEKTAR_SOLVERUTILS_FILTERS_FILTERCHECKPOINT_H */
