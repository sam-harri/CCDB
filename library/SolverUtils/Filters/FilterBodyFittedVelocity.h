///////////////////////////////////////////////////////////////////////////////
//
// File FilterBodyFittedVelocity.h
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
// Description: Read the body-fitted coordinate system from the file output by
//              FieldConvert module bodyFittedVelocity and compute local max/
//              min/original body-fitted velocity compinents.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_FILTERS_FILTERBODYFITTEDVELOCITY_H
#define NEKTAR_SOLVERUTILS_FILTERS_FILTERBODYFITTEDVELOCITY_H

#include <SolverUtils/Filters/FilterFieldConvert.h>

namespace Nektar::SolverUtils
{
class FilterBodyFittedVelocity : public FilterFieldConvert
{
public:
    friend class MemoryManager<FilterBodyFittedVelocity>;

    /// Creates an instance of this class
    static FilterSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const std::shared_ptr<EquationSystem> &pEquation,
        const std::map<std::string, std::string> &pParams)
    {
        FilterSharedPtr p =
            MemoryManager<FilterBodyFittedVelocity>::AllocateSharedPtr(
                pSession, pEquation, pParams);
        return p;
    }

    /// Name of the class
    static std::string className;

    SOLVER_UTILS_EXPORT FilterBodyFittedVelocity(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const std::shared_ptr<EquationSystem> &pEquation,
        const ParamMap &pParams);
    SOLVER_UTILS_EXPORT ~FilterBodyFittedVelocity() override;

protected:
    // Enumerate types as flags
    enum FilterType
    {
        eOriginal,
        eMax,
        eMin
    };
    enum ProblemType
    {
        eCompressible,
        eIncompressible,
        eOthers
    };

    std::string m_bodyFittedCooriateFile;

    FilterType m_filterType;
    ProblemType m_problemType;
    unsigned int m_spaceDim;
    unsigned int m_nFields;
    unsigned int m_nAddFields;
    unsigned int m_nVars;

    std::vector<std::string> m_bfsVars;
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> m_bfcsDir;

    std::vector<Array<OneD, NekDouble>> m_curFieldsVels_Car;
    std::vector<Array<OneD, NekDouble>> m_curFieldsVels;
    std::vector<Array<OneD, NekDouble>> m_outFieldsVels;

    std::vector<Array<OneD, NekDouble>> m_curFieldsThermalVars;
    std::vector<Array<OneD, NekDouble>> m_outFieldsThermalVars;

    void v_Initialise(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time) override;
    void v_FillVariablesName(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields)
        override;
    void v_ProcessSample(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        std::vector<Array<OneD, NekDouble>> &fieldcoeffs,
        const NekDouble &time) override;
    void v_PrepareOutput(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time) override;
    NekDouble v_GetScale() override;
    std::string v_GetFileSuffix() override
    {
        if (m_filterType == eMax)
        {
            return "_max";
        }
        else if (m_filterType == eMin)
        {
            return "_min";
        }
        else
        {
            return "_original";
        }
    }

private:
    bool m_initialized;
};
} // namespace Nektar::SolverUtils

#endif /* NEKTAR_SOLVERUTILS_FILTERS_FILTERCHECKPOINT_H */
