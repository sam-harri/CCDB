///////////////////////////////////////////////////////////////////////////////
//
// File: SkewSymmetricAdvection.h
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
// Description: Skew-Symmetric non linear convective term for Navier-Stokes
// equations
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_NAVIERSTOKESADVECTION_H
#define NEKTAR_SOLVERS_NAVIERSTOKESADVECTION_H

#include <SolverUtils/Advection/Advection.h>

namespace Nektar
{

class SkewSymmetricAdvection : public SolverUtils::Advection

{
public:
    friend class MemoryManager<SkewSymmetricAdvection>;

    /// Creates an instance of this class
    static SolverUtils::AdvectionSharedPtr create(std::string)
    {
        return MemoryManager<SkewSymmetricAdvection>::AllocateSharedPtr();
    }

    /// Name of class
    static std::string className;
    static std::string className2;

protected:
    SkewSymmetricAdvection();

    ~SkewSymmetricAdvection() override;

    void v_InitObject(
        LibUtilities::SessionReaderSharedPtr pSession,
        Array<OneD, MultiRegions::ExpListSharedPtr> pFields) override;

    void v_Advect(const int nConvectiveFields,
                  const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                  const Array<OneD, Array<OneD, NekDouble>> &advVel,
                  const Array<OneD, Array<OneD, NekDouble>> &inarray,
                  Array<OneD, Array<OneD, NekDouble>> &outarray,
                  const NekDouble &time,
                  const Array<OneD, Array<OneD, NekDouble>> &pFwd =
                      NullNekDoubleArrayOfArray,
                  const Array<OneD, Array<OneD, NekDouble>> &pBwd =
                      NullNekDoubleArrayOfArray) override;

private:
    bool m_homogen_dealiasing;
    bool m_SingleMode;
    bool m_HalfMode;
};

} // namespace Nektar

#endif // NEKTAR_SOLVERS_INCNAVIERSTOKES_H
