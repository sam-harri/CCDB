///////////////////////////////////////////////////////////////////////////////
//
// File: NodalPrismSPI.h
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
// Description: Header file of 3D Nodal prism SPI points
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NODALPRISMSPI_H
#define NODALPRISMSPI_H

#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/Foundations/NodalUtil.h>
#include <memory>

namespace Nektar::LibUtilities
{

class NodalPrismSPI : public Points<NekDouble>
{
public:
    ~NodalPrismSPI() override
    {
    }

    LIB_UTILITIES_EXPORT static std::shared_ptr<PointsBaseType> Create(
        const PointsKey &key);

    NodalPrismSPI(const PointsKey &key) : PointsBaseType(key)
    {
    }

private:
    static bool initPointsManager[];

    NodalPrismSPI()                            = delete;
    NodalPrismSPI(const NodalPrismSPI &points) = delete;

    Array<OneD, NekDouble> m_t0, m_t1, m_tw, m_e0, m_ew;
    size_t m_numtri;

    void v_CalculatePoints() final;
    void v_CalculateWeights() final;
    void v_CalculateDerivMatrix() final;
};
} // namespace Nektar::LibUtilities

#endif // NODALTRIELEC_H
