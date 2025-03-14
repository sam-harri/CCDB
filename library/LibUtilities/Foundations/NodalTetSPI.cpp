///////////////////////////////////////////////////////////////////////////////
//
// File: NodalTetSPI.cpp
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
// Description: 3D nodal tetrahedral SPI points
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Foundations/NodalTetSPI.h>
#include <LibUtilities/Foundations/NodalTetSPIData.h>

namespace Nektar::LibUtilities
{

bool NodalTetSPI::initPointsManager[] = {PointsManager().RegisterCreator(
    PointsKey(0, eNodalTetSPI), NodalTetSPI::Create)};

void NodalTetSPI::v_CalculatePoints()
{
    // Allocate the storage for points
    size_t numPoints = GetNumPoints();

    for (size_t i = 0; i < 3; i++)
    {
        m_points[i] = Array<OneD, DataType>(NodalTetSPINPTS[numPoints - 2]);
    }

    size_t index = 0;

    // initialize values
    for (size_t i = 0; i < numPoints - 2; ++i)
    {
        index += NodalTetSPINPTS[i];
    }

    for (size_t i = 0; i < NodalTetSPINPTS[numPoints - 2]; i++)
    {
        m_points[0][i] = NodalTetSPIData[index][0];
        m_points[1][i] = NodalTetSPIData[index][1];
        m_points[2][i] = NodalTetSPIData[index][2];
        index++;
    }
}

void NodalTetSPI::v_CalculateWeights()
{
    size_t numPoints = GetNumPoints();

    m_weights = Array<OneD, DataType>(NodalTetSPINPTS[numPoints - 2]);

    size_t index = 0;

    // initialize values
    for (size_t i = 0; i < numPoints - 2; ++i)
    {
        index += NodalTetSPINPTS[i];
    }

    for (size_t i = 0; i < NodalTetSPINPTS[numPoints - 2]; i++)
    {
        m_weights[i] = NodalTetSPIData[index][3];
        index++;
    }
}

void NodalTetSPI::v_CalculateDerivMatrix()
{
}

std::shared_ptr<PointsBaseType> NodalTetSPI::Create(const PointsKey &key)
{
    std::shared_ptr<PointsBaseType> returnval(
        MemoryManager<NodalTetSPI>::AllocateSharedPtr(key));
    returnval->Initialize();
    return returnval;
}

} // namespace Nektar::LibUtilities
