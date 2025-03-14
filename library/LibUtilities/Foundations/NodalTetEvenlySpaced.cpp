///////////////////////////////////////////////////////////////////////////////
//
// File: NodalTetEvenlySpaced.cpp
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
// Description: 3D Nodal Tetrahedron Evenly Spaced Point Definitions
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Foundations/NodalTetEvenlySpaced.h>
#include <vector>

namespace Nektar::LibUtilities
{
bool NodalTetEvenlySpaced::initPointsManager[] = {
    PointsManager().RegisterCreator(PointsKey(0, eNodalTetEvenlySpaced),
                                    NodalTetEvenlySpaced::Create)};

namespace
{
// construct the geometory and set the coordinate of tetrahedron
// edges and vertices are ordered as anticlockwise

bool isVertex(size_t x, size_t y, size_t z, size_t npts)
{
    return (x == 0 && y == 0 && z == 0) ||
           (x == (npts - 1) && y == 0 && z == 0) ||
           (x == 0 && y == (npts - 1) && z == 0) ||
           (x == 0 && y == 0 && z == (npts - 1));
}

bool isEdge_01([[maybe_unused]] size_t x, size_t y, size_t z,
               [[maybe_unused]] size_t npts)
{ // edge 0
    return y == 0 && z == 0;
}

bool isEdge_12(size_t x, size_t y, size_t z, size_t npts)
{ // edge 1
    return z == 0 && x + y == npts - 1;
}

bool isEdge_20(size_t x, [[maybe_unused]] size_t y, size_t z,
               [[maybe_unused]] size_t npts)
{ // edge 2
    return x == 0 && z == 0;
}

bool isEdge_03(size_t x, size_t y, [[maybe_unused]] size_t z,
               [[maybe_unused]] size_t npts)
{ // edge 3
    return x == 0 && y == 0;
}

bool isEdge_13(size_t x, size_t y, size_t z, size_t npts)
{ // edge 4
    return y == 0 && x + z == npts - 1;
}

bool isEdge_23(size_t x, size_t y, size_t z, size_t npts)
{ // edge 5
    return x == 0 && y + z == npts - 1;
}

bool isEdge(size_t x, size_t y, size_t z, size_t npts)
{
    return isEdge_01(x, y, z, npts) || isEdge_12(x, y, z, npts) ||
           isEdge_20(x, y, z, npts) || isEdge_03(x, y, z, npts) ||
           isEdge_13(x, y, z, npts) || isEdge_23(x, y, z, npts);
}

bool isFace_012([[maybe_unused]] size_t x, [[maybe_unused]] size_t y, size_t z,
                [[maybe_unused]] size_t npts)
{ // bottom face (face 0)
    return z == 0;
}

bool isFace_013([[maybe_unused]] size_t x, size_t y, [[maybe_unused]] size_t z,
                [[maybe_unused]] size_t npts)
{ // face 1
    return y == 0;
}

bool isFace_123(size_t x, size_t y, size_t z, size_t npts)
{ // face 2
    return x + y + z == npts - 1;
}

bool isFace_203(size_t x, [[maybe_unused]] size_t y, [[maybe_unused]] size_t z,
                [[maybe_unused]] size_t npts)
{ // face 3
    return x == 0;
}

bool isFace(size_t x, size_t y, size_t z, size_t npts)
{
    return isFace_012(x, y, z, npts) || isFace_013(x, y, z, npts) ||
           isFace_123(x, y, z, npts) || isFace_203(x, y, z, npts);
}
} // namespace

// Calculate evenly spaced number of points
void NodalTetEvenlySpaced::v_CalculatePoints()
{
    // Allocate the storage for points
    PointsBaseType::v_CalculatePoints();

    // Populate m_points
    size_t npts     = GetNumPoints();
    NekDouble delta = 2.0 / (npts - 1.0);
    for (size_t z = 0, index = 0; z < npts; ++z)
    {
        for (size_t y = 0; y < npts - z; ++y)
        {
            for (size_t x = 0; x < npts - z - y; ++x, ++index)
            {
                NekDouble xi = -1.0 + x * delta;
                NekDouble yi = -1.0 + y * delta;
                NekDouble zi = -1.0 + z * delta;

                m_points[0][index] = xi;
                m_points[1][index] = yi;
                m_points[2][index] = zi;
            }
        }
    }

    NodalPointReorder3d();
    m_util = MemoryManager<NodalUtilTetrahedron>::AllocateSharedPtr(
        npts - 1, m_points[0], m_points[1], m_points[2]);
}

void NodalTetEvenlySpaced::NodalPointReorder3d()
{
    size_t npts = GetNumPoints();
    using std::vector;
    vector<int> vertex;
    vector<int> iEdge_01;             // interior edge 0
    vector<int> iEdge_12;             // interior edge 1
    vector<int> iEdge_20;             // interior edge 2
    vector<int> iEdge_03;             // interior edge 3
    vector<int> iEdge_13;             // interior edge 4
    vector<int> iEdge_23;             // interior edge 5
    vector<int> iFace_012;            // interior face 0
    vector<int> iFace_013;            // interior face 1
    vector<int> iFace_123;            // interior face 2
    vector<int> iFace_203;            // interior face 3
    vector<int> interiorVolumePoints; // interior volume points
    vector<int> map;

    // Build the lattice tetrahedron left to right - bottom to top
    for (size_t z = 0, index = 0; z < npts; ++z)
    {
        for (size_t y = 0; y < npts - z; ++y)
        {
            for (size_t x = 0; x < npts - z - y; ++x, ++index)
            {

                if (isVertex(x, y, z, npts))
                { // vertex

                    vertex.push_back(index);
                }
                else if (isEdge(x, y, z, npts))
                { // interior edge

                    if (isEdge_01(x, y, z, npts))
                    { // interior edge 0

                        iEdge_01.push_back(index);
                    }
                    else if (isEdge_12(x, y, z, npts))
                    { // interior edge 1

                        iEdge_12.push_back(index);
                    }
                    else if (isEdge_20(x, y, z, npts))
                    { // interior edge 2

                        iEdge_20.insert(iEdge_20.begin(), index);
                    }
                    else if (isEdge_03(x, y, z, npts))
                    { // interior edge 3

                        iEdge_03.push_back(index);
                    }
                    else if (isEdge_13(x, y, z, npts))
                    { // interior edge 4

                        iEdge_13.push_back(index);
                    }
                    else if (isEdge_23(x, y, z, npts))
                    { // interior edge 5

                        iEdge_23.push_back(index);
                    }
                }
                else if (isFace(x, y, z, npts))
                { // interior face

                    if (isFace_012(x, y, z, npts))
                    { // interior face 0

                        iFace_012.push_back(index);
                    }
                    else if (isFace_013(x, y, z, npts))
                    { // interior face 1

                        iFace_013.push_back(index);
                    }
                    else if (isFace_123(x, y, z, npts))
                    { // interior face 2

                        iFace_123.push_back(index);
                    }
                    else if (isFace_203(x, y, z, npts))
                    { // interior face 3

                        iFace_203.push_back(index);
                    }
                }
                else
                { // interior volume points

                    interiorVolumePoints.push_back(index);
                }
            }
        }
    }

    // Mapping the vertex, edges, faces, interior volume points using the
    // permutation matrix, so the points are ordered anticlockwise.
    for (size_t n = 0; n < vertex.size(); ++n)
    {

        map.push_back(vertex[n]);
    }

    for (size_t n = 0; n < iEdge_01.size(); ++n)
    {

        map.push_back(iEdge_01[n]);
    }

    for (size_t n = 0; n < iEdge_12.size(); ++n)
    {

        map.push_back(iEdge_12[n]);
    }

    for (size_t n = 0; n < iEdge_20.size(); ++n)
    {

        map.push_back(iEdge_20[n]);
    }

    for (size_t n = 0; n < iEdge_03.size(); ++n)
    {

        map.push_back(iEdge_03[n]);
    }

    for (size_t n = 0; n < iEdge_13.size(); ++n)
    {

        map.push_back(iEdge_13[n]);
    }

    for (size_t n = 0; n < iEdge_23.size(); ++n)
    {

        map.push_back(iEdge_23[n]);
    }

    for (size_t n = 0; n < iFace_012.size(); ++n)
    {

        map.push_back(iFace_012[n]);
    }

    for (size_t n = 0; n < iFace_013.size(); ++n)
    {

        map.push_back(iFace_013[n]);
    }

    for (size_t n = 0; n < iFace_123.size(); ++n)
    {

        map.push_back(iFace_123[n]);
    }

    for (size_t n = 0; n < iFace_203.size(); ++n)
    {

        map.push_back(iFace_203[n]);
    }

    for (size_t n = 0; n < interiorVolumePoints.size(); ++n)
    {

        map.push_back(interiorVolumePoints[n]);
    }

    Array<OneD, NekDouble> points[3];
    points[0] = Array<OneD, NekDouble>(GetTotNumPoints());
    points[1] = Array<OneD, NekDouble>(GetTotNumPoints());
    points[2] = Array<OneD, NekDouble>(GetTotNumPoints());
    for (size_t index = 0; index < map.size(); ++index)
    {

        points[0][index] = m_points[0][index];
        points[1][index] = m_points[1][index];
        points[2][index] = m_points[2][index];
    }

    for (size_t index = 0; index < map.size(); ++index)
    {

        m_points[0][index] = points[0][map[index]];
        m_points[1][index] = points[1][map[index]];
        m_points[2][index] = points[2][map[index]];
    }
}

void NodalTetEvenlySpaced::v_CalculateWeights()
{
    // Allocate storage for points
    PointsBaseType::v_CalculateWeights();

    typedef DataType T;

    // Solve the Vandermonde system of integrals for the weight vector
    NekVector<T> w = m_util->GetWeights();
    m_weights      = Array<OneD, T>(w.GetRows(), w.GetPtr());
}

// ////////////////////////////////////////
//        CalculateInterpMatrix()
void NodalTetEvenlySpaced::CalculateInterpMatrix(
    const Array<OneD, const NekDouble> &xia,
    const Array<OneD, const NekDouble> &yia,
    const Array<OneD, const NekDouble> &zia, Array<OneD, NekDouble> &interp)
{
    Array<OneD, Array<OneD, NekDouble>> xi(3);
    xi[0] = xia;
    xi[1] = yia;
    xi[2] = zia;

    std::shared_ptr<NekMatrix<NekDouble>> mat =
        m_util->GetInterpolationMatrix(xi);
    Vmath::Vcopy(mat->GetRows() * mat->GetColumns(), mat->GetRawPtr(), 1,
                 &interp[0], 1);
}

// ////////////////////////////////////////
//        CalculateDerivMatrix()
void NodalTetEvenlySpaced::v_CalculateDerivMatrix()
{
    // Allocate the derivative matrix.
    PointsBaseType::v_CalculateDerivMatrix();

    m_derivmatrix[0] = m_util->GetDerivMatrix(0);
    m_derivmatrix[1] = m_util->GetDerivMatrix(1);
    m_derivmatrix[2] = m_util->GetDerivMatrix(2);
}

std::shared_ptr<PointsBaseType> NodalTetEvenlySpaced::Create(
    const PointsKey &key)
{
    std::shared_ptr<PointsBaseType> returnval(
        MemoryManager<NodalTetEvenlySpaced>::AllocateSharedPtr(key));

    returnval->Initialize();

    return returnval;
}

} // namespace Nektar::LibUtilities
