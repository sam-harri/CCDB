////////////////////////////////////////////////////////////////////////////////
//
//  File: Prism.h
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description: Mesh prism object.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKMESH_MESHELEMENTS_PRISM
#define NEKMESH_MESHELEMENTS_PRISM

#include <NekMesh/MeshElements/Element.h>
#include <NekMesh/NekMeshDeclspec.h>

namespace Nektar::NekMesh
{
/**
 * @brief A 3-dimensional five-faced element (2 triangles, 3
 * quadrilaterals).
 */
class Prism : public Element
{
public:
    /// Creates an instance of this class
    static ElementSharedPtr create(ElmtConfig pConf,
                                   std::vector<NodeSharedPtr> pNodeList,
                                   std::vector<int> pTagList)
    {
        return std::make_shared<Prism>(pConf, pNodeList, pTagList);
    }
    /// Element type
    static LibUtilities::ShapeType m_type;

    NEKMESH_EXPORT Prism(ElmtConfig pConf, std::vector<NodeSharedPtr> pNodeList,
                         std::vector<int> pTagList);
    NEKMESH_EXPORT Prism(const Prism &pSrc);
    NEKMESH_EXPORT ~Prism() override
    {
    }

    NEKMESH_EXPORT SpatialDomains::GeometrySharedPtr GetGeom(
        int coordDim) override;
    NEKMESH_EXPORT void GetCurvedNodes(
        std::vector<NodeSharedPtr> &nodeList) const override;
    NEKMESH_EXPORT StdRegions::Orientation GetEdgeOrient(
        int edgeId, EdgeSharedPtr edge) override;
    NEKMESH_EXPORT void MakeOrder(int order,
                                  SpatialDomains::GeometrySharedPtr geom,
                                  LibUtilities::PointsType pType, int coordDim,
                                  int &id, bool justConfig = false) override;

    NEKMESH_EXPORT static unsigned int GetNumNodes(ElmtConfig pConf);
    NEKMESH_EXPORT int GetFaceVertex(int i, int j) override
    {
        return m_faceIds[i][j];
    }

    /**
     * Orientation of prism; unchanged = 0; clockwise = 1;
     * counter-clockwise = 2. This is set by OrientPrism.
     */
    unsigned int m_orientation;

protected:
    void OrientPrism();

private:
    static int m_faceIds[5][4];
};
} // namespace Nektar::NekMesh

#endif
