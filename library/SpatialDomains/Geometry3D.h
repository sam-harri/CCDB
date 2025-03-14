////////////////////////////////////////////////////////////////////////////////
//
//  File: Geometry3D.h
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
//  Description: 3D geometry information.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SPATIALDOMAINS_GEOMETRY3D_H
#define NEKTAR_SPATIALDOMAINS_GEOMETRY3D_H

#include <StdRegions/StdExpansion3D.h>
#include <StdRegions/StdRegions.hpp>

#include <SpatialDomains/Geometry.h>
#include <SpatialDomains/Geometry1D.h>
#include <SpatialDomains/SpatialDomainsDeclspec.h>

namespace Nektar::SpatialDomains
{

class Geometry2D;
class Geometry3D;
typedef std::shared_ptr<Geometry3D> Geometry3DSharedPtr;
typedef std::shared_ptr<Geometry2D> Geometry2DSharedPtr;
typedef std::vector<Geometry3DSharedPtr> Geometry3DVector;
typedef std::vector<Geometry2DSharedPtr> Geometry2DVector;

class PointGeom;
typedef std::shared_ptr<PointGeom> PointGeomSharedPtr;
typedef std::vector<PointGeomSharedPtr> PointGeomVector;

class SegGeom;
typedef std::shared_ptr<SegGeom> SegGeomSharedPtr;
typedef std::vector<SegGeomSharedPtr> SegGeomVector;

/// 3D geometry information
class Geometry3D : public Geometry
{
public:
    SPATIAL_DOMAINS_EXPORT Geometry3D();
    SPATIAL_DOMAINS_EXPORT Geometry3D(const int coordim);
    SPATIAL_DOMAINS_EXPORT ~Geometry3D() override;

    //---------------------------------------
    // Helper functions
    //---------------------------------------
    SPATIAL_DOMAINS_EXPORT static const int kDim = 3;

protected:
    PointGeomVector m_verts;
    SegGeomVector m_edges;
    Geometry2DVector m_faces;
    std::vector<StdRegions::Orientation> m_eorient;
    std::vector<StdRegions::Orientation> m_forient;
    int m_eid;
    bool m_ownverts;

    //---------------------------------------
    // 3D Geometry Methods
    //---------------------------------------

    SPATIAL_DOMAINS_EXPORT NekDouble
    v_GetLocCoords(const Array<OneD, const NekDouble> &coords,
                   Array<OneD, NekDouble> &Lcoords) override;

    void NewtonIterationForLocCoord(const Array<OneD, const NekDouble> &coords,
                                    const Array<OneD, const NekDouble> &ptsx,
                                    const Array<OneD, const NekDouble> &ptsy,
                                    const Array<OneD, const NekDouble> &ptsz,
                                    Array<OneD, NekDouble> &Lcoords,
                                    NekDouble &dist);
    void NewtonIterationForLocCoord(const Array<OneD, const NekDouble> &coords,
                                    Array<OneD, NekDouble> &Lcoords);

    void v_FillGeom() override;
    NekDouble v_GetCoord(const int i,
                         const Array<OneD, const NekDouble> &Lcoord) override;
    void v_CalculateInverseIsoParam() override;
    int v_AllLeftCheck(const Array<OneD, const NekDouble> &gloCoord) override;

    //---------------------------------------
    // Helper functions
    //---------------------------------------
    int v_GetShapeDim() const override;
    int v_GetNumVerts() const override;
    int v_GetNumEdges() const override;
    int v_GetNumFaces() const override;
    PointGeomSharedPtr v_GetVertex(int i) const override;
    Geometry1DSharedPtr v_GetEdge(int i) const override;
    Geometry2DSharedPtr v_GetFace(int i) const override;
    StdRegions::Orientation v_GetEorient(const int i) const override;
    StdRegions::Orientation v_GetForient(const int i) const override;
};

} // namespace Nektar::SpatialDomains

#endif // NEKTAR_SPATIALDOMAINS_GEOMETRY3D_H
