////////////////////////////////////////////////////////////////////////////////
//
//  File: TriGeom.h
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
//  Description:
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SPATIALDOMAINS_TRIGEOM_H
#define NEKTAR_SPATIALDOMAINS_TRIGEOM_H

#include <SpatialDomains/Geometry2D.h>
#include <SpatialDomains/PointGeom.h>
#include <SpatialDomains/SegGeom.h>
#include <SpatialDomains/SpatialDomainsDeclspec.h>
#include <StdRegions/StdRegions.hpp>

#include <SpatialDomains/GeomFactors.h>
#include <StdRegions/StdTriExp.h>

namespace Nektar::SpatialDomains
{

class TriGeom;
class SegGeom;
struct Curve;

typedef std::shared_ptr<Curve> CurveSharedPtr;
typedef std::shared_ptr<SegGeom> SegGeomSharedPtr;
typedef std::shared_ptr<TriGeom> TriGeomSharedPtr;
typedef std::map<int, TriGeomSharedPtr> TriGeomMap;

class TriGeom : public Geometry2D
{
public:
    SPATIAL_DOMAINS_EXPORT TriGeom();
    SPATIAL_DOMAINS_EXPORT TriGeom(const TriGeom &in);
    SPATIAL_DOMAINS_EXPORT TriGeom(
        const int id, const SegGeomSharedPtr edges[],
        const CurveSharedPtr curve = CurveSharedPtr());
    SPATIAL_DOMAINS_EXPORT ~TriGeom() override;

    /// Get the orientation of face1.
    SPATIAL_DOMAINS_EXPORT static const int kNedges = 3;
    SPATIAL_DOMAINS_EXPORT static const int kNverts = 3;

    SPATIAL_DOMAINS_EXPORT static StdRegions::Orientation GetFaceOrientation(
        const TriGeom &face1, const TriGeom &face2, bool doRot, int dir,
        NekDouble angle, NekDouble tol);

    SPATIAL_DOMAINS_EXPORT static StdRegions::Orientation GetFaceOrientation(
        const PointGeomVector &face1, const PointGeomVector &face2, bool doRot,
        int dir, NekDouble angle, NekDouble tol);

protected:
    SPATIAL_DOMAINS_EXPORT NekDouble v_GetCoord(
        const int i, const Array<OneD, const NekDouble> &Lcoord) override;
    SPATIAL_DOMAINS_EXPORT void v_GenGeomFactors() override;
    SPATIAL_DOMAINS_EXPORT void v_FillGeom() override;
    SPATIAL_DOMAINS_EXPORT int v_GetDir(const int faceidx,
                                        const int facedir) const override;
    SPATIAL_DOMAINS_EXPORT void v_Reset(CurveMap &curvedEdges,
                                        CurveMap &curvedFaces) override;
    SPATIAL_DOMAINS_EXPORT void v_Setup() override;

private:
    void SetUpXmap();
};

} // namespace Nektar::SpatialDomains

#endif
