////////////////////////////////////////////////////////////////////////////////
//
//  File: Geometry1D.h
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
//  Description:  1D geometry information
//
//
////////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_SPATIALDOMAINS_GEOMETRY1D_H
#define NEKTAR_SPATIALDOMAINS_GEOMETRY1D_H

#include <SpatialDomains/Geometry.h>
#include <SpatialDomains/PointGeom.h>
#include <SpatialDomains/SpatialDomainsDeclspec.h>
#include <StdRegions/StdExpansion1D.h>

namespace Nektar::SpatialDomains
{

class Geometry1D;

typedef std::shared_ptr<Geometry1D> Geometry1DSharedPtr;
typedef std::vector<Geometry1DSharedPtr> Geometry1DVector;

/// 1D geometry information
class Geometry1D : public Geometry
{
public:
    SPATIAL_DOMAINS_EXPORT Geometry1D();
    SPATIAL_DOMAINS_EXPORT Geometry1D(const int coordim);
    SPATIAL_DOMAINS_EXPORT ~Geometry1D() override;

    SPATIAL_DOMAINS_EXPORT static const int kDim = 1;

protected:
    int v_GetShapeDim() const override;
    NekDouble v_GetLocCoords(const Array<OneD, const NekDouble> &coords,
                             Array<OneD, NekDouble> &Lcoords) override;
};

} // namespace Nektar::SpatialDomains

#endif // NEKTAR_SPATIALDOMAINS_GEOMETRY1D_H
