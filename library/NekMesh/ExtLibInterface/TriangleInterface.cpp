////////////////////////////////////////////////////////////////////////////////
//
//  File: TriangleInterface.cpp
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
//  Description: Interface to triangle mesher
//
////////////////////////////////////////////////////////////////////////////////

#include <NekMesh/ExtLibInterface/TriangleInterface.h>

#include <sstream>

using namespace std;
namespace Nektar::NekMesh
{

void TriangleInterface::Mesh(bool Quality)
{
    SetUp();

    int numPoints = 0;
    int numSeg    = 0;
    for (int i = 0; i < m_boundingloops.size(); i++)
    {
        numSeg += m_boundingloops[i].size();
    }
    numPoints = numSeg + m_stienerpoints.size();

    stringstream ss;
    ss << "3 points required for triangulation, " << numPoints << " provided";

    ASSERTL0(numPoints > 2, ss.str());

    dt.in.numberofpoints          = numPoints;
    dt.in.numberofpointattributes = 0;
    dt.in.pointlist               = new double[dt.in.numberofpoints * 2];

    int pointc = 0;

    for (int i = 0; i < m_boundingloops.size(); i++)
    {
        for (int j = 0; j < m_boundingloops[i].size(); j++, pointc++)
        {
            nodemap[pointc] = m_boundingloops[i][j];

            auto uv = m_boundingloops[i][j]->GetCADSurfInfo(sid);
            dt.in.pointlist[pointc * 2 + 0] = uv[0] * m_str;
            dt.in.pointlist[pointc * 2 + 1] = uv[1];
        }
    }

    for (int i = 0; i < m_stienerpoints.size(); i++, pointc++)
    {
        nodemap[pointc] = m_stienerpoints[i];

        auto uv = m_stienerpoints[i]->GetCADSurfInfo(sid);
        dt.in.pointlist[pointc * 2 + 0] = uv[0] * m_str;
        dt.in.pointlist[pointc * 2 + 1] = uv[1];
    }

    dt.in.numberofsegments = numSeg;
    dt.in.segmentlist      = new int[dt.in.numberofsegments * 2];
    pointc                 = 0;
    for (int i = 0; i < m_boundingloops.size(); i++, pointc++)
    {
        int first = pointc;
        for (int j = 0; j < m_boundingloops[i].size() - 1; j++, pointc++)
        {
            dt.in.segmentlist[pointc * 2 + 0] = pointc;
            dt.in.segmentlist[pointc * 2 + 1] = pointc + 1;
        }
        dt.in.segmentlist[pointc * 2 + 0] = pointc;
        dt.in.segmentlist[pointc * 2 + 1] = first;
    }

    dt.in.numberofregions = 0;
    dt.in.numberofholes   = m_centers.size() - 1;
    dt.in.holelist        = new double[dt.in.numberofholes * 2];

    for (int i = 1; i < m_centers.size(); i++)
    {
        dt.in.holelist[(i - 1) * 2 + 0] = m_centers[i][0] * m_str;
        dt.in.holelist[(i - 1) * 2 + 1] = m_centers[i][1];
    }

    string cmd;
    if (Quality)
    {
        cmd = "pqzQY";
    }
    else if (!Quality)
    {
        cmd = "pzQY";
    }
    char *cstr = new char[cmd.length() + 1];
    strcpy(cstr, cmd.c_str());

    dt.Run(cstr);
}

void TriangleInterface::SetUp()
{
    dt.in.pointlist               = (double *)nullptr;
    dt.in.pointattributelist      = (double *)nullptr;
    dt.in.pointmarkerlist         = (int *)nullptr;
    dt.in.numberofpoints          = 0;
    dt.in.numberofpointattributes = 0;
    //
    dt.in.trianglelist               = (int *)nullptr;
    dt.in.triangleattributelist      = (double *)nullptr;
    dt.in.trianglearealist           = (double *)nullptr;
    dt.in.neighborlist               = (int *)nullptr;
    dt.in.numberoftriangles          = 0;
    dt.in.numberofcorners            = 0;
    dt.in.numberoftriangleattributes = 0;
    //
    dt.in.segmentlist       = (int *)nullptr;
    dt.in.segmentmarkerlist = (int *)nullptr;
    dt.in.numberofsegments  = 0;
    //
    dt.in.holelist      = (double *)nullptr;
    dt.in.numberofholes = 0;
    //
    dt.in.regionlist      = (double *)nullptr;
    dt.in.numberofregions = 0;
    //
    dt.in.edgelist       = (int *)nullptr;
    dt.in.edgemarkerlist = (int *)nullptr;
    dt.in.normlist       = (double *)nullptr;
    dt.in.numberofedges  = 0;
    //
    dt.out.pointlist               = (double *)nullptr;
    dt.out.pointattributelist      = (double *)nullptr;
    dt.out.pointmarkerlist         = (int *)nullptr;
    dt.out.numberofpoints          = 0;
    dt.out.numberofpointattributes = 0;
    //
    dt.out.trianglelist               = (int *)nullptr;
    dt.out.triangleattributelist      = (double *)nullptr;
    dt.out.trianglearealist           = (double *)nullptr;
    dt.out.neighborlist               = (int *)nullptr;
    dt.out.numberoftriangles          = 0;
    dt.out.numberofcorners            = 0;
    dt.out.numberoftriangleattributes = 0;
    //
    dt.out.segmentlist       = (int *)nullptr;
    dt.out.segmentmarkerlist = (int *)nullptr;
    dt.out.numberofsegments  = 0;
    //
    dt.out.holelist      = (double *)nullptr;
    dt.out.numberofholes = 0;
    //
    dt.out.regionlist      = (double *)nullptr;
    dt.out.numberofregions = 0;
    //
    dt.out.edgelist       = (int *)nullptr;
    dt.out.edgemarkerlist = (int *)nullptr;
    dt.out.normlist       = (double *)nullptr;
    dt.out.numberofedges  = 0;
}

void TriangleInterface::Extract(std::vector<std::vector<NodeSharedPtr>> &Connec)
{
    Connec.clear();
    for (int i = 0; i < dt.out.numberoftriangles; i++)
    {
        map<int, NodeSharedPtr>::iterator n1, n2, n3;
        n1 = nodemap.find(dt.out.trianglelist[i * 3 + 0]);
        n2 = nodemap.find(dt.out.trianglelist[i * 3 + 1]);
        n3 = nodemap.find(dt.out.trianglelist[i * 3 + 2]);

        ASSERTL0(n1 != nodemap.end() && n2 != nodemap.end() &&
                     n3 != nodemap.end(),
                 "node index error");

        vector<NodeSharedPtr> tri(3);
        tri[0] = n1->second;
        tri[1] = n2->second;
        tri[2] = n3->second;
        Connec.push_back(tri);
    }
}
} // namespace Nektar::NekMesh
