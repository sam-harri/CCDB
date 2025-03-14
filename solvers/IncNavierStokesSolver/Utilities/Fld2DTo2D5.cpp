///////////////////////////////////////////////////////////////////////////////
//
// File: Fld2DTo2D5.cpp
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
// Description:
//
///////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <vector>

#include <MultiRegions/ExpList.h>
#include <MultiRegions/ExpList2DHomogeneous1D.h>
#include <MultiRegions/ExpList3DHomogeneous1D.h>
#include <MultiRegions/ExpList3DHomogeneous2D.h>
#include <SpatialDomains/MeshGraphIO.h>

using namespace std;
using namespace Nektar;

int main(int argc, char *argv[])
{
    int i, j, k;

    if (argc != 6)
    {
        fprintf(stderr, "Usage: Fld2DTo2D5 2dmeshfile 2dfieldfile 3dmeshfile "
                        "3dfieldfile outfield\n");
        exit(1);
    }
    string datasave(argv[5]);

    string mesh2d(argv[1]);
    string mesh3d(argv[3]);

    // create 2d session
    LibUtilities::SessionReaderSharedPtr vSession2d =
        LibUtilities::SessionReader::CreateInstance(2, argv);
    std::vector<std::string> filenames;
    filenames.push_back(mesh3d);
    // create 3D session
    LibUtilities::SessionReaderSharedPtr vSession3d =
        LibUtilities::SessionReader::CreateInstance(2, argv, filenames,
                                                    vSession2d->GetComm());

    SpatialDomains::MeshGraphSharedPtr graphShPt2d =
        SpatialDomains::MeshGraphIO::Read(vSession2d);
    SpatialDomains::MeshGraphSharedPtr graphShPt3d =
        SpatialDomains::MeshGraphIO::Read(vSession3d);
    // 2D
    string field2dfile(argv[2]);
    vector<LibUtilities::FieldDefinitionsSharedPtr> field2ddef;
    vector<vector<NekDouble>> field2ddata;
    LibUtilities::Import(field2dfile, field2ddef, field2ddata);
    // 3D
    string field3dfile(argv[4]);
    vector<LibUtilities::FieldDefinitionsSharedPtr> field3ddef;
    vector<vector<NekDouble>> field3ddata;
    LibUtilities::Import(field3dfile, field3ddef, field3ddata);
    vector<vector<NekDouble>> field3ddatanew(field3ddef.size());
    // Set up Expansion information
    vector<vector<LibUtilities::PointsType>> pointstype2d;
    vector<vector<LibUtilities::PointsType>> pointstype3d;
    for (i = 0; i < field2ddef.size(); ++i)
    {
        vector<LibUtilities::PointsType> ptype2d;
        for (j = 0; j < 2; ++j)
        {
            ptype2d.push_back(LibUtilities::ePolyEvenlySpaced);
        }
        pointstype2d.push_back(ptype2d);
    }
    graphShPt2d->SetExpansionInfo(field2ddef, pointstype2d);
    for (i = 0; i < field3ddef.size(); ++i)
    {
        vector<LibUtilities::PointsType> ptype3d;
        for (j = 0; j < 2; ++j)
        {
            ptype3d.push_back(LibUtilities::ePolyEvenlySpaced);
        }
        pointstype3d.push_back(ptype3d);
    }
    graphShPt3d->SetExpansionInfo(field3ddef, pointstype3d);
    bool useFFT     = false;
    bool dealiasing = false;
    // Define Expansion
    // int expdim2d  = graphShPt2d->GetMeshDimension();
    int nfields2d = field2ddef[0]->m_fields.size();
    // int expdim3d  = graphShPt3d->GetMeshDimension();
    int nfields3d = field3ddef[0]->m_fields.size();
    // Gen 2d
    Array<OneD, MultiRegions::ExpListSharedPtr> Exp2d(nfields2d);
    MultiRegions::ExpListSharedPtr Exp2D;
    Exp2D = MemoryManager<MultiRegions::ExpList>::AllocateSharedPtr(
        vSession2d, graphShPt2d);
    Exp2d[0] = Exp2D;
    for (i = 1; i < nfields2d; ++i)
    {
        Exp2d[i] =
            MemoryManager<MultiRegions::ExpList>::AllocateSharedPtr(*Exp2D);
    }
    // Gen 3d
    Array<OneD, MultiRegions::ExpListSharedPtr> Exp3d(nfields3d);
    MultiRegions::ExpList3DHomogeneous1DSharedPtr Exp3DH1;
    // Define Homogeneous expansion
    int nplanes;
    // vSession3d->LoadParameter("HomModesZ",nplanes,field3ddef[0]->m_numModes[2]);
    nplanes = field3ddef[0]->m_numModes[2];
    cout << nplanes << endl;
    // nplanes + 1 points
    const LibUtilities::PointsKey Pkey(nplanes,
                                       LibUtilities::ePolyEvenlySpaced);
    const LibUtilities::BasisKey Bkey(field3ddef[0]->m_basis[2], nplanes, Pkey);
    NekDouble lz = field3ddef[0]->m_homogeneousLengths[0];
    Exp3DH1 =
        MemoryManager<MultiRegions::ExpList3DHomogeneous1D>::AllocateSharedPtr(
            vSession3d, Bkey, lz, useFFT, dealiasing, graphShPt3d);
    Exp3d[0] = Exp3DH1;
    for (j = 1; j < nfields3d; ++j)
    {
        Exp3d[j] = MemoryManager<
            MultiRegions::ExpList3DHomogeneous1D>::AllocateSharedPtr(*Exp3DH1);
    }

    k = 0;
    for (j = 0; j < nfields2d; ++j)
    {
        if (j < nfields2d - 1)
        {
            for (int i = 0; i < field2ddata.size(); ++i)
            {
                Exp2d[j]->ExtractDataToCoeffs(
                    field2ddef[i], field2ddata[i], field2ddef[i]->m_fields[j],
                    Exp3d[j]->GetPlane(k)->UpdateCoeffs());
            }
        }
        if (j == nfields2d - 1)
        {
            for (int i = 0; i < field2ddata.size(); ++i)
            {
                Exp2d[j]->ExtractDataToCoeffs(
                    field2ddef[i], field2ddata[i], field2ddef[i]->m_fields[j],
                    Exp3d[j + 1]->GetPlane(k)->UpdateCoeffs());
            }
        }
    }
    Array<OneD, Array<OneD, NekDouble>> fieldcoeffs(
        vSession3d->GetVariables().size());
    for (j = 0; j < fieldcoeffs.size(); ++j)
    {
        fieldcoeffs[j] = Exp3d[j]->UpdateCoeffs();
        for (int i = 0; i < field3ddef.size(); ++i)
        {
            Exp3d[0]->AppendFieldData(field3ddef[i], field3ddatanew[i],
                                      fieldcoeffs[j]);
        }
    }
    LibUtilities::Write(datasave, field3ddef, field3ddatanew);
    return 0;
}
// Only for 2d to 2d5 Hui Xu 23 Aug 2013
