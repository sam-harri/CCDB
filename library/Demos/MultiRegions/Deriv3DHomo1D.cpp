///////////////////////////////////////////////////////////////////////////////
//
// File: Deriv3DHomo1D.cpp
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

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <MultiRegions/ContField3DHomogeneous1D.h>
#include <SpatialDomains/MeshGraphIO.h>

using namespace std;
using namespace Nektar;

int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr vSession =
        LibUtilities::SessionReader::CreateInstance(argc, argv);

    LibUtilities::CommSharedPtr vComm = vSession->GetComm();

    MultiRegions::ContField3DHomogeneous1DSharedPtr Exp_u, Exp_v, Exp_w;

    StdRegions::ConstFactorMap factors;
    FlagList flags;

    if (argc < 2)
    {
        fprintf(stderr, "Usage: Deriv3DHomo2D meshfile [SysSolnType]   \n");
        exit(1);
    }

    //----------------------------------------------
    // Read in mesh from input file
    SpatialDomains::MeshGraphSharedPtr graph2D =
        SpatialDomains::MeshGraphIO::Read(vSession);
    //----------------------------------------------

    //----------------------------------------------
    // Define Expansion
    int nzpoints;
    NekDouble lz;
    int FFT;

    vSession->LoadParameter("HomModesZ", nzpoints);
    vSession->LoadParameter("LZ", lz);
    vSession->LoadParameter("USEFFT", FFT);

    bool useFFT = false;
    bool deal   = false;
    if (FFT == 1)
    {
        useFFT = true;
    }

    const LibUtilities::PointsKey PkeyZ(nzpoints,
                                        LibUtilities::eFourierEvenlySpaced);
    const LibUtilities::BasisKey BkeyZ(LibUtilities::eFourier, nzpoints, PkeyZ);

    Exp_u = MemoryManager<MultiRegions::ContField3DHomogeneous1D>::
        AllocateSharedPtr(vSession, BkeyZ, lz, useFFT, deal, graph2D,
                          vSession->GetVariable(0));
    Exp_v = MemoryManager<MultiRegions::ContField3DHomogeneous1D>::
        AllocateSharedPtr(vSession, BkeyZ, lz, useFFT, deal, graph2D,
                          vSession->GetVariable(1));
    Exp_w = MemoryManager<MultiRegions::ContField3DHomogeneous1D>::
        AllocateSharedPtr(vSession, BkeyZ, lz, useFFT, deal, graph2D,
                          vSession->GetVariable(2));

    //----------------------------------------------

    //----------------------------------------------
    // Print summary of solution details
    flags.set(eUseGlobal, false);

    const SpatialDomains::ExpansionInfoMap &expansions =
        graph2D->GetExpansionInfo();

    LibUtilities::BasisKey bkey0 =
        expansions.begin()->second->m_basisKeyVector[0];

    cout << "Calculating Derivatives (Homogeneous in z-plane):" << endl;
    cout << "         Lz              : " << lz << endl;
    cout << "         N.modes         : " << bkey0.GetNumModes() << endl;
    cout << "         N.Z homo modes  : " << BkeyZ.GetNumModes() << endl;
    cout << endl;
    //----------------------------------------------

    //----------------------------------------------
    // Set up coordinates of mesh for Forcing function evaluation
    int nq = Exp_u->GetTotPoints();

    Array<OneD, NekDouble> xc0, xc1, xc2;

    xc0 = Array<OneD, NekDouble>(nq, 0.0);
    xc1 = Array<OneD, NekDouble>(nq, 0.0);
    xc2 = Array<OneD, NekDouble>(nq, 0.0);

    Exp_u->GetCoords(xc0, xc1, xc2);
    //----------------------------------------------
    Array<OneD, NekDouble> dudx, dvdy, dwdz;
    Array<OneD, NekDouble> dump;
    dump = Array<OneD, NekDouble>(nq, 0.0);
    dudx = Array<OneD, NekDouble>(nq, 0.0);
    dvdy = Array<OneD, NekDouble>(nq, 0.0);
    dwdz = Array<OneD, NekDouble>(nq, 0.0);
    //----------------------------------------------
    // Define initial fields
    LibUtilities::EquationSharedPtr ffunc_u =
        vSession->GetFunction("InitialCondition", 0);
    LibUtilities::EquationSharedPtr ffunc_v =
        vSession->GetFunction("InitialCondition", 1);
    LibUtilities::EquationSharedPtr ffunc_w =
        vSession->GetFunction("InitialCondition", 2);

    LibUtilities::EquationSharedPtr exac_u =
        vSession->GetFunction("ExactSolution", 0);
    LibUtilities::EquationSharedPtr exac_v =
        vSession->GetFunction("ExactSolution", 1);
    LibUtilities::EquationSharedPtr exac_w =
        vSession->GetFunction("ExactSolution", 2);

    ffunc_u->Evaluate(xc0, xc1, xc2, Exp_u->UpdatePhys());
    ffunc_v->Evaluate(xc0, xc1, xc2, Exp_v->UpdatePhys());
    ffunc_w->Evaluate(xc0, xc1, xc2, Exp_w->UpdatePhys());

    exac_u->Evaluate(xc0, xc1, xc2, dudx);
    exac_v->Evaluate(xc0, xc1, xc2, dvdy);
    exac_w->Evaluate(xc0, xc1, xc2, dwdz);

    //----------------------------------------------
    // Taking derivative and printing the error

    Exp_u->PhysDeriv(Exp_u->GetPhys(), Exp_u->UpdatePhys(), dump, dump);

    cout << "L infinity error (variable dudx): "
         << Exp_u->Linf(Exp_u->GetPhys(), dudx) << endl;
    cout << "L 2 error (variable dudx)       : "
         << Exp_u->L2(Exp_u->GetPhys(), dudx) << endl;

    Exp_v->PhysDeriv(Exp_v->GetPhys(), dump, Exp_v->UpdatePhys(), dump);

    cout << "L infinity error (variable dvdy): "
         << Exp_v->Linf(Exp_v->GetPhys(), dvdy) << endl;
    cout << "L 2 error (variable dvdy)       : "
         << Exp_v->L2(Exp_v->GetPhys(), dvdy) << endl;

    Exp_w->PhysDeriv(Exp_w->GetPhys(), dump, dump, Exp_w->UpdatePhys());

    cout << "L infinity error (variable dwdz): "
         << Exp_w->Linf(Exp_w->GetPhys(), dwdz) << endl;
    cout << "L 2 error (variable dwdz)       : "
         << Exp_w->L2(Exp_w->GetPhys(), dwdz) << endl;

    return 0;
}
