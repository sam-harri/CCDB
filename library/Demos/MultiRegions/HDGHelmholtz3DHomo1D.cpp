///////////////////////////////////////////////////////////////////////////////
//
// File: HDGHelmholtz3DHomo1D.cpp
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
#include <MultiRegions/DisContField3DHomogeneous1D.h>
#include <SpatialDomains/MeshGraphIO.h>

// #define TIMING
#ifdef TIMING
#include <time.h>
#define Timing(s)                                                              \
    fprintf(stdout, "%s Took %g seconds\n", s,                                 \
            (clock() - st) / (double)CLOCKS_PER_SEC);                          \
    st = clock();
#else
#define Timing(s) /* Nothing */
#endif

using namespace std;
using namespace Nektar;

int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr vSession =
        LibUtilities::SessionReader::CreateInstance(argc, argv);

    LibUtilities::CommSharedPtr vComm = vSession->GetComm();
    string meshfile(argv[1]);

    MultiRegions::DisContField3DHomogeneous1DSharedPtr Exp, Fce;
    MultiRegions::ExpListSharedPtr DerExp1, DerExp2, DerExp3;
    int i, nq;
    Array<OneD, NekDouble> fce;
    Array<OneD, NekDouble> xc0, xc1, xc2;
    StdRegions::ConstFactorMap factors;
    NekDouble lz;

    if (argc < 2)
    {
        fprintf(stderr, "Usage: Helmholtz2D  meshfile\n");
        exit(1);
    }

    LibUtilities::FieldIOSharedPtr fld =
        LibUtilities::FieldIO::CreateDefault(vSession);

    //----------------------------------------------
    // Read in mesh from input file
    SpatialDomains::MeshGraphSharedPtr graph2D =
        SpatialDomains::MeshGraphIO::Read(vSession);
    //----------------------------------------------

    //----------------------------------------------
    // Define Expansion
    int nplanes = vSession->GetParameter("HomModesZ");
    lz          = vSession->GetParameter("LZ");
    bool useFFT = false;
    bool deal   = false;
    const LibUtilities::PointsKey Pkey(nplanes,
                                       LibUtilities::eFourierEvenlySpaced);
    const LibUtilities::BasisKey Bkey(LibUtilities::eFourier, nplanes, Pkey);
    Exp = MemoryManager<MultiRegions::DisContField3DHomogeneous1D>::
        AllocateSharedPtr(vSession, Bkey, lz, useFFT, deal, graph2D,
                          vSession->GetVariable(0));
    //----------------------------------------------
    Timing("Read files and define exp ..");

    //----------------------------------------------
    // Print summary of solution details
    factors[StdRegions::eFactorLambda] = vSession->GetParameter("Lambda");
    factors[StdRegions::eFactorTau]    = 1.0;

    const SpatialDomains::ExpansionInfoMap &expansions =
        graph2D->GetExpansionInfo();
    LibUtilities::BasisKey bkey0 =
        expansions.begin()->second->m_basisKeyVector[0];
    cout << "Solving 3D Helmholtz (Homogeneous in z-direction):" << endl;
    cout << "         Lambda         : " << factors[StdRegions::eFactorLambda]
         << endl;
    cout << "         Lz             : " << lz << endl;
    cout << "         No. modes      : " << bkey0.GetNumModes() << endl;
    cout << "         No. hom. modes : " << Bkey.GetNumModes() << endl;
    cout << endl;
    //----------------------------------------------

    //----------------------------------------------
    // Set up coordinates of mesh for Forcing function evaluation
    nq  = Exp->GetTotPoints();
    xc0 = Array<OneD, NekDouble>(nq, 0.0);
    xc1 = Array<OneD, NekDouble>(nq, 0.0);
    xc2 = Array<OneD, NekDouble>(nq, 0.0);

    Exp->GetCoords(xc0, xc1, xc2);
    //----------------------------------------------

    //----------------------------------------------
    // Define forcing function for first variable defined in file
    fce                                   = Array<OneD, NekDouble>(nq);
    LibUtilities::EquationSharedPtr ffunc = vSession->GetFunction("Forcing", 0);

    ffunc->Evaluate(xc0, xc1, xc2, fce);

    //----------------------------------------------

    //----------------------------------------------
    // Setup expansion containing the  forcing function
    Fce = MemoryManager<
        MultiRegions::DisContField3DHomogeneous1D>::AllocateSharedPtr(*Exp);
    Fce->SetPhys(fce);
    //----------------------------------------------
    Timing("Define forcing ..");

    //----------------------------------------------
    // Helmholtz solution taking physical forcing
    Exp->HelmSolve(Fce->GetPhys(), Exp->UpdateCoeffs(), factors);
    //----------------------------------------------

    Timing("Helmholtz Solve ..");

#ifdef TIMING
    for (i = 0; i < 100; ++i)
    {
        Exp->HelmSolve(Fce->GetPhys(), Exp->UpdateCoeffs(), NullFlagList,
                       factors);
    }

    Timing("100 Helmholtz Solves:... ");
#endif

    //-----------------------------------------------
    // Backward Transform Solution to get solved values at
    Exp->BwdTrans(Exp->GetCoeffs(), Exp->UpdatePhys());
    //-----------------------------------------------
    Timing("Backard Transform ..");

    //-----------------------------------------------
    // Write solution to file
    string out = meshfile.substr(0, meshfile.find_last_of(".")) + ".fld";
    std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef =
        Exp->GetFieldDefinitions();
    std::vector<std::vector<NekDouble>> FieldData(FieldDef.size());

    for (i = 0; i < FieldDef.size(); ++i)
    {
        FieldDef[i]->m_fields.push_back("u");
        Exp->AppendFieldData(FieldDef[i], FieldData[i]);
    }
    fld->Write(out, FieldDef, FieldData);

    //-----------------------------------------------

    //-----------------------------------------------
    // See if there is an exact solution, if so
    // evaluate and plot errors
    LibUtilities::EquationSharedPtr ex_sol =
        vSession->GetFunction("ExactSolution", 0);

    if (ex_sol)
    {
        //----------------------------------------------
        // evaluate exact solution

        ex_sol->Evaluate(xc0, xc1, xc2, fce);

        //----------------------------------------------

        //--------------------------------------------
        // Calculate error
        Fce->SetPhys(fce);
        Fce->SetPhysState(true);

        cout << "L infinity error:  "
             << Exp->Linf(Exp->GetPhys(), Fce->GetPhys()) << endl;
        cout << "L 2 error  :       " << Exp->L2(Exp->GetPhys(), Fce->GetPhys())
             << endl;
        //--------------------------------------------
    }

    Timing("Output ..");
    //----------------------------------------------
    return 0;
}
