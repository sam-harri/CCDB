///////////////////////////////////////////////////////////////////////////////
//
// File: SteadyAdvectionDiffusionReaction2D.cpp
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
#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <MultiRegions/ContField.h>
#include <SpatialDomains/MeshGraphIO.h>

using namespace std;
using namespace Nektar;

#ifdef TIMING
#include <time.h>
#define Timing(s)                                                              \
    fprintf(stdout, "%s Took %g seconds\n", s,                                 \
            (clock() - st) / (double)CLOCKS_PER_SEC);                          \
    st = clock();
#else
#define Timing(s) /* Nothing */
#endif

int NoCaseStringCompare(const string &s1, const string &s2);

int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr vSession =
        LibUtilities::SessionReader::CreateInstance(argc, argv);

    MultiRegions::ContFieldSharedPtr Exp, Fce;
    int nq, coordim;
    Array<OneD, NekDouble> fce;
    Array<OneD, NekDouble> xc0, xc1, xc2;
    NekDouble lambda;
    NekDouble ax, ay;

    if (argc < 2)
    {
        fprintf(stderr, "Usage: SteadyAdvectionDiffusionReaction2D  meshfile "
                        "[SysSolnType]\n");
        exit(1);
    }

    //----------------------------------------------
    // Read in mesh from input file
    SpatialDomains::MeshGraphSharedPtr graph2D =
        SpatialDomains::MeshGraphIO::Read(vSession);
    //----------------------------------------------

    //----------------------------------------------
    // Get Advection Velocity
    ax = vSession->GetParameter("Advection_x");
    ay = vSession->GetParameter("Advection_y");
    //----------------------------------------------

    //----------------------------------------------
    // Print summary of solution details
    lambda = vSession->GetParameter("Lambda");
    cout << "            Lambda         : " << lambda << endl;
    const SpatialDomains::ExpansionInfoMap &expansions =
        graph2D->GetExpansionInfo();
    LibUtilities::BasisKey bkey0 =
        expansions.begin()->second->m_basisKeyVector[0];
    LibUtilities::BasisKey bkey1 =
        expansions.begin()->second->m_basisKeyVector[1];
    cout << "Solving Steady 2D LinearAdvection :" << endl;
    cout << "            Advection_x    : " << ax << endl;
    cout << "            Advection_y    : " << ay << endl;
    cout << "            Expansion      : ("
         << LibUtilities::BasisTypeMap[bkey0.GetBasisType()] << ","
         << LibUtilities::BasisTypeMap[bkey1.GetBasisType()] << ")" << endl;
    cout << "            No. modes      : " << bkey0.GetNumModes() << endl;
    cout << endl;
    //----------------------------------------------

    //----------------------------------------------
    // Define Expansion
    Exp = MemoryManager<MultiRegions::ContField>::AllocateSharedPtr(
        vSession, graph2D, vSession->GetVariable(0));
    //----------------------------------------------

    Timing("Read files and define exp ..");

    //----------------------------------------------
    // Set up coordinates of mesh for Forcing function evaluation
    coordim = Exp->GetCoordim(0);
    nq      = Exp->GetTotPoints();

    xc0 = Array<OneD, NekDouble>(nq, 0.0);
    xc1 = Array<OneD, NekDouble>(nq, 0.0);
    xc2 = Array<OneD, NekDouble>(nq, 0.0);

    switch (coordim)
    {
        case 1:
            Exp->GetCoords(xc0);
            break;
        case 2:
            Exp->GetCoords(xc0, xc1);
            break;
        case 3:
            Exp->GetCoords(xc0, xc1, xc2);
            break;
    }

    Array<OneD, Array<OneD, NekDouble>> Vel(2);
    Vel[0] = Array<OneD, NekDouble>(nq, ax);
    Vel[1] = Array<OneD, NekDouble>(nq, ay);

    StdRegions::ConstFactorMap factors;
    StdRegions::VarCoeffMap varcoeffs;

    factors[StdRegions::eFactorLambda] = lambda;

    // Set advection velocities
    StdRegions::VarCoeffType varcoefftypes[] = {StdRegions::eVarCoeffVelX,
                                                StdRegions::eVarCoeffVelY};
    for (int i = 0; i < 2; i++)
    {
        varcoeffs[varcoefftypes[i]] = Vel[i];
    }
    //----------------------------------------------

    //----------------------------------------------
    // Define forcing function for first variable defined in file
    fce                                   = Array<OneD, NekDouble>(nq);
    LibUtilities::EquationSharedPtr ffunc = vSession->GetFunction("Forcing", 0);

    ffunc->Evaluate(xc0, xc1, xc2, fce);

    //----------------------------------------------

    //----------------------------------------------
    // Setup expansion containing the  forcing function
    Fce = MemoryManager<MultiRegions::ContField>::AllocateSharedPtr(*Exp);
    Fce->SetPhys(fce);
    //----------------------------------------------
    Timing("Define forcing ..");

    //----------------------------------------------
    // ADR solution taking physical forcing
    Exp->LinearAdvectionDiffusionReactionSolve(
        Fce->GetPhys(), Exp->UpdateCoeffs(), factors, varcoeffs);
    //----------------------------------------------
    Timing("Linear Advection Solve ..");

    //----------------------------------------------
    // Backward Transform Solution to get solved values
    Exp->BwdTrans(Exp->GetCoeffs(), Exp->UpdatePhys());
    // Exp->BwdTrans(Exp->GetContCoeffs(), Exp->UpdatePhys(), true);
    //----------------------------------------------

    //----------------------------------------------
    // See if there is an exact solution, if so
    // evaluate and plot errors
    LibUtilities::EquationSharedPtr ex_sol =
        vSession->GetFunction("ExactSolution", 0);

    if (ex_sol)
    {
        //----------------------------------------------
        // evaluate exact solution
        ex_sol->Evaluate(xc0, xc1, xc2, fce);

        //--------------------------------------------
        // Calculate L_inf error
        Fce->SetPhys(fce);
        Fce->SetPhysState(true);

        cout << "L infinity error: "
             << Exp->Linf(Exp->GetPhys(), Fce->GetPhys()) << endl;
        cout << "L 2 error:        " << Exp->L2(Exp->GetPhys(), Fce->GetPhys())
             << endl;
        //--------------------------------------------
    }
    //----------------------------------------------
    return 0;
}

/**
 * Performs a case-insensitive string comparison (from web).
 * @param   s1          First string to compare.
 * @param   s2          Second string to compare.
 * @returns             0 if the strings match.
 */
int NoCaseStringCompare(const string &s1, const string &s2)
{
    string::const_iterator it1 = s1.begin();
    string::const_iterator it2 = s2.begin();

    // stop when either string's end has been reached
    while ((it1 != s1.end()) && (it2 != s2.end()))
    {
        if (::toupper(*it1) != ::toupper(*it2)) // letters differ?
        {
            // return -1 to indicate smaller than, 1 otherwise
            return (::toupper(*it1) < ::toupper(*it2)) ? -1 : 1;
        }

        // proceed to the next character in each string
        ++it1;
        ++it2;
    }

    size_t size1 = s1.size();
    size_t size2 = s2.size(); // cache lengths

    // return -1,0 or 1 according to strings' lengths
    if (size1 == size2)
    {
        return 0;
    }

    return (size1 < size2) ? -1 : 1;
}
