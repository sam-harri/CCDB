///////////////////////////////////////////////////////////////////////////////
//
// File: Aliasing.cpp
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
#include <SolverUtils/Driver.h>
#include <SpatialDomains/MeshGraphIO.h>

#include <IncNavierStokesSolver/AdvectionTerms/NavierStokesAdvection.h>
#include <IncNavierStokesSolver/EquationSystems/IncNavierStokes.h>

using namespace std;
using namespace Nektar;
using namespace Nektar::SolverUtils;

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        fprintf(stderr, "Usage: ./Aliasing file.xml \n");
        fprintf(stderr, "\t Method will read intiial conditions section of "
                        ".xml file for input \n");
        exit(1);
    }

    LibUtilities::SessionReaderSharedPtr session;
    SpatialDomains::MeshGraphSharedPtr graph;
    string vDriverModule;
    DriverSharedPtr drv;
    try
    {
        // Create session reader.
        session = LibUtilities::SessionReader::CreateInstance(argc, argv);

        // Create MeshGraph.
        graph = SpatialDomains::MeshGraphIO::Read(session);

        // Create driver
        session->LoadSolverInfo("Driver", vDriverModule, "Standard");
        drv = GetDriverFactory().CreateInstance(vDriverModule, session, graph);

        EquationSystemSharedPtr EqSys   = drv->GetEqu()[0];
        IncNavierStokesSharedPtr IncNav = EqSys->as<IncNavierStokes>();

        IncNav->SetInitialConditions(0.0, false);
        Array<OneD, MultiRegions::ExpListSharedPtr> fields =
            IncNav->UpdateFields();

        int i;
        int nConvectiveFields = IncNav->GetNConvectiveFields();
        int nphys             = fields[0]->GetTotPoints();
        Array<OneD, Array<OneD, NekDouble>> VelFields(nConvectiveFields);
        Array<OneD, Array<OneD, NekDouble>> NonLinear(nConvectiveFields);
        Array<OneD, Array<OneD, NekDouble>> NonLinearDealiased(
            nConvectiveFields);

        for (i = 0; i < nConvectiveFields; ++i)
        {
            VelFields[i]          = fields[i]->UpdatePhys();
            NonLinear[i]          = Array<OneD, NekDouble>(nphys);
            NonLinearDealiased[i] = Array<OneD, NekDouble>(nphys);
        }

        std::shared_ptr<NavierStokesAdvection> A =
            std::dynamic_pointer_cast<NavierStokesAdvection>(
                IncNav->GetAdvObject());

        if (!A)
        {
            cout << "Must use non-linear Navier-Stokes advection" << endl;
            exit(-1);
        }

        // calculate non-linear terms without dealiasing
        A->SetSpecHPDealiasing(false);
        A->Advect(nConvectiveFields, fields, VelFields, VelFields, NonLinear,
                  0.0);

        // calculate non-linear terms with dealiasing
        A->SetSpecHPDealiasing(true);
        A->Advect(nConvectiveFields, fields, VelFields, VelFields,
                  NonLinearDealiased, 0.0);

        // Evaulate Difference and put into fields;
        for (i = 0; i < nConvectiveFields; ++i)
        {
            Vmath::Vsub(nphys, NonLinearDealiased[i], 1, NonLinear[i], 1,
                        NonLinear[i], 1);
            fields[i]->FwdTransLocalElmt(NonLinear[i],
                                         fields[i]->UpdateCoeffs());
            // Need to reset varibale name for output
            string name = "NL_Aliasing_" + session->GetVariable(i);
            session->SetVariable(i, name.c_str());
        }

        // Reset session name for output file
        std::string outname = IncNav->GetSessionName();

        outname += "_NonLinear_Aliasing";
        IncNav->ResetSessionName(outname);
        IncNav->Output();
    }
    catch (const std::runtime_error &)
    {
        return 1;
    }
    catch (const std::string &eStr)
    {
        cout << "Error: " << eStr << endl;
    }

    return 0;
}
