////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessInterpPtsToPts.cpp
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
//  Description: Interpolate  field to a series of specified points.
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
using namespace std;

#include <boost/lexical_cast.hpp>

#include <FieldUtils/Interpolator.h>
#include <LibUtilities/BasicUtils/CsvIO.h>
#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <LibUtilities/BasicUtils/Progressbar.hpp>
#include <LibUtilities/BasicUtils/PtsIO.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include "ProcessInterpPtsToPts.h"

namespace Nektar::FieldUtils
{

ModuleKey ProcessInterpPtsToPts::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "interpptstopts"),
        ProcessInterpPtsToPts::create,
        "Interpolates a set of points to another, requires fromfld and "
        "fromxml to be defined, a line, plane or block of points can be "
        "defined");

ProcessInterpPtsToPts::ProcessInterpPtsToPts(FieldSharedPtr f)
    : ProcessModule(f)
{
    m_config["topts"] =
        ConfigOption(false, "NotSet", "Pts file to which interpolate field");
    m_config["line"] = ConfigOption(false, "NotSet",
                                    "Specify a line of N points using "
                                    "line=N,x0,y0,z0,z1,y1,z1");
    m_config["plane"] =
        ConfigOption(false, "NotSet",
                     "Specify a plane of N1 x N2 points using "
                     "plane=N1,N2,x0,y0,z0,z1,y1,z1,x2,y2,z2,x3,y3,z3");
    m_config["box"] =
        ConfigOption(false, "NotSet",
                     "Specify a rectangular box of N1 x N2 x N3 points "
                     "using a box of points limited by box="
                     "N1,N2,N3,xmin,xmax,ymin,ymax,zmin,zmax");

    m_config["clamptolowervalue"] =
        ConfigOption(false, "-10000000", "Lower bound for interpolation value");
    m_config["clamptouppervalue"] =
        ConfigOption(false, "10000000", "Upper bound for interpolation value");
    m_config["defaultvalue"] =
        ConfigOption(false, "0", "Default value if point is outside domain");

    m_config["cp"] =
        ConfigOption(false, "NotSet",
                     "Parameters p0 and q to determine pressure coefficients");
}

ProcessInterpPtsToPts::~ProcessInterpPtsToPts()
{
}

void ProcessInterpPtsToPts::v_Process(po::variables_map &vm)
{
    ASSERTL0(m_f->m_fieldPts != LibUtilities::NullPtsField,
             "Should have a PtsField for ProcessInterpPtsToPts.");
    ASSERTL0(m_f->m_comm->GetSpaceComm()->GetSize() == 1,
             "ProcessInterpPtsToPts not implemented in parallel.");

    // Move m_f->m_fieldPts
    LibUtilities::PtsFieldSharedPtr oldPts = m_f->m_fieldPts;
    m_f->m_fieldPts                        = LibUtilities::NullPtsField;

    // Create new fieldPts
    CreateFieldPts(vm);

    int nfields = m_f->m_variables.size();
    for (int j = 0; j < nfields; ++j)
    {
        Array<OneD, NekDouble> newPts(m_f->m_fieldPts->GetNpoints());
        m_f->m_fieldPts->AddField(newPts, m_f->m_variables[j]);
    }

    NekDouble clamp_low = m_config["clamptolowervalue"].as<NekDouble>();
    NekDouble clamp_up  = m_config["clamptouppervalue"].as<NekDouble>();
    NekDouble def_value = m_config["defaultvalue"].as<NekDouble>();

    InterpolatePtsToPts(oldPts, m_f->m_fieldPts, clamp_low, clamp_up,
                        def_value);

    if (!boost::iequals(m_config["cp"].as<string>(), "NotSet"))
    {
        calcCp0();
    }
}

void ProcessInterpPtsToPts::CreateFieldPts(
    [[maybe_unused]] po::variables_map &vm)
{
    int rank   = m_f->m_comm->GetSpaceComm()->GetRank();
    int nprocs = m_f->m_comm->GetSpaceComm()->GetSize();
    // Check for command line point specification
    if (m_config["topts"].as<string>().compare("NotSet") != 0)
    {
        string inFile = m_config["topts"].as<string>();

        if (fs::path(inFile).extension() == ".pts")
        {
            LibUtilities::PtsIOSharedPtr ptsIO =
                MemoryManager<LibUtilities::PtsIO>::AllocateSharedPtr(
                    m_f->m_comm);

            ptsIO->Import(inFile, m_f->m_fieldPts);
        }
        else if (fs::path(inFile).extension() == ".csv")
        {
            LibUtilities::CsvIOSharedPtr csvIO =
                MemoryManager<LibUtilities::CsvIO>::AllocateSharedPtr(
                    m_f->m_comm);

            csvIO->Import(inFile, m_f->m_fieldPts);
        }
        else
        {
            ASSERTL0(false, "unknown topts file type");
        }
    }
    else if (m_config["line"].as<string>().compare("NotSet") != 0)
    {
        vector<NekDouble> values;
        ASSERTL0(
            ParseUtils::GenerateVector(m_config["line"].as<string>(), values),
            "Failed to interpret line string");

        ASSERTL0(values.size() > 2, "line string should contain 2*Dim+1 values "
                                    "N,x0,y0,z0,x1,y1,z1");

        double tmp;
        ASSERTL0(std::modf(values[0], &tmp) == 0.0, "N is not an integer");
        ASSERTL0(values[0] > 1, "N is not a valid number");

        int dim  = (values.size() - 1) / 2;
        int npts = values[0];

        // Information for partitioning
        int ptsPerProc = npts / nprocs;
        int extraPts   = (rank < nprocs - 1) ? 0 : npts % nprocs;
        int locPts     = ptsPerProc + extraPts;
        int start      = rank * ptsPerProc;
        int end        = start + locPts;

        Array<OneD, Array<OneD, NekDouble>> pts(dim);
        Array<OneD, NekDouble> delta(dim);
        for (int i = 0; i < dim; ++i)
        {
            pts[i]   = Array<OneD, NekDouble>(locPts);
            delta[i] = (values[dim + i + 1] - values[i + 1]) / (npts - 1);
        }

        for (int i = 0, cntLoc = 0; i < npts; ++i)
        {
            if (i >= start && i < end)
            {
                for (int n = 0; n < dim; ++n)
                {
                    pts[n][cntLoc] = values[n + 1] + i * delta[n];
                }
                ++cntLoc;
            }
        }

        vector<size_t> ppe;
        ppe.push_back(npts);
        m_f->m_fieldPts =
            MemoryManager<LibUtilities::PtsField>::AllocateSharedPtr(dim, pts);
        m_f->m_fieldPts->SetPtsType(LibUtilities::ePtsLine);
        m_f->m_fieldPts->SetPointsPerEdge(ppe);
    }
    else if (m_config["plane"].as<string>().compare("NotSet") != 0)
    {
        vector<NekDouble> values;
        ASSERTL0(
            ParseUtils::GenerateVector(m_config["plane"].as<string>(), values),
            "Failed to interpret plane string");

        ASSERTL0(values.size() > 9,
                 "plane string should contain 4 Dim+2 values "
                 "N1,N2,x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3");

        double tmp;
        ASSERTL0(std::modf(values[0], &tmp) == 0.0, "N1 is not an integer");
        ASSERTL0(std::modf(values[1], &tmp) == 0.0, "N2 is not an integer");

        ASSERTL0(values[0] > 1, "N1 is not a valid number");
        ASSERTL0(values[1] > 1, "N2 is not a valid number");

        int dim = (values.size() - 2) / 4;

        Array<OneD, int> npts(2);
        npts[0] = values[0];
        npts[1] = values[1];

        int totpts = npts[0] * npts[1];

        // Information for partitioning
        int ptsPerProc = totpts / nprocs;
        int extraPts   = (rank < nprocs - 1) ? 0 : totpts % nprocs;
        int locPts     = ptsPerProc + extraPts;
        int start      = rank * ptsPerProc;
        int end        = start + locPts;

        Array<OneD, Array<OneD, NekDouble>> pts(dim);
        Array<OneD, NekDouble> delta1(dim);
        Array<OneD, NekDouble> delta2(dim);
        for (int i = 0; i < dim; ++i)
        {
            pts[i]    = Array<OneD, NekDouble>(locPts);
            delta1[i] = (values[2 + 1 * dim + i] - values[2 + 0 * dim + i]) /
                        (npts[0] - 1);
            delta2[i] = (values[2 + 2 * dim + i] - values[2 + 3 * dim + i]) /
                        (npts[0] - 1);
        }

        for (int j = 0, cnt = 0, cntLoc = 0; j < npts[1]; ++j)
        {
            for (int i = 0; i < npts[0]; ++i, ++cnt)
            {
                if (cnt >= start && cnt < end)
                {
                    for (int n = 0; n < dim; ++n)
                    {
                        pts[n][cntLoc] =
                            (values[2 + n] + i * delta1[n]) *
                                (1.0 - j / ((NekDouble)(npts[1] - 1))) +
                            (values[2 + 3 * dim + n] + i * delta2[n]) *
                                (j / ((NekDouble)(npts[1] - 1)));
                    }
                    ++cntLoc;
                }
            }
        }

        vector<size_t> ppe;
        ppe.push_back(npts[0]);
        ppe.push_back(npts[1]);
        m_f->m_fieldPts =
            MemoryManager<LibUtilities::PtsField>::AllocateSharedPtr(dim, pts);
        m_f->m_fieldPts->SetPtsType(LibUtilities::ePtsPlane);
        m_f->m_fieldPts->SetPointsPerEdge(ppe);
    }
    else if (m_config["box"].as<string>().compare("NotSet") != 0)
    {
        vector<NekDouble> values;
        ASSERTL0(
            ParseUtils::GenerateVector(m_config["box"].as<string>(), values),
            "Failed to interpret box string");

        ASSERTL0(values.size() == 9, "box string should contain 9 values "
                                     "N1,N2,N3,xmin,xmax,ymin,ymax,zmin,zmax");

        int dim = 3;
        Array<OneD, int> npts(3);
        npts[0] = values[0];
        npts[1] = values[1];
        npts[2] = values[2];

        int totpts = npts[0] * npts[1] * npts[2];

        Array<OneD, Array<OneD, NekDouble>> pts(dim);
        Array<OneD, NekDouble> delta(dim);

        // Information for partitioning
        int ptsPerProc = totpts / nprocs;
        int extraPts   = (rank < nprocs - 1) ? 0 : totpts % nprocs;
        int locPts     = ptsPerProc + extraPts;
        int start      = rank * ptsPerProc;
        int end        = start + locPts;

        for (int i = 0; i < dim; ++i)
        {
            pts[i]   = Array<OneD, NekDouble>(locPts);
            delta[i] = (values[4 + 2 * i] - values[3 + 2 * i]) / (npts[i] - 1);
        }

        for (int k = 0, cnt = 0, cntLoc = 0; k < npts[2]; ++k)
        {
            for (int j = 0; j < npts[1]; ++j)
            {
                for (int i = 0; i < npts[0]; ++i, ++cnt)
                {
                    if (cnt >= start && cnt < end)
                    {
                        pts[0][cntLoc] = values[3] + i * delta[0];
                        pts[1][cntLoc] = values[5] + j * delta[1];
                        pts[2][cntLoc] = values[7] + k * delta[2];
                        ++cntLoc;
                    }
                }
            }
        }

        vector<size_t> ppe;
        ppe.push_back(npts[0]);
        ppe.push_back(npts[1]);
        ppe.push_back(npts[2]);
        m_f->m_fieldPts =
            MemoryManager<LibUtilities::PtsField>::AllocateSharedPtr(dim, pts);
        m_f->m_fieldPts->SetPtsType(LibUtilities::ePtsBox);
        m_f->m_fieldPts->SetPointsPerEdge(ppe);
        vector<NekDouble> boxdim;
        boxdim.assign(&values[3], &values[3] + 6);
        m_f->m_fieldPts->SetBoxSize(boxdim);
    }
    else
    {
        ASSERTL0(false,
                 "ProcessInterpPtsToPts requires line, plane or box option.");
    }
}

void ProcessInterpPtsToPts::InterpolatePtsToPts(
    LibUtilities::PtsFieldSharedPtr &fromPts,
    LibUtilities::PtsFieldSharedPtr &toPts, NekDouble clamp_low,
    NekDouble clamp_up, [[maybe_unused]] NekDouble def_value)
{
    ASSERTL0(toPts->GetNFields() >= fromPts->GetNFields(),
             "ptField has too few fields");

    int nfields = fromPts->GetNFields();

    Interpolator<std::vector<MultiRegions::ExpListSharedPtr>> interp;
    if (m_f->m_comm->GetRank() == 0)
    {
        interp.SetProgressCallback(&ProcessInterpPtsToPts::PrintProgressbar,
                                   this);
    }
    interp.Interpolate(fromPts, toPts);
    if (m_f->m_comm->GetRank() == 0)
    {
        cout << endl;
    }

    for (int f = 0; f < nfields; ++f)
    {
        for (int i = 0; i < toPts->GetNpoints(); ++i)
        {
            if (toPts->GetPointVal(f, i) > clamp_up)
            {
                toPts->SetPointVal(f, i, clamp_up);
            }
            else if (toPts->GetPointVal(f, i) < clamp_low)
            {
                toPts->SetPointVal(f, i, clamp_low);
            }
        }
    }
}

void ProcessInterpPtsToPts::calcCp0()
{
    LibUtilities::PtsFieldSharedPtr pts = m_f->m_fieldPts;
    int dim                             = pts->GetDim();
    int nq1                             = pts->GetNpoints();
    int r, f;
    int pfield = -1;
    NekDouble p0, qinv;
    vector<int> velid;

    vector<NekDouble> values;
    ASSERTL0(ParseUtils::GenerateVector(m_config["cp"].as<string>(), values),
             "Failed to interpret cp string");

    ASSERTL0(values.size() == 2, "cp string should contain 2 values "
                                 "p0 and q (=1/2 rho u^2)");

    p0   = values[0];
    qinv = 1.0 / values[1];

    for (int i = 0; i < pts->GetNFields(); ++i)
    {
        if (boost::iequals(pts->GetFieldName(i), "p"))
        {
            pfield = i;
        }

        if (boost::iequals(pts->GetFieldName(i), "u") ||
            boost::iequals(pts->GetFieldName(i), "v") ||
            boost::iequals(pts->GetFieldName(i), "w"))
        {
            velid.push_back(i);
        }
    }

    if (pfield != -1)
    {
        if (!velid.size())
        {
            WARNINGL0(false, "Did not find velocity components for Cp0");
        }
    }
    else
    {
        WARNINGL0(false, "Failed to find 'p' field to determine cp0");
    }

    // Allocate data storage
    Array<OneD, Array<OneD, NekDouble>> data(2);

    for (f = 0; f < 2; ++f)
    {
        data[f] = Array<OneD, NekDouble>(nq1, 0.0);
    }

    for (r = 0; r < nq1; r++)
    {
        if (pfield != -1) // calculate cp
        {
            data[0][r] = qinv * (pts->GetPointVal(dim + pfield, r) - p0);

            if (velid.size()) // calculate cp0
            {
                NekDouble q = 0;
                for (int i = 0; i < velid.size(); ++i)
                {
                    q += 0.5 * pts->GetPointVal(dim + velid[i], r) *
                         pts->GetPointVal(dim + velid[i], r);
                }
                data[1][r] =
                    qinv * (pts->GetPointVal(dim + pfield, r) + q - p0);
            }
        }
    }

    if (pfield != -1)
    {
        pts->AddField(data[0], "Cp");
        m_f->m_variables.push_back("Cp");
        if (velid.size())
        {
            pts->AddField(data[1], "Cp0");
            m_f->m_variables.push_back("Cp0");
        }
    }
}

void ProcessInterpPtsToPts::PrintProgressbar(const int position,
                                             const int goal) const
{
    LibUtilities::PrintProgressbar(position, goal, "Interpolating");
}
} // namespace Nektar::FieldUtils
