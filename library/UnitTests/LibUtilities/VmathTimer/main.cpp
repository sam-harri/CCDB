///////////////////////////////////////////////////////////////////////////////
//
// File: main.cpp
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

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/VmathArray.hpp>
#include <LibUtilities/SimdLib/tinysimd.hpp>

#include <LibUtilities/BasicUtils/Likwid.hpp>
#include <LibUtilities/BasicUtils/Timer.h>

using namespace Nektar;

int main(int argc, char const *argv[])
{

    LIKWID_MARKER_INIT;
    LIKWID_MARKER_THREADINIT;
    LIKWID_MARKER_REGISTER("GathrArray");
    LIKWID_MARKER_REGISTER("GathrScalar");
    LIKWID_MARKER_REGISTER("GathrSimd");

    size_t nPts;
    [[maybe_unused]] int count   = 0;
    [[maybe_unused]] double time = 0.0;

    if (argc < 2)
    {
        nPts = 100;
    }
    else
    {
        nPts = std::stoi(argv[1]);
    }

    std::cout << "number of points\t" << nPts << '\n';

    // number of experiments
    constexpr size_t experiments = 1 << 18;
    {
        // data should be randomized
        Array<OneD, NekDouble> data{nPts * 10, 0.0};
        Array<OneD, NekDouble> dataTrace{nPts, 0.0};
        size_t zero = 0;
        Array<OneD, size_t> indexTrace{nPts, zero};

        // inizialize index (should be randomized)
        for (size_t i = 0; i < nPts; ++i)
        {
            indexTrace[i] = i;
        }

        LIKWID_MARKER_START("GathrArray");
        for (size_t j = 0; j < experiments; ++j)
        {
            // time
            Vmath::Gathr(nPts, data, indexTrace, dataTrace);
        }
        LIKWID_MARKER_STOP("GathrArray");
        // get likwid events
        constexpr short CPU_CLK_UNHALTED_REF_id = 2;
        int nevents{20};
        std::vector<double> events(nevents);
        //
        LIKWID_MARKER_GET("GathrArray", &nevents, events.data(), &time, &count);
        // print out CPE
        double cpeGathrArray =
            events[CPU_CLK_UNHALTED_REF_id] / nPts / experiments;
        std::cout << "GathrArray likwid CPE\t" << cpeGathrArray << '\n';
        std::cout << dataTrace[0] << std::endl;
    }

    {
        // data should be randomized
        Array<OneD, NekDouble> data{nPts * 10, 0.0};
        Array<OneD, NekDouble> dataTrace{nPts, 0.0};
        size_t zero = 0;
        Array<OneD, size_t> indexTrace{nPts, zero};

        // inizialize index (should be randomized)
        for (size_t i = 0; i < nPts; ++i)
        {
            indexTrace[i] = i;
        }

        LIKWID_MARKER_START("GathrSimd");
        for (size_t j = 0; j < experiments; ++j)
        {
            // time
            Vmath::SIMD::Gathr(nPts, data.data(), indexTrace.data(),
                               dataTrace.data());
        }
        LIKWID_MARKER_STOP("GathrSimd");
        // get likwid events
        constexpr short CPU_CLK_UNHALTED_REF_id = 2;
        int nevents{20};
        std::vector<double> events(nevents);
        //
        LIKWID_MARKER_GET("GathrSimd", &nevents, events.data(), &time, &count);
        // print out CPE
        double cpeGathrSimd =
            events[CPU_CLK_UNHALTED_REF_id] / nPts / experiments;
        std::cout << "GathrSimd likwid CPE\t" << cpeGathrSimd << '\n';
        std::cout << dataTrace[0] << std::endl;
    }

    {
        // data should be randomized
        Array<OneD, NekDouble> data{nPts * 10, 0.0};
        Array<OneD, NekDouble> dataTrace{nPts, 0.0};
        size_t zero = 0;
        Array<OneD, size_t> indexTrace{nPts, zero};

        // inizialize index (should be randomized)
        for (size_t i = 0; i < nPts; ++i)
        {
            indexTrace[i] = i;
        }

        LIKWID_MARKER_START("GathrScalar");
        for (size_t j = 0; j < experiments; ++j)
        {
            // time
            Vmath::Gathr(nPts, data.data(), indexTrace.data(),
                         dataTrace.data());
        }
        LIKWID_MARKER_STOP("GathrScalar");
        // get likwid events
        constexpr short CPU_CLK_UNHALTED_REF_id = 2;
        int nevents{20};
        std::vector<double> events(nevents);
        //
        LIKWID_MARKER_GET("GathrScalar", &nevents, events.data(), &time,
                          &count);
        // print out CPE
        double cpeGathrScalar =
            events[CPU_CLK_UNHALTED_REF_id] / nPts / experiments;
        std::cout << "GathrScalar likwid CPE\t" << cpeGathrScalar << '\n';
        std::cout << dataTrace[0] << std::endl;
    }

    LIKWID_MARKER_CLOSE;
}
