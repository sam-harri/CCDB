///////////////////////////////////////////////////////////////////////////////
//
// File: TimeIntegrationScheme.cpp
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
// Description: implementation of time integration scheme class
// this class is the parent class for the TimeIntegrationSchemeGLM,
// TimeIntegrationSchemeFIT, TimeIntegrationSchemeGEM, and
// TimeIntegrationSchemeSDC classes.
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>

namespace Nektar::LibUtilities
{

TimeIntegrationSchemeFactory &GetTimeIntegrationSchemeFactory()
{
    static TimeIntegrationSchemeFactory instance;
    return instance;
}

TimeIntegrationScheme::TimeIntegrationScheme(std::string variant, size_t order,
                                             std::vector<NekDouble> freeParams)
{
    boost::ignore_unused(variant, order, freeParams);
}

// Access Methods
std::string TimeIntegrationScheme::v_GetFullName() const
{
    return GetName() + GetVariant() + "Order" + std::to_string(GetOrder());
}

// Friend Operators
std::ostream &operator<<(std::ostream &os, const TimeIntegrationScheme &rhs)
{
    rhs.print(os);

    return os;
}

std::ostream &operator<<(std::ostream &os,
                         const TimeIntegrationSchemeSharedPtr &rhs)
{
    os << *rhs.get();

    return os;
}

} // namespace Nektar::LibUtilities
