///////////////////////////////////////////////////////////////////////////////
//
// File: CsvIO.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2017 Kilian Lackhove
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
// Description: Csv IO
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_BASIC_UTILS_CSVIO_H
#define NEKTAR_LIB_UTILITIES_BASIC_UTILS_CSVIO_H

#include <memory>
#include <string>
#include <vector>

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/BasicUtils/PtsField.h>
#include <LibUtilities/BasicUtils/PtsIO.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Communication/Comm.h>

namespace Nektar::LibUtilities
{

class CsvIO : public PtsIO
{
public:
    LIB_UTILITIES_EXPORT CsvIO(LibUtilities::CommSharedPtr pComm,
                               bool sharedFilesystem = false);

    LIB_UTILITIES_EXPORT ~CsvIO() override
    {
    }

    LIB_UTILITIES_EXPORT void Write(const std::string &outFile,
                                    const PtsFieldSharedPtr &ptsField,
                                    const bool backup = false);

protected:
    LIB_UTILITIES_EXPORT void v_ImportPtsFieldData(
        const std::string inFile, PtsFieldSharedPtr &ptsField,
        DomainRangeShPtr &Range) override;

    LIB_UTILITIES_EXPORT std::string v_GetFileEnding() const override;
};

typedef std::shared_ptr<CsvIO> CsvIOSharedPtr;
} // namespace Nektar::LibUtilities
#endif
