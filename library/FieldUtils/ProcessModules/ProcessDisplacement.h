////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessDisplacement.h
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
//  Description: Deforms a mesh given input field defining displacement.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef FIELDUTILS_PROCESSDISPLACEMENT
#define FIELDUTILS_PROCESSDISPLACEMENT

#include "ProcessBoundaryExtract.h"

namespace Nektar::FieldUtils
{

class ProcessDisplacement : public ProcessBoundaryExtract
{
public:
    /// Creates an instance of this class
    static std::shared_ptr<Module> create(FieldSharedPtr f)
    {
        return MemoryManager<ProcessDisplacement>::AllocateSharedPtr(f);
    }
    static ModuleKey className;

    ProcessDisplacement(FieldSharedPtr f);
    ~ProcessDisplacement() override;

protected:
    /// Write mesh to output file.
    void v_Process(po::variables_map &vm) override;

    std::string v_GetModuleName() override
    {
        return "ProcessDisplacement";
    }

    std::string v_GetModuleDescription() override
    {
        return "Calculating displacement";
    }
};
} // namespace Nektar::FieldUtils

#endif
