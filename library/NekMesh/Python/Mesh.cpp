///////////////////////////////////////////////////////////////////////////////
//
// File: Mesh.cpp
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
// Description: Python wrapper for Mesh.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Python/NekPyConfig.hpp>
#include <NekMesh/MeshElements/Mesh.h>
#include <NekMesh/Python/NekMesh.h>

using namespace Nektar;
using namespace Nektar::NekMesh;

void export_Mesh(py::module &m)
{
    // Wrapper for ElementMap based loosely on stl_bind.h.
    py::class_<ElementMap>(m, "ElementMap")
        .def(py::init<>())
        .def("__len__", [](const ElementMap &v) { return v.size(); })
        .def(
            "__iter__",
            [](ElementMap &v) { return py::make_iterator(v.begin(), v.end()); },
            py::keep_alive<0, 1>())
        .def(
            "__getitem__",
            [](ElementMap &v, int &key) -> std::vector<ElementSharedPtr> & {
                if (key < 0 || key > 3)
                {
                    throw py::index_error();
                }
                return v[key];
            },
            py::return_value_policy::reference_internal);

    py::class_<Mesh, std::shared_ptr<Mesh>>(m, "Mesh")
        .def(py::init<>())
        .def_readwrite("node", &Mesh::m_vertexSet)
        .def_readwrite("expDim", &Mesh::m_expDim)
        .def_readwrite("spaceDim", &Mesh::m_spaceDim)
        .def_readwrite("nummode", &Mesh::m_nummode)
        .def_readonly("element", &Mesh::m_element);
}
