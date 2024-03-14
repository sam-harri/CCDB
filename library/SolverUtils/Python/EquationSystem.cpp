///////////////////////////////////////////////////////////////////////////////
//
// File: EquationSystem.cpp
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
// Description: Python wrapper for the EquationSystem.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Python/NekPyConfig.hpp>
#include <SolverUtils/EquationSystem.h>

using namespace Nektar;
using namespace Nektar::SolverUtils;

/**
 * @brief Dummy equation system that can be used for Python testing.
 *
 * This class simply evaluates a known function as the implementation as part of
 * the DoInitialise method.
 */
class DummyEquationSystem : public EquationSystem
{
public:
    friend class MemoryManager<DummyEquationSystem>;

    /// Creates an instance of this class
    static SolverUtils::EquationSystemSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const SpatialDomains::MeshGraphSharedPtr &pGraph)
    {
        EquationSystemSharedPtr p =
            MemoryManager<DummyEquationSystem>::AllocateSharedPtr(pSession,
                                                                  pGraph);
        p->InitObject();
        return p;
    }
    /// Name of class
    static std::string className;

    ~DummyEquationSystem() override
    {
    }

protected:
    DummyEquationSystem(const LibUtilities::SessionReaderSharedPtr &pSession,
                        const SpatialDomains::MeshGraphSharedPtr &pGraph)
        : EquationSystem(pSession, pGraph)
    {
    }

    void v_InitObject(bool DeclareField) override
    {
        EquationSystem::v_InitObject(DeclareField);

        for (int i = 0; i < m_fields.size(); ++i)
        {
            int nq = m_fields[0]->GetNpoints();

            Array<OneD, NekDouble> x(nq), y(nq), z(nq);
            m_fields[i]->GetCoords(x, y, z);

            for (int j = 0; j < nq; ++j)
            {
                m_fields[i]->UpdatePhys()[j] = sin(x[j]) * sin(y[j]);
            }

            m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),
                                  m_fields[i]->UpdateCoeffs());
        }
    }
};

std::string DummyEquationSystem::className =
    SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
        "Dummy", DummyEquationSystem::create, "Dummy equation system.");

EquationSystemSharedPtr EquationSystem_Create(
    std::string eqnSysName, LibUtilities::SessionReaderSharedPtr session,
    SpatialDomains::MeshGraphSharedPtr mesh)
{
    return GetEquationSystemFactory().CreateInstance(eqnSysName, session, mesh);
}

py::list EquationSystem_GetFields(EquationSystemSharedPtr eqSys)
{
    py::list expLists;
    auto &fields = eqSys->UpdateFields();

    for (int i = 0; i < fields.size(); ++i)
    {
        expLists.append(py::object(fields[i]));
    }

    return expLists;
}

void export_EquationSystem()
{
    py::class_<EquationSystem, std::shared_ptr<EquationSystem>,
               boost::noncopyable>("EquationSystem", py::no_init)
        .def("InitObject", &EquationSystem::InitObject)
        .def("DoInitialise", &EquationSystem::DoInitialise)
        .def("DoSolve", &EquationSystem::DoSolve)

        .def("GetFields", &EquationSystem_GetFields)

        // Factory functions.
        .def("Create", &EquationSystem_Create)
        .staticmethod("Create");
}
