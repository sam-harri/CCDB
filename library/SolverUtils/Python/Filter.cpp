///////////////////////////////////////////////////////////////////////////////
//
// File: Filter.cpp
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
// Description: Python wrapper for Filter.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Python/NekPyConfig.hpp>
#include <SolverUtils/Filters/Filter.h>

using namespace Nektar;
using namespace Nektar::SolverUtils;

/**
 * @brief Filter wrapper to handle virtual function calls in @c Filter and its
 * subclasses.
 */
struct FilterWrap : public Filter, public py::wrapper<Filter>
{
    /**
     * @brief Constructor, which is identical to Filter.
     *
     * @param field  Input field.
     */
    FilterWrap(LibUtilities::SessionReaderSharedPtr session,
               std::weak_ptr<EquationSystem> eqn)
        : Filter(session, eqn), py::wrapper<Filter>()
    {
    }

    py::list ArrayOneDToPyList(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields)
    {
        py::list expLists;

        for (int i = 0; i < pFields.size(); ++i)
        {
            expLists.append(py::object(pFields[i]));
        }

        return expLists;
    }

    void v_Initialise(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time) override
    {
        this->get_override("Initialise")(ArrayOneDToPyList(pFields), time);
    }

    void v_Update(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time) override
    {
        this->get_override("Update")(ArrayOneDToPyList(pFields), time);
    }

    void v_Finalise(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time) override
    {
        this->get_override("Finalise")(ArrayOneDToPyList(pFields), time);
    }

    bool v_IsTimeDependent() override
    {
        return this->get_override("IsTimeDependent")();
    }
};

/**
 * @brief Lightweight wrapper for Filter factory creation function.
 */
FilterSharedPtr Filter_Create(py::tuple args, py::dict kwargs)
{
    using NekError = ErrorUtil::NekError;

    if (py::len(args) != 3)
    {
        throw NekError("Filter.Create() requires three arguments: "
                       "filter name, a SessionReader object and an "
                       "EquationSystem object.");
    }

    std::string filterName = py::extract<std::string>(args[0]);

    if (!py::extract<LibUtilities::SessionReaderSharedPtr>(args[1]).check())
    {
        throw NekError("Second argument to Create() should be a SessionReader "
                       "object.");
    }

    if (!py::extract<EquationSystemSharedPtr>(args[2]).check())
    {
        throw NekError("Second argument to Create() should be a EquationSystem "
                       "object.");
    }

    LibUtilities::SessionReaderSharedPtr session =
        py::extract<LibUtilities::SessionReaderSharedPtr>(args[1]);
    EquationSystemSharedPtr eqn = py::extract<EquationSystemSharedPtr>(args[2]);

    // Process keyword arguments.
    py::list items = kwargs.items();
    Filter::ParamMap params;

    for (int i = 0; i < py::len(items); ++i)
    {
        std::string arg = py::extract<std::string>(items[i][0]);
        std::string val =
            py::extract<std::string>(items[i][1].attr("__str__")());
        params[arg] = val;
    }

    std::shared_ptr<Filter> filter =
        GetFilterFactory().CreateInstance(filterName, session, eqn, params);

    return filter;
}

class FilterRegisterHelper
{
public:
    FilterRegisterHelper(py::object obj) : m_obj(obj)
    {
        py::incref(obj.ptr());
    }

    ~FilterRegisterHelper()
    {
        py::decref(m_obj.ptr());
    }

    FilterSharedPtr create(const LibUtilities::SessionReaderSharedPtr &session,
                           const std::weak_ptr<EquationSystem> &eqSys,
                           const Filter::ParamMap &params)
    {
        py::dict args;
        for (auto &param : params)
        {
            args[py::str(param.first)] = py::str(param.second);
        }
        py::object inst = m_obj(session, eqSys.lock(), args);
        return py::extract<FilterSharedPtr>(inst);
    }

protected:
    py::object m_obj;
};

#if PY_MAJOR_VERSION == 2
void FilterCapsuleDestructor(void *ptr)
{
    FilterRegisterHelper *tmp = (FilterRegisterHelper *)ptr;
    delete tmp;
}
#else
void FilterCapsuleDestructor(PyObject *ptr)
{
    FilterRegisterHelper *tmp =
        (FilterRegisterHelper *)PyCapsule_GetPointer(ptr, nullptr);
    delete tmp;
}
#endif

/**
 * @brief Lightweight wrapper for the Filter factory RegisterCreatorFunction, to
 * support the ability for Python subclasses of Filter to register themselves
 * with the Nektar++ Filter factory.
 *
 * This function wraps the NekFactory RegisterCreatorFunction. This function
 * expects a function pointer to a C++ object that will construct a Filter. In
 * this case we therefore need to construct a function call that will construct
 * our Python object (which is a subclass of Filter), and then pass this back to
 * Boost.Python to give the Python object back.
 *
 * We have to do some indirection here to get this to work, but we can
 * achieve this with the following strategy:
 *
 * - Create a @c FilterRegisterHelper object, which as an argument will store
 *   the Python class instance that will be instantiated from the Python side.
 * - Using std::bind, construct a function pointer to the helper's creation
 *   function, FilterRegisterHelper::create.
 * - Create a Python capsule that will contain the @c FilterRegisterHelper
 *   instance, and register this in the global namespace of the current
 *   module. This then ties the capsule to the lifetime of the module.
 */
void Filter_Register(std::string const &filterName, py::object &obj)
{
    // Create a filter register helper, which will call the C++ function to
    // create the filter.
    FilterRegisterHelper *helper = new FilterRegisterHelper(obj);

    // Register this with the filter factory using std::bind to grab a function
    // pointer to that particular object's function.
    GetFilterFactory().RegisterCreatorFunction(
        filterName,
        std::bind(&FilterRegisterHelper::create, helper, std::placeholders::_1,
                  std::placeholders::_2, std::placeholders::_3));

    // Create a capsule that will be embedded in the __main__ namespace. So
    // deallocation will occur, but only once Python ends or the Python module
    // is deallocated.
    std::string filterkey = "_" + filterName;

#if PY_MAJOR_VERSION == 2
    py::object capsule(
        py::handle<>(PyCObject_FromVoidPtr(helper, FilterCapsuleDestructor)));
#else
    py::object capsule(
        py::handle<>(PyCapsule_New(helper, nullptr, FilterCapsuleDestructor)));
#endif

    // Embed this in __main__.
    py::import("__main__").attr(filterkey.c_str()) = capsule;
}

struct FilterWrapConverter
{
    FilterWrapConverter()
    {
        // An important bit of code which will register allow shared_ptr<Filter>
        // as something that boost::python recognises, otherwise modules
        // constructed from the factory will not work from Python.
        py::objects::class_value_wrapper<
            std::shared_ptr<Filter>,
            py::objects::make_ptr_instance<
                Filter, py::objects::pointer_holder<std::shared_ptr<Filter>,
                                                    Filter>>>();
    }
};

inline Array<OneD, MultiRegions::ExpListSharedPtr> PyListToOneDArray(
    py::list &pyExpList)
{
    using NekError = ErrorUtil::NekError;
    Array<OneD, MultiRegions::ExpListSharedPtr> expLists(py::len(pyExpList));

    for (int i = 0; i < expLists.size(); ++i)
    {
        if (!py::extract<MultiRegions::ExpListSharedPtr>(pyExpList[i]).check())
        {
            throw NekError("List should contain only ExpList objects.");
        }

        expLists[i] = py::extract<MultiRegions::ExpListSharedPtr>(pyExpList[i]);
    }

    return expLists;
}

void Filter_Initialise(std::shared_ptr<Filter> filter, py::list expLists,
                       NekDouble time)
{
    filter->Initialise(PyListToOneDArray(expLists), time);
}

void Filter_Update(std::shared_ptr<Filter> filter, py::list expLists,
                   NekDouble time)
{
    filter->Update(PyListToOneDArray(expLists), time);
}

void Filter_Finalise(std::shared_ptr<Filter> filter, py::list expLists,
                     NekDouble time)
{
    filter->Finalise(PyListToOneDArray(expLists), time);
}

bool Filter_IsTimeDependent(std::shared_ptr<Filter> filter)
{
    return filter->IsTimeDependent();
}

void export_Filter()
{
    // Define FilterWrap to be implicitly convertible to a Filter, since it
    // seems that doesn't sometimes get picked up.
    py::implicitly_convertible<std::shared_ptr<FilterWrap>,
                               std::shared_ptr<Filter>>();

    // Wrapper for the Filter class. Note that since Filter contains a pure
    // virtual function, we need the FilterWrap helper class to handle this for
    // us. In the lightweight wrappers above, we therefore need to ensure we're
    // passing std::shared_ptr<Filter> as the first argument, otherwise they
    // won't accept objects constructed from Python.
    py::class_<FilterWrap, std::shared_ptr<FilterWrap>, boost::noncopyable>(
        "Filter", py::init<LibUtilities::SessionReaderSharedPtr,
                           EquationSystemSharedPtr>())

        .def("Initialise", &Filter_Initialise)
        .def("Update", &Filter_Update)
        .def("Finalise", &Filter_Finalise)
        .def("IsTimeDependent", &Filter_IsTimeDependent)

        // Factory functions.
        .def("Create", py::raw_function(Filter_Create))
        .staticmethod("Create")
        .def("Register", &Filter_Register)
        .staticmethod("Register");

    FilterWrapConverter();
}
