////////////////////////////////////////////////////////////////////////////////
//
//  File: GeomElements.cpp
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
//  Description: Python wrapper for GeomElements.
//
////////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Python/NekPyConfig.hpp>
#include <SpatialDomains/HexGeom.h>
#include <SpatialDomains/PointGeom.h>
#include <SpatialDomains/PrismGeom.h>
#include <SpatialDomains/PyrGeom.h>
#include <SpatialDomains/QuadGeom.h>
#include <SpatialDomains/SegGeom.h>
#include <SpatialDomains/TetGeom.h>
#include <SpatialDomains/TriGeom.h>

using namespace Nektar;
using namespace Nektar::SpatialDomains;

template <class T, class S>
std::shared_ptr<T> Geometry_Init(int id, py::list &facets)
{
    std::vector<std::shared_ptr<S>> geomVec;

    for (int i = 0; i < py::len(facets); ++i)
    {
        geomVec.push_back(py::extract<std::shared_ptr<S>>(facets[i]));
    }

    return std::make_shared<T>(id, &geomVec[0]);
}

template <class T, class S>
std::shared_ptr<T> Geometry_Init_Curved(int id, py::list &facets,
                                        CurveSharedPtr curve)
{
    std::vector<std::shared_ptr<S>> geomVec;

    for (int i = 0; i < py::len(facets); ++i)
    {
        geomVec.push_back(py::extract<std::shared_ptr<S>>(facets[i]));
    }

    return std::make_shared<T>(id, &geomVec[0], curve);
}

template <class T, class S> void export_Geom_2d(const char *name)
{
    py::class_<T, py::bases<Geometry2D>, std::shared_ptr<T>>(name, py::init<>())
        .def("__init__", py::make_constructor(
                             &Geometry_Init<T, S>, py::default_call_policies(),
                             (py::arg("id"), py::arg("segments") = py::list())))
        .def("__init__",
             py::make_constructor(
                 &Geometry_Init_Curved<T, S>, py::default_call_policies(),
                 (py::arg("id"), py::arg("segments"), py::arg("curve"))));
}

template <class T, class S> void export_Geom_3d(const char *name)
{
    py::class_<T, py::bases<Geometry3D>, std::shared_ptr<T>>(name, py::init<>())
        .def("__init__",
             py::make_constructor(
                 &Geometry_Init<T, S>, py::default_call_policies(),
                 (py::arg("id"), py::arg("segments") = py::list())));
}

SegGeomSharedPtr SegGeom_Init(int id, int coordim, py::list &points,
                              py::object &curve)
{
    std::vector<PointGeomSharedPtr> geomVec;

    for (int i = 0; i < py::len(points); ++i)
    {
        geomVec.push_back(py::extract<PointGeomSharedPtr>(points[i]));
    }

    if (curve.is_none())
    {
        return std::make_shared<SegGeom>(id, coordim, &geomVec[0]);
    }
    else
    {
        return std::make_shared<SegGeom>(id, coordim, &geomVec[0],
                                         py::extract<CurveSharedPtr>(curve));
    }
}

py::tuple PointGeom_GetCoordinates(const PointGeom &self)
{
    return py::make_tuple(self.x(), self.y(), self.z());
}

void export_GeomElements()
{
    // Geometry dimensioned base classes
    py::class_<Geometry1D, py::bases<Geometry>, std::shared_ptr<Geometry1D>,
               boost::noncopyable>("Geometry1D", py::no_init);
    py::class_<Geometry2D, py::bases<Geometry>, std::shared_ptr<Geometry2D>,
               boost::noncopyable>("Geometry2D", py::no_init)
        .def("GetCurve", &Geometry2D::GetCurve);
    py::class_<Geometry3D, py::bases<Geometry>, std::shared_ptr<Geometry3D>,
               boost::noncopyable>("Geometry3D", py::no_init);

    // Point geometries
    py::class_<PointGeom, py::bases<Geometry>, std::shared_ptr<PointGeom>>(
        "PointGeom", py::init<>())
        .def(py::init<int, int, NekDouble, NekDouble, NekDouble>())
        .def("GetCoordinates", &PointGeom_GetCoordinates);

    // Segment geometries
    py::class_<SegGeom, py::bases<Geometry>, std::shared_ptr<SegGeom>>(
        "SegGeom", py::init<>())
        .def("__init__",
             py::make_constructor(&SegGeom_Init, py::default_call_policies(),
                                  (py::arg("id"), py::arg("coordim"),
                                   py::arg("points") = py::list(),
                                   py::arg("curve")  = py::object())))
        .def("GetCurve", &SegGeom::GetCurve);

    export_Geom_2d<TriGeom, SegGeom>("TriGeom");
    export_Geom_2d<QuadGeom, SegGeom>("QuadGeom");
    export_Geom_3d<TetGeom, TriGeom>("TetGeom");
    export_Geom_3d<PrismGeom, Geometry2D>("PrismGeom");
    export_Geom_3d<PyrGeom, Geometry2D>("PyrGeom");
    export_Geom_3d<HexGeom, QuadGeom>("HexGeom");
}
