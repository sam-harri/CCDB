///////////////////////////////////////////////////////////////////////////////
//
// File: NodalTriElec.h
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
// Description: Header file of 2D Nodal Triangle Fekete Points
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NODALTRIELEC_H
#define NODALTRIELEC_H

#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/Foundations/NodalUtil.h>
#include <memory>

namespace Nektar::LibUtilities
{

class NodalTriElec : public Points<NekDouble>
{
public:
    ~NodalTriElec() override
    {
    }

    LIB_UTILITIES_EXPORT static std::shared_ptr<PointsBaseType> Create(
        const PointsKey &key);

    NodalTriElec(const PointsKey &key) : PointsBaseType(key)
    {
    }

protected:
    const MatrixSharedPtrType v_GetI(const PointsKey &pkey) override
    {
        ASSERTL0(pkey.GetPointsDim() == 2,
                 "NodalTriElec Points can only interpolate to other 2d "
                 "point distributions");
        Array<OneD, const NekDouble> x, y;
        PointsManager()[pkey]->GetPoints(x, y);
        return GetI(x, y);
    }

    const MatrixSharedPtrType v_GetI(
        const Array<OneD, const NekDouble> &x,
        const Array<OneD, const NekDouble> &y) override
    {
        size_t numpoints = x.size();
        Array<OneD, NekDouble> interp(GetTotNumPoints() * numpoints);
        CalculateInterpMatrix(x, y, interp);

        size_t np    = GetTotNumPoints();
        NekDouble *d = interp.data();
        return MemoryManager<NekMatrix<NekDouble>>::AllocateSharedPtr(numpoints,
                                                                      np, d);
    }

private:
    static bool initPointsManager[];

    std::shared_ptr<NodalUtilTriangle> m_util;

    NodalTriElec()                           = delete;
    NodalTriElec(const NodalTriElec &points) = delete;

    void NodalPointReorder2d();

    void v_CalculatePoints() final;
    void v_CalculateWeights() final;
    void v_CalculateDerivMatrix() final;

    void CalculateInterpMatrix(const Array<OneD, const NekDouble> &xia,
                               const Array<OneD, const NekDouble> &yia,
                               Array<OneD, NekDouble> &interp);
};
} // namespace Nektar::LibUtilities

#endif // NODALTRIELEC_H
