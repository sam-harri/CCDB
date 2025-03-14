///////////////////////////////////////////////////////////////////////////////
//
// File: NodalTriFekete.h
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

#ifndef NODALTRIFEKETE_H
#define NODALTRIFEKETE_H

#include <iostream>
#include <memory>

#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/Foundations/NodalUtil.h>

namespace Nektar::LibUtilities
{

class NodalTriFekete : public Points<NekDouble>
{
public:
    ~NodalTriFekete() override
    {
    }

    NodalTriFekete(const PointsKey &key) : PointsBaseType(key)
    {
    }

    LIB_UTILITIES_EXPORT static std::shared_ptr<PointsBaseType> Create(
        const PointsKey &key);

protected:
    const MatrixSharedPtrType v_GetI(const PointsKey &pkey) override
    {
        ASSERTL0(pkey.GetPointsDim() == 2,
                 "Fekete Points can only interp to other 2d "
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
        size_t np        = GetTotNumPoints();

        Array<OneD, NekDouble> interp(GetTotNumPoints() * numpoints);
        CalculateInterpMatrix(x, y, interp);

        NekDouble *d = interp.data();
        return MemoryManager<NekMatrix<NekDouble>>::AllocateSharedPtr(numpoints,
                                                                      np, d);
    }

private:
    static bool initPointsManager[];

    std::shared_ptr<NodalUtilTriangle> m_util;

    NodalTriFekete()                             = delete;
    NodalTriFekete(const NodalTriFekete &points) = delete;

    void NodalPointReorder2d();

    void v_CalculatePoints() final;
    void v_CalculateWeights() final;
    void v_CalculateDerivMatrix() final;

    void CalculateInterpMatrix(const Array<OneD, const NekDouble> &xi,
                               const Array<OneD, const NekDouble> &yi,
                               Array<OneD, NekDouble> &interp);
}; // end of NodalTriFekete
} // namespace Nektar::LibUtilities

#endif // NODALTRIFEKETE_H
