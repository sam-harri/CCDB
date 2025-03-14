///////////////////////////////////////////////////////////////////////////////
//
// File: Points.hpp
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
// Description: Header file of Points definition
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_FOUNDATIONS_POINTS_H
#define NEKTAR_LIB_UTILITIES_FOUNDATIONS_POINTS_H

#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Foundations/Foundations.hpp>
#include <LibUtilities/Foundations/FoundationsFwd.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrixFwd.hpp>

namespace Nektar::LibUtilities
{

/// Defines a specification for a set of points.
class PointsKey
{
public:
    // Used for looking up the creator. The creator for number of points
    // can generate for any number, so we want the same creator called
    // for all number.
    struct opLess
    {
        LIB_UTILITIES_EXPORT bool operator()(const PointsKey &lhs,
                                             const PointsKey &rhs) const;
    };

    /// Default constructor.
    PointsKey(void)
        : m_numpoints(0), m_pointstype(eNoPointsType),
          m_factor(NekConstants::kNekUnsetDouble)
    {
    }

    /// Constructor defining the number and distribution of points.
    PointsKey(const size_t &numpoints, const PointsType &pointstype,
              const NekDouble factor = NekConstants::kNekUnsetDouble)
        : m_numpoints(numpoints), m_pointstype(pointstype), m_factor(factor)
    {
    }

    /// Destructor.
    virtual ~PointsKey()
    {
    }

    /// Copy constructor.
    PointsKey(const PointsKey &key) = default;

    PointsKey &operator=(const PointsKey &key) = default;

    inline size_t GetNumPoints() const
    {
        return m_numpoints;
    }

    inline PointsType GetPointsType() const
    {
        return m_pointstype;
    }

    inline NekDouble GetFactor() const
    {
        return m_factor;
    }

    inline bool operator==(const PointsKey &key)
    {

        if (fabs(m_factor - key.m_factor) < NekConstants::kNekZeroTol)
        {
            return (m_numpoints == key.m_numpoints &&
                    m_pointstype == key.m_pointstype);
        }

        return false;
    }

    inline bool operator==(const PointsKey *y)
    {
        return (*this == *y);
    }

    inline bool operator!=(const PointsKey &y)
    {
        return (!(*this == y));
    }

    inline bool operator!=(const PointsKey *y)
    {
        return (!(*this == *y));
    }

    // If new points are added, this function must be modified
    inline size_t GetPointsDim() const
    {
        size_t dimpoints = 1;

        switch (m_pointstype)
        {
            case eNodalTriElec:
            case eNodalTriFekete:
            case eNodalTriEvenlySpaced:
            case eNodalTriSPI:
            case eNodalQuadElec:
                dimpoints = 2;
                break;

            case eNodalTetElec:
            case eNodalTetEvenlySpaced:
            case eNodalPrismEvenlySpaced:
            case eNodalPrismElec:
            case eNodalHexElec:
                dimpoints = 3;
                break;

            default:
                break;
        }

        return dimpoints;
    }

    // If new points are added, this function must be modified
    inline size_t GetTotNumPoints() const
    {
        size_t totpoints = m_numpoints;

        switch (m_pointstype)
        {
            case eNodalTriElec:
            case eNodalTriFekete:
            case eNodalTriEvenlySpaced:
                totpoints = m_numpoints * (m_numpoints + 1) / 2;
                break;
            case eNodalTriSPI:
                NEKERROR(ErrorUtil::efatal,
                         "This method cannot be implemented");
                break;

            case eNodalQuadElec:
                totpoints = m_numpoints * m_numpoints;
                break;

            case eNodalTetElec:
            case eNodalTetEvenlySpaced:
                totpoints =
                    m_numpoints * (m_numpoints + 1) * (m_numpoints + 2) / 6;
                break;
            case eNodalTetSPI:
                NEKERROR(ErrorUtil::efatal,
                         "This method cannot be implemented");
                break;

            case eNodalPrismEvenlySpaced:
            case eNodalPrismElec:
                totpoints = m_numpoints * m_numpoints * (m_numpoints + 1) / 2;
                break;
            case eNodalPrismSPI:
                NEKERROR(ErrorUtil::efatal,
                         "This method cannot be implemented");
                break;

            case eNodalHexElec:
                totpoints = m_numpoints * m_numpoints * m_numpoints;
                break;

            default:
                break;
        }

        return totpoints;
    }

    LIB_UTILITIES_EXPORT friend bool operator==(const PointsKey &lhs,
                                                const PointsKey &rhs);
    LIB_UTILITIES_EXPORT friend bool operator<(const PointsKey &lhs,
                                               const PointsKey &rhs);
    LIB_UTILITIES_EXPORT friend bool opLess::operator()(
        const PointsKey &lhs, const PointsKey &rhs) const;

protected:
    size_t m_numpoints;      //!< number of the points (as appropriately
                             //!< defined for PointsType)
    PointsType m_pointstype; //!< Type of Points
    NekDouble m_factor;      //!< optional factor
private:
};

static const PointsKey NullPointsKey(0, eNoPointsType);

LIB_UTILITIES_EXPORT bool operator==(const PointsKey &lhs,
                                     const PointsKey &rhs);
LIB_UTILITIES_EXPORT bool operator<(const PointsKey &lhs, const PointsKey &rhs);
LIB_UTILITIES_EXPORT std::ostream &operator<<(std::ostream &os,
                                              const PointsKey &rhs);

typedef std::vector<PointsKey> PointsKeyVector;

/// Stores a set of points of datatype DataT, defined by a PointKey.
template <typename DataT> class Points
{
public:
    typedef DataT DataType;
    typedef std::shared_ptr<NekMatrix<DataType>> MatrixSharedPtrType;

    virtual ~Points()
    {
    }

    inline void Initialize(void)
    {
        v_Initialize();
    }

    inline size_t GetPointsDim() const
    {
        return m_pointsKey.GetPointsDim();
    }

    inline size_t GetNumPoints() const
    {
        return m_pointsKey.GetNumPoints();
    }

    inline size_t GetTotNumPoints() const
    {
        return m_pointsKey.GetTotNumPoints();
    }

    inline PointsType GetPointsType() const
    {
        return m_pointsKey.GetPointsType();
    }

    inline const Array<OneD, const DataType> &GetZ() const
    {
        return m_points[0];
    }

    inline const Array<OneD, const DataType> &GetW() const
    {
        return m_weights;
    }

    inline void GetZW(Array<OneD, const DataType> &z,
                      Array<OneD, const DataType> &w) const
    {
        z = m_points[0];
        w = m_weights;
    }

    inline const Array<OneD, const NekDouble> &GetBaryWeights() const
    {
        return m_bcweights;
    }

    inline void GetPoints(Array<OneD, const DataType> &x) const
    {
        x = m_points[0];
    }

    inline void GetPoints(Array<OneD, const DataType> &x,
                          Array<OneD, const DataType> &y) const
    {
        x = m_points[0];
        y = m_points[1];
    }

    inline void GetPoints(Array<OneD, const DataType> &x,
                          Array<OneD, const DataType> &y,
                          Array<OneD, const DataType> &z) const
    {
        x = m_points[0];
        y = m_points[1];
        z = m_points[2];
    }

    inline const MatrixSharedPtrType &GetD(Direction dir = xDir) const
    {
        return m_derivmatrix[(int)dir];
    }

    const MatrixSharedPtrType GetI(const PointsKey &key)
    {
        return v_GetI(key);
    }

    const MatrixSharedPtrType GetI(const Array<OneD, const DataType> &x)
    {
        return v_GetI(x);
    }

    const MatrixSharedPtrType GetI(size_t uint,
                                   const Array<OneD, const DataType> &x)
    {
        return v_GetI(uint, x);
    }

    const MatrixSharedPtrType GetI(const Array<OneD, const DataType> &x,
                                   const Array<OneD, const DataType> &y)
    {
        return v_GetI(x, y);
    }

    const MatrixSharedPtrType GetI(const Array<OneD, const DataType> &x,
                                   const Array<OneD, const DataType> &y,
                                   const Array<OneD, const DataType> &z)
    {
        return v_GetI(x, y, z);
    }

    const MatrixSharedPtrType GetGalerkinProjection(const PointsKey &pkey)
    {
        return v_GetGalerkinProjection(pkey);
    }

protected:
    /// Points type for this points distributions.
    PointsKey m_pointsKey;
    /// Storage for the point locations, allowing for up to a 3D points
    /// storage.
    Array<OneD, DataType> m_points[3];
    /// Quadrature weights for the weights.
    Array<OneD, DataType> m_weights;
    /// Barycentric weights.
    Array<OneD, DataType> m_bcweights;
    /// Derivative matrices.
    MatrixSharedPtrType m_derivmatrix[3];
    NekManager<PointsKey, NekMatrix<DataType>, PointsKey::opLess>
        m_InterpManager;
    NekManager<PointsKey, NekMatrix<DataType>, PointsKey::opLess>
        m_GalerkinProjectionManager;

    virtual void v_Initialize(void)
    {
        v_CalculatePoints();
        v_CalculateWeights();
        v_CalculateBaryWeights();
        v_CalculateDerivMatrix();
    }

    virtual void v_CalculatePoints()
    {
        size_t pointsDim    = GetPointsDim();
        size_t totNumPoints = GetTotNumPoints();

        for (size_t i = 0; i < pointsDim; ++i)
        {
            m_points[i] = Array<OneD, DataType>(totNumPoints);
        }
    }

    virtual void v_CalculateWeights()
    {
        m_weights = Array<OneD, DataType>(GetTotNumPoints());
    }

    /**
     * @brief This function calculates the barycentric weights used for
     * enhanced interpolation speed.
     *
     * For the points distribution \f$ z_i \f$ with \f% 1\leq z_i \leq N
     * \f$, the barycentric weights are computed as:
     *
     * \f[
     * b_i=\prod_{\substack{1\leq j\leq N\\ i\neq j}} \frac{1}{z_i-z_j}
     * \f]
     */
    virtual void v_CalculateBaryWeights()
    {
        const size_t totNumPoints = m_pointsKey.GetNumPoints();
        m_bcweights               = Array<OneD, DataType>(totNumPoints, 1.0);

        Array<OneD, DataType> z = m_points[0];

        for (size_t i = 0; i < totNumPoints; ++i)
        {
            for (size_t j = 0; j < totNumPoints; ++j)
            {
                if (i == j)
                {
                    continue;
                }

                m_bcweights[i] *= (z[i] - z[j]);
            }

            m_bcweights[i] = 1.0 / m_bcweights[i];
        }
    }

    virtual void v_CalculateDerivMatrix()
    {
        size_t totNumPoints = GetTotNumPoints();
        for (size_t i = 0; i < m_pointsKey.GetPointsDim(); ++i)
        {
            m_derivmatrix[i] =
                MemoryManager<NekMatrix<DataType>>::AllocateSharedPtr(
                    totNumPoints, totNumPoints);
        }
    }

    Points(const PointsKey &key) : m_pointsKey(key)
    {
    }

    virtual const MatrixSharedPtrType v_GetI(
        [[maybe_unused]] const PointsKey &key)
    {
        NEKERROR(ErrorUtil::efatal, "Method not implemented ");
        std::shared_ptr<NekMatrix<NekDouble>> returnval(
            MemoryManager<NekMatrix<NekDouble>>::AllocateSharedPtr());
        return returnval;
    }

    virtual const MatrixSharedPtrType v_GetI(
        [[maybe_unused]] const Array<OneD, const DataType> &x)
    {
        NEKERROR(ErrorUtil::efatal, "Method not implemented");
        std::shared_ptr<NekMatrix<NekDouble>> returnval(
            MemoryManager<NekMatrix<NekDouble>>::AllocateSharedPtr());
        return returnval;
    }

    virtual const MatrixSharedPtrType v_GetI(
        size_t, [[maybe_unused]] const Array<OneD, const DataType> &x)
    {
        NEKERROR(ErrorUtil::efatal, "Method not implemented");
        std::shared_ptr<NekMatrix<NekDouble>> returnval(
            MemoryManager<NekMatrix<NekDouble>>::AllocateSharedPtr());
        return returnval;
    }

    virtual const MatrixSharedPtrType v_GetI(
        [[maybe_unused]] const Array<OneD, const DataType> &x,
        [[maybe_unused]] const Array<OneD, const DataType> &y)
    {
        NEKERROR(ErrorUtil::efatal, "Method not implemented");
        std::shared_ptr<NekMatrix<NekDouble>> returnval(
            MemoryManager<NekMatrix<NekDouble>>::AllocateSharedPtr());
        return returnval;
    }

    virtual const MatrixSharedPtrType v_GetI(
        [[maybe_unused]] const Array<OneD, const DataType> &x,
        [[maybe_unused]] const Array<OneD, const DataType> &y,
        [[maybe_unused]] const Array<OneD, const DataType> &z)
    {
        NEKERROR(ErrorUtil::efatal, "Method not implemented");
        std::shared_ptr<NekMatrix<NekDouble>> returnval(
            MemoryManager<NekMatrix<NekDouble>>::AllocateSharedPtr());
        return returnval;
    }

    virtual const MatrixSharedPtrType v_GetGalerkinProjection(
        [[maybe_unused]] const PointsKey &pkey)
    {
        NEKERROR(ErrorUtil::efatal, "Method not implemented ");
        std::shared_ptr<NekMatrix<NekDouble>> returnval(
            MemoryManager<NekMatrix<NekDouble>>::AllocateSharedPtr());
        return returnval;
    }

private:
    Points(const Points &pts) = delete;
    Points()                  = delete;
};

} // namespace Nektar::LibUtilities

#endif // NEKTAR_LIB_UTILITIES_FOUNDATIONS_POINTS_H
