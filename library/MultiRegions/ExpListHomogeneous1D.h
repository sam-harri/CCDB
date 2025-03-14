///////////////////////////////////////////////////////////////////////////////
//
// File: ExpListHomogeneous1D.h
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
// Description: Base class for expansions which are homogeneous in 1
// direction
//
///////////////////////////////////////////////////////////////////////////////

#ifndef EXPLISTHOMO1D_H
#define EXPLISTHOMO1D_H
#include <LibUtilities/Communication/Transposition.h>
#include <LibUtilities/FFT/NektarFFT.h>
#include <MultiRegions/ExpList.h>
#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/MultiRegionsDeclspec.h>
#include <vector>

namespace Nektar::MultiRegions
{

enum Homogeneous1DMatType
{
    eForwardsCoeffSpace1D,
    eBackwardsCoeffSpace1D,
    eForwardsPhysSpace1D,
    eBackwardsPhysSpace1D
};

/// A map between homo matrix keys and their associated block
/// matrices.
typedef std::map<Homogeneous1DMatType, DNekBlkMatSharedPtr>
    Homo1DBlockMatrixMap;
/// A shared pointer to a BlockMatrixMap.
typedef std::shared_ptr<Homo1DBlockMatrixMap> Homo1DBlockMatrixMapShPtr;

// Forward declaration for typedefs
class ExpListHomogeneous1D;

/// Shared pointer to an ExpList3DHomogeneous1D object.
typedef std::shared_ptr<ExpListHomogeneous1D> ExpListHomogeneous1DSharedPtr;
/// Vector of pointers to ExpList3DHomogeneous1D objects.
typedef std::vector<ExpListHomogeneous1DSharedPtr> ExpListHomogeneous1DVector;

/// Abstraction of a two-dimensional multi-elemental expansion which
/// is merely a collection of local expansions.
class ExpListHomogeneous1D : public ExpList
{
public:
    /// Default constructor.
    MULTI_REGIONS_EXPORT ExpListHomogeneous1D(const ExpansionType type);

    MULTI_REGIONS_EXPORT ExpListHomogeneous1D(
        const ExpansionType type,
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const LibUtilities::BasisKey &HomoBasis, const NekDouble lz,
        const bool useFFT, const bool dealiasing);

    /// Copy constructor.
    MULTI_REGIONS_EXPORT ExpListHomogeneous1D(const ExpListHomogeneous1D &In);

    MULTI_REGIONS_EXPORT ExpListHomogeneous1D(
        const ExpListHomogeneous1D &In, const std::vector<unsigned int> &eIDs,
        const Collections::ImplementationType ImpType =
            Collections::eNoImpType);

    /// Destructor.
    MULTI_REGIONS_EXPORT ~ExpListHomogeneous1D() override;

    MULTI_REGIONS_EXPORT void Homogeneous1DTrans(
        const int npts, const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray, bool IsForwards, bool Shuff = true,
        bool UnShuff = true);

    LibUtilities::BasisSharedPtr GetHomogeneousBasis(void)
    {
        return m_homogeneousBasis;
    }

    MULTI_REGIONS_EXPORT void PhysDeriv(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &out_d0, Array<OneD, NekDouble> &out_d1,
        Array<OneD, NekDouble> &out_d2);

    MULTI_REGIONS_EXPORT void PhysDeriv(
        Direction edir, const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &out_d);

    ExpListSharedPtr &GetPlane(int n)
    {
        return m_planes[n];
    }

    LibUtilities::TranspositionSharedPtr m_transposition;
    LibUtilities::CommSharedPtr m_StripZcomm;

protected:
    /// FFT variables
    bool m_useFFT;
    LibUtilities::NektarFFTSharedPtr m_FFT;

    LibUtilities::NektarFFTSharedPtr m_FFT_deal;

    Array<OneD, NekDouble> m_tmpIN;
    Array<OneD, NekDouble> m_tmpOUT;

    /// Definition of the total number of degrees of freedom and
    /// quadrature points. Sets up the storage for \a m_coeff and \a
    ///  m_phys.
    LibUtilities::BasisSharedPtr m_homogeneousBasis;
    NekDouble m_lhom; ///< Width of homogeneous direction
    Homo1DBlockMatrixMapShPtr m_homogeneous1DBlockMat;
    Array<OneD, ExpListSharedPtr> m_planes;

    DNekBlkMatSharedPtr GenHomogeneous1DBlockMatrix(
        Homogeneous1DMatType mattype) const;

    DNekBlkMatSharedPtr GetHomogeneous1DBlockMatrix(
        Homogeneous1DMatType mattype) const;

    NekDouble GetSpecVanVisc(const int k)
    {
        NekDouble returnval = 0.0;

        if (m_specVanVisc.size())
        {
            returnval = m_specVanVisc[k];
        }

        return returnval;
    }

    //  virtual functions
    void v_SetHomo1DSpecVanVisc(Array<OneD, NekDouble> visc) override
    {
        m_specVanVisc = visc;
    }

    size_t v_GetNumElmts(void) override
    {
        return m_planes[0]->GetExpSize();
    }

    LibUtilities::BasisSharedPtr v_GetHomogeneousBasis(void) override
    {
        return GetHomogeneousBasis();
    }

    void v_FwdTrans(const Array<OneD, const NekDouble> &inarray,
                    Array<OneD, NekDouble> &outarray) override;

    void v_FwdTransLocalElmt(const Array<OneD, const NekDouble> &inarray,
                             Array<OneD, NekDouble> &outarray) override;

    void v_FwdTransBndConstrained(const Array<OneD, const NekDouble> &inarray,
                                  Array<OneD, NekDouble> &outarray) override;

    void v_BwdTrans(const Array<OneD, const NekDouble> &inarray,
                    Array<OneD, NekDouble> &outarray) override;

    void v_IProductWRTBase(const Array<OneD, const NekDouble> &inarray,
                           Array<OneD, NekDouble> &outarray) override;

    void v_IProductWRTDerivBase(const int dir,
                                const Array<OneD, const NekDouble> &inarray,
                                Array<OneD, NekDouble> &outarray) override;

    void v_IProductWRTDerivBase(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, NekDouble> &outarray) override;

    std::vector<LibUtilities::FieldDefinitionsSharedPtr> v_GetFieldDefinitions(
        void) override;

    void v_GetFieldDefinitions(
        std::vector<LibUtilities::FieldDefinitionsSharedPtr> &fielddef)
        override;

    void v_AppendFieldData(LibUtilities::FieldDefinitionsSharedPtr &fielddef,
                           std::vector<NekDouble> &fielddata) override;

    void v_AppendFieldData(LibUtilities::FieldDefinitionsSharedPtr &fielddef,
                           std::vector<NekDouble> &fielddata,
                           Array<OneD, NekDouble> &coeffs) override;

    void v_ExtractDataToCoeffs(
        LibUtilities::FieldDefinitionsSharedPtr &fielddef,
        std::vector<NekDouble> &fielddata, std::string &field,
        Array<OneD, NekDouble> &coeffs,
        std::unordered_map<int, int> zIdToPlane) override;

    void v_ExtractCoeffsToCoeffs(const std::shared_ptr<ExpList> &fromExpList,
                                 const Array<OneD, const NekDouble> &fromCoeffs,
                                 Array<OneD, NekDouble> &toCoeffs) override;

    void v_WriteVtkPieceData(std::ostream &outfile, int expansion,
                             std::string var) override;

    void v_PhysInterp1DScaled(const NekDouble scale,
                              const Array<OneD, NekDouble> &inarray,
                              Array<OneD, NekDouble> &outarray) override;

    void v_PhysGalerkinProjection1DScaled(
        const NekDouble scale, const Array<OneD, NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray) override;

    void v_HomogeneousFwdTrans(const int npts,
                               const Array<OneD, const NekDouble> &inarray,
                               Array<OneD, NekDouble> &outarray,
                               bool Shuff = true, bool UnShuff = true) override;

    void v_HomogeneousBwdTrans(const int npts,
                               const Array<OneD, const NekDouble> &inarray,
                               Array<OneD, NekDouble> &outarray,
                               bool Shuff = true, bool UnShuff = true) override;

    void v_DealiasedProd(const int num_dofs,
                         const Array<OneD, NekDouble> &inarray1,
                         const Array<OneD, NekDouble> &inarray2,
                         Array<OneD, NekDouble> &outarray) override;

    void v_DealiasedDotProd(
        const int num_dofs, const Array<OneD, Array<OneD, NekDouble>> &inarray1,
        const Array<OneD, Array<OneD, NekDouble>> &inarray2,
        Array<OneD, Array<OneD, NekDouble>> &outarray) override;

    void v_PhysDeriv(const Array<OneD, const NekDouble> &inarray,
                     Array<OneD, NekDouble> &out_d0,
                     Array<OneD, NekDouble> &out_d1,
                     Array<OneD, NekDouble> &out_d2) override;

    void v_PhysDeriv(Direction edir,
                     const Array<OneD, const NekDouble> &inarray,
                     Array<OneD, NekDouble> &out_d) override;

    LibUtilities::TranspositionSharedPtr v_GetTransposition(void) override;

    Array<OneD, const unsigned int> v_GetZIDs(void) override;

    ExpListSharedPtr &v_GetPlane(int n) override
    {
        return GetPlane(n);
    }

    NekDouble v_GetHomoLen(void) override;

    void v_SetHomoLen(const NekDouble lhom) override;

    NekDouble v_Integral(const Array<OneD, const NekDouble> &inarray) override;

private:
    // Padding operations variables
    bool m_dealiasing;
    int m_padsize;

    /// Spectral vanishing Viscosity coefficient for stabilisation
    Array<OneD, NekDouble> m_specVanVisc;
};

} // namespace Nektar::MultiRegions

#endif // EXPLISTHOMO1D_H
