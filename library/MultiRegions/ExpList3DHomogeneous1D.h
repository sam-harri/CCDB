///////////////////////////////////////////////////////////////////////////////
//
// File: ExpList3DHomogeneous1D.h
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
// Description: A 3D field which is homogeneous in 1 direction and so
// uses much of the functionality from a ExpList2D and its daughters
//
///////////////////////////////////////////////////////////////////////////////

#ifndef EXPLIST3DHOMO1D_H
#define EXPLIST3DHOMO1D_H

#include <MultiRegions/ExpListHomogeneous1D.h>
#include <MultiRegions/MultiRegionsDeclspec.h>
#include <vector>

namespace Nektar::MultiRegions
{

// Forward declaration for typedefs
class ExpList3DHomogeneous1D;

/// Shared pointer to an ExpList3DHomogeneous1D object.
typedef std::shared_ptr<ExpList3DHomogeneous1D> ExpList3DHomogeneous1DSharedPtr;
/// Vector of pointers to ExpList3DHomogeneous1D objects.
typedef std::vector<ExpList3DHomogeneous1DSharedPtr>
    ExpList3DHomogeneous1DVector;

/// Abstraction of a two-dimensional multi-elemental expansion which
/// is merely a collection of local expansions.
class ExpList3DHomogeneous1D : public ExpListHomogeneous1D
{
public:
    /// Default constructor.
    MULTI_REGIONS_EXPORT ExpList3DHomogeneous1D();

    MULTI_REGIONS_EXPORT ExpList3DHomogeneous1D(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const LibUtilities::BasisKey &HomoBasis, const NekDouble lhom,
        const bool useFFT, const bool dealiasing);

    /// Sets up a list of local expansions based on an input mesh.
    MULTI_REGIONS_EXPORT ExpList3DHomogeneous1D(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const LibUtilities::BasisKey &HomoBasis, const NekDouble lhom,
        const bool useFFT, const bool dealiasing,
        const SpatialDomains::MeshGraphSharedPtr &graph2D,
        const std::string &var = "DefaultVar",
        const Collections::ImplementationType ImpType =
            Collections::eNoImpType);

    /// Sets up a list of local expansions based on an mesh expansion
    MULTI_REGIONS_EXPORT ExpList3DHomogeneous1D(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const LibUtilities::BasisKey &HomoBasis, const NekDouble lhom,
        const bool useFFT, const bool dealiasing,
        const SpatialDomains::ExpansionInfoMap &expansions,
        const Collections::ImplementationType ImpType =
            Collections::eNoImpType);

    /// Copy constructor.
    MULTI_REGIONS_EXPORT ExpList3DHomogeneous1D(
        const ExpList3DHomogeneous1D &In,
        const bool DeclarePlanesSetCoeffPhys = true);

    /// Constructor copying only elements defined in eIds.
    MULTI_REGIONS_EXPORT ExpList3DHomogeneous1D(
        const ExpList3DHomogeneous1D &In, const std::vector<unsigned int> &eIDs,
        const bool DeclarePlanesSetCoeffPhys = true,
        const Collections::ImplementationType ImpType =
            Collections::eNoImpType);

    /// Destructor.
    MULTI_REGIONS_EXPORT ~ExpList3DHomogeneous1D() override;

protected:
    /// Definition of the total number of degrees of freedom and
    /// quadrature points. Sets up the storage for \a m_coeff and \a
    ///  m_phys.
    void SetCoeffPhys(void);

    //  virtual functions
    void v_GetCoords(Array<OneD, NekDouble> &coord_0,
                     Array<OneD, NekDouble> &coord_1,
                     Array<OneD, NekDouble> &coord_2) override;

    void v_GetCoords(const int eid, Array<OneD, NekDouble> &xc0,
                     Array<OneD, NekDouble> &xc1,
                     Array<OneD, NekDouble> &xc2) final;

    void v_WriteTecplotConnectivity(std::ostream &outfile,
                                    int expansion) override;

    void v_WriteVtkPieceHeader(std::ostream &outfile, int expansion,
                               int istrip) override;

    NekDouble v_L2(const Array<OneD, const NekDouble> &inarray,
                   const Array<OneD, const NekDouble> &soln =
                       NullNekDouble1DArray) override;

    Array<OneD, const NekDouble> v_HomogeneousEnergy(void) override;

    void v_GetPeriodicEntities(
        PeriodicMap &periodicVerts, PeriodicMap &periodicEdges,
        [[maybe_unused]] PeriodicMap &periodicFaces) override
    {
        m_planes[0]->GetPeriodicEntities(periodicVerts, periodicEdges);
    }

private:
    MULTI_REGIONS_EXPORT void GenExpList3DHomogeneous1D(
        const SpatialDomains::ExpansionInfoMap &expansions,
        const Collections::ImplementationType ImpType);
};

} // namespace Nektar::MultiRegions

#endif // EXPLIST3DHOMO1D_H
