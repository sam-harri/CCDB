///////////////////////////////////////////////////////////////////////////////
//
// File: ContField.h
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
// Description: Field definition in tow-dimensions
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBS_MULTIREGIONS_CONTFIELD2D_H
#define NEKTAR_LIBS_MULTIREGIONS_CONTFIELD2D_H

#include <MultiRegions/AssemblyMap/AssemblyMapCG.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <MultiRegions/DisContField.h>
#include <MultiRegions/GJPStabilisation.h>
#include <MultiRegions/GlobalLinSys.h>
#include <MultiRegions/GlobalMatrix.h>
#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/MultiRegionsDeclspec.h>
#include <SpatialDomains/Conditions.h>

namespace Nektar::MultiRegions
{

/// This class is the abstraction of a global continuous two-
/// dimensional spectral/hp element expansion which approximates the
/// solution of a set of partial differential equations.
class ContField : public DisContField
{
public:
    /// The default constructor.
    MULTI_REGIONS_EXPORT ContField();

    /// This constructor sets up global continuous field based on an
    /// input mesh and boundary conditions.
    MULTI_REGIONS_EXPORT ContField(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const SpatialDomains::MeshGraphSharedPtr &graph2D,
        const std::string &variable       = "DefaultVar",
        const bool DeclareCoeffPhysArrays = true,
        const bool CheckIfSingularSystem  = false,
        const Collections::ImplementationType ImpType =
            Collections::eNoImpType);

    /// Construct a global continuous field with solution type based on
    /// another field but using a separate input mesh and boundary
    /// conditions.
    MULTI_REGIONS_EXPORT ContField(
        const ContField &In, const SpatialDomains::MeshGraphSharedPtr &graph2D,
        const std::string &variable, const bool DeclareCoeffPhysArrays = true,
        const bool CheckIfSingularSystem = false);

    /// The copy constructor.
    MULTI_REGIONS_EXPORT ContField(const ContField &In,
                                   bool DeclareCoeffPhysArrays = true);

    /// Copy constructor.
    MULTI_REGIONS_EXPORT ContField(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const ExpList &In);

    /// The default destructor.
    MULTI_REGIONS_EXPORT ~ContField() override;

    /// Assembles the global coefficients \f$\boldsymbol{\hat{u}}_g\f$
    /// from the local coefficients \f$\boldsymbol{\hat{u}}_l\f$.
    inline void Assemble();

    /// Assembles the global coefficients \f$\boldsymbol{\hat{u}}_g\f$
    /// from the local coefficients \f$\boldsymbol{\hat{u}}_l\f$.
    inline void Assemble(const Array<OneD, const NekDouble> &inarray,
                         Array<OneD, NekDouble> &outarray) const;

    /// Returns the map from local to global level.
    inline const AssemblyMapCGSharedPtr &GetLocalToGlobalMap() const;

    /// Solves the two-dimensional Laplace equation, subject to the
    /// boundary conditions specified.
    MULTI_REGIONS_EXPORT void LaplaceSolve(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray,
        const Array<OneD, const NekDouble> &dirForcing = NullNekDouble1DArray,
        const Array<OneD, Array<OneD, NekDouble>> &variablecoeffs =
            NullNekDoubleArrayOfArray,
        NekDouble time = 0.0);

    /// Compute the eigenvalues of the linear advection operator.
    MULTI_REGIONS_EXPORT void LinearAdvectionEigs(
        const NekDouble ax, const NekDouble ay, Array<OneD, NekDouble> &Real,
        Array<OneD, NekDouble> &Imag,
        Array<OneD, NekDouble> &Evecs = NullNekDouble1DArray);

    inline int GetGlobalMatrixNnz(const GlobalMatrixKey &gkey);

    /// Solves the linear system specified by the key \a key.
    MULTI_REGIONS_EXPORT void GlobalSolve(
        const GlobalLinSysKey &key, const Array<OneD, const NekDouble> &rhs,
        Array<OneD, NekDouble> &inout,
        const Array<OneD, const NekDouble> &dirForcing = NullNekDouble1DArray);

    MULTI_REGIONS_EXPORT const GJPStabilisationSharedPtr GetGJPForcing()
    {
        // initialize if required
        if (!m_GJPData)
        {
            m_GJPData = MemoryManager<GJPStabilisation>::AllocateSharedPtr(
                GetSharedThisPtr());
        }

        return m_GJPData;
    }

    MULTI_REGIONS_EXPORT void SetGJPForcing(
        const GJPStabilisationSharedPtr &GJPData)
    {
        m_GJPData = GJPData;
    }

protected:
    // private:
    /// (A shared pointer to) the object which contains all the
    /// required information for the transformation from local to
    /// global degrees of freedom.
    AssemblyMapCGSharedPtr m_locToGloMap;

    /// (A shared pointer to) a list which collects all the global
    /// matrices being assembled, such that they should be constructed
    /// only once.
    GlobalMatrixMapShPtr m_globalMat;

    /// A manager which collects all the global
    /// linear systems being assembled, such that they should be
    /// constructed only once.
    LibUtilities::NekManager<GlobalLinSysKey, GlobalLinSys>
        m_globalLinSysManager;

    /// Data for Gradient Jump Penalisation (GJP) stabilisaiton
    GJPStabilisationSharedPtr m_GJPData;

    /// Returns the global matrix specified by \a mkey.
    MULTI_REGIONS_EXPORT GlobalMatrixSharedPtr
    GetGlobalMatrix(const GlobalMatrixKey &mkey);

    /// Returns the linear system specified by the key \a mkey.
    MULTI_REGIONS_EXPORT GlobalLinSysSharedPtr
    GetGlobalLinSys(const GlobalLinSysKey &mkey);

    MULTI_REGIONS_EXPORT GlobalLinSysSharedPtr
    GenGlobalLinSys(const GlobalLinSysKey &mkey);

    /// Impose the Dirichlet Boundary Conditions on outarray
    MULTI_REGIONS_EXPORT void v_ImposeDirichletConditions(
        Array<OneD, NekDouble> &outarray) override;

    MULTI_REGIONS_EXPORT void v_FillBndCondFromField(
        const Array<OneD, NekDouble> coeffs) override;

    MULTI_REGIONS_EXPORT void v_FillBndCondFromField(
        const int nreg, const Array<OneD, NekDouble> coeffs) override;

    /// Gathers the global coefficients \f$\boldsymbol{\hat{u}}_g\f$
    /// from the local coefficients \f$\boldsymbol{\hat{u}}_l\f$.
    MULTI_REGIONS_EXPORT void v_LocalToGlobal(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray, bool useComm) override;

    MULTI_REGIONS_EXPORT void v_LocalToGlobal(bool useComm) override;

    /// Scatters from the global coefficients
    /// \f$\boldsymbol{\hat{u}}_g\f$ to the local coefficients
    /// \f$\boldsymbol{\hat{u}}_l\f$.
    MULTI_REGIONS_EXPORT void v_GlobalToLocal(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray) override;

    MULTI_REGIONS_EXPORT void v_GlobalToLocal(void) override;

    /// Template method virtual forwarder for FwdTrans().
    MULTI_REGIONS_EXPORT void v_FwdTrans(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray) override;

    /// Template method virtual forwarded for SmoothField().
    MULTI_REGIONS_EXPORT void v_SmoothField(
        Array<OneD, NekDouble> &field) override;

    /// Template method virtual forwarder for MultiplyByInvMassMatrix().
    MULTI_REGIONS_EXPORT void v_MultiplyByInvMassMatrix(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray) override;

    /// Solves the two-dimensional Helmholtz equation, subject to the
    /// boundary conditions specified.
    MULTI_REGIONS_EXPORT GlobalLinSysKey
    v_HelmSolve(const Array<OneD, const NekDouble> &inarray,
                Array<OneD, NekDouble> &outarray,
                const StdRegions::ConstFactorMap &factors,
                const StdRegions::VarCoeffMap &varcoeff,
                const MultiRegions::VarFactorsMap &varfactors,
                const Array<OneD, const NekDouble> &dirForcing,
                const bool PhysSpaceForcing) override;

    // Solve the linear advection problem assuming that m_coeffs
    // vector contains an intial estimate for solution
    MULTI_REGIONS_EXPORT GlobalLinSysKey
    v_LinearAdvectionDiffusionReactionSolve(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray,
        const StdRegions::ConstFactorMap &factors,
        const StdRegions::VarCoeffMap &varcoeff,
        const MultiRegions::VarFactorsMap &varfactors,
        const Array<OneD, const NekDouble> &dirForcing,
        const bool PhysSpaceForcing) override;

    // Solve the linear advection problem assuming that m_coeff
    // vector contains an intial estimate for solution
    MULTI_REGIONS_EXPORT void v_LinearAdvectionReactionSolve(
        const Array<OneD, Array<OneD, NekDouble>> &velocity,
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray, const NekDouble lambda,
        const Array<OneD, const NekDouble> &dirForcing =
            NullNekDouble1DArray) override;

    /// Returns the boundary conditions expansion.
    inline const Array<OneD, const MultiRegions::ExpListSharedPtr> &
    v_GetBndCondExpansions() override;

    /// Template method virtual forwarder for GetBndConditions().
    MULTI_REGIONS_EXPORT const Array<
        OneD, const SpatialDomains ::BoundaryConditionShPtr> &
    v_GetBndConditions() override;
    MULTI_REGIONS_EXPORT void v_ClearGlobalLinSysManager(void) override;

    // Get manager pool count; intended for unit tests
    MULTI_REGIONS_EXPORT int v_GetPoolCount(std::string) override;

    // Remove GlobalLinSys, StaticCond Blocks and LocalMatrix Blocks
    MULTI_REGIONS_EXPORT void v_UnsetGlobalLinSys(GlobalLinSysKey,
                                                  bool) override;
};

typedef std::shared_ptr<ContField> ContFieldSharedPtr;

/**
 * This operation is evaluated as:
 * \f{tabbing}
 * \hspace{1cm}  \= Do \= $e=$  $1, N_{\mathrm{el}}$ \\
 * \> \> Do \= $i=$  $0,N_m^e-1$ \\
 * \> \> \> $\boldsymbol{\hat{u}}_g[\mbox{map}[e][i]] =
 * \boldsymbol{\hat{u}}_g[\mbox{map}[e][i]]+\mbox{sign}[e][i] \cdot
 * \boldsymbol{\hat{u}}^{e}[i]$\\
 * \> \> continue\\
 * \> continue
 * \f}
 * where \a map\f$[e][i]\f$ is the mapping array and \a
 * sign\f$[e][i]\f$ is an array of similar dimensions ensuring the
 * correct modal connectivity between the different elements (both
 * these arrays are contained in the data member #m_locToGloMap). This
 * operation is equivalent to the gather operation
 * \f$\boldsymbol{\hat{u}}_g=\mathcal{A}^{T}\boldsymbol{\hat{u}}_l\f$,
 * where \f$\mathcal{A}\f$ is the
 * \f$N_{\mathrm{eof}}\times N_{\mathrm{dof}}\f$ permutation matrix.
 *
 * @note    The array #m_coeffs should be filled with the local
 *          coefficients \f$\boldsymbol{\hat{u}}_l\f$ and that the
 *          resulting global coefficients \f$\boldsymbol{\hat{u}}_g\f$
 *          will be stored in #m_coeffs.
 */
inline void ContField::Assemble()
{
    m_locToGloMap->Assemble(m_coeffs, m_coeffs);
}

/**
 * This operation is evaluated as:
 * \f{tabbing}
 * \hspace{1cm}  \= Do \= $e=$  $1, N_{\mathrm{el}}$ \\
 * \> \> Do \= $i=$  $0,N_m^e-1$ \\
 * \> \> \> $\boldsymbol{\hat{u}}_g[\mbox{map}[e][i]] =
 * \boldsymbol{\hat{u}}_g[\mbox{map}[e][i]]+\mbox{sign}[e][i] \cdot
 * \boldsymbol{\hat{u}}^{e}[i]$\\
 *  \> \> continue\\
 * \> continue
 * \f}
 * where \a map\f$[e][i]\f$ is the mapping array and \a
 * sign\f$[e][i]\f$ is an array of similar dimensions ensuring the
 * correct modal connectivity between the different elements (both
 * these arrays are contained in the data member #m_locToGloMap). This
 * operation is equivalent to the gather operation
 * \f$\boldsymbol{\hat{u}}_g=\mathcal{A}^{T}\boldsymbol{\hat{u}}_l\f$,
 * where \f$\mathcal{A}\f$ is the
 * \f$N_{\mathrm{eof}}\times N_{\mathrm{dof}}\f$ permutation matrix.
 *
 * @param   inarray     An array of size \f$N_\mathrm{eof}\f$
 *                      containing the local degrees of freedom
 *                      \f$\boldsymbol{x}_l\f$.
 * @param   outarray    The resulting global degrees of freedom
 *                      \f$\boldsymbol{x}_g\f$ will be stored in this
 *                      array of size \f$N_\mathrm{dof}\f$.
 */
inline void ContField::Assemble(const Array<OneD, const NekDouble> &inarray,
                                Array<OneD, NekDouble> &outarray) const
{
    m_locToGloMap->Assemble(inarray, outarray);
}

inline const AssemblyMapCGSharedPtr &ContField::GetLocalToGlobalMap() const
{
    return m_locToGloMap;
}

inline const Array<OneD, const MultiRegions::ExpListSharedPtr> &ContField::
    v_GetBndCondExpansions()
{
    return m_bndCondExpansions;
}

inline const Array<OneD, const SpatialDomains::BoundaryConditionShPtr> &ContField::
    v_GetBndConditions()
{
    return m_bndConditions;
}

inline int ContField::GetGlobalMatrixNnz(const GlobalMatrixKey &gkey)
{
    ASSERTL1(gkey.LocToGloMapIsDefined(),
             "To use method must have a AssemblyMap "
             "attached to key");

    auto matrixIter = m_globalMat->find(gkey);

    if (matrixIter == m_globalMat->end())
    {
        return 0;
    }
    else
    {
        return matrixIter->second->GetNumNonZeroEntries();
    }

    return 0;
}

} // namespace Nektar::MultiRegions

#endif // MULTIERGIONS_CONTFIELD2D_H
