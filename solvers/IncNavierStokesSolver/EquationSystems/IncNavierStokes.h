///////////////////////////////////////////////////////////////////////////////
//
// File: IncNavierStokes.h
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
// Description: Basic Advection Diffusion Reaction Field definition in
// two-dimensions
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_INCNAVIERSTOKES_H
#define NEKTAR_SOLVERS_INCNAVIERSTOKES_H
#include <IncNavierStokesSolver/BoundaryConditions/IncBoundaryConditions.h>
#include <IncNavierStokesSolver/EquationSystems/Extrapolate.h>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/VmathArray.hpp>
#include <SolverUtils/AdvectionSystem.h>
#include <SolverUtils/Filters/FilterInterfaces.hpp>
#include <SolverUtils/Forcing/Forcing.h>
#include <SolverUtils/UnsteadySystem.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <complex>

namespace Nektar
{
namespace bnu = boost::numeric::ublas;

enum SolverType
{
    eNoSolverType,
    eCoupledLinearisedNS,
    eSmoothedProfileMethod,
    eVCSMapping,
    eVelocityCorrectionScheme,
    eVCSWeakPressure,
    eWeakPressure,
    eVCSImplicit,
    eImplicit,
};

enum EquationType
{
    eNoEquationType,
    eSteadyStokes,
    eSteadyOseen,
    eSteadyLinearisedNS,
    eUnsteadyStokes,
    eUnsteadyLinearisedNS,
    eUnsteadyNavierStokes,
    eSteadyNavierStokes,
    eEquationTypeSize
};

// Keep this consistent with the enums in EquationType.
const std::string kEquationTypeStr[] = {
    "NoType",
    "SteadyStokes",
    "SteadyOseen",
    "SteadyLinearisedNS",
    "UnsteadyStokes",
    "UnsteadyLinearisedNS",
    "UnsteadyNavierStokes",
    "SteadyNavierStokes",
};

enum AdvectionForm
{
    eNoAdvectionForm,
    eConvective,
    eNonConservative,
    eLinearised,
    eAdjoint,
    eSkewSymmetric,
    eNoAdvection,
    eAdvectionFormSize
};

// Keep this consistent with the enums in EquationType.
const std::string kAdvectionFormStr[] = {"NoType",
                                         "Convective",
                                         "NonConservative",
                                         "Linearised",
                                         "Adjoint",
                                         "SkewSymmetric"
                                         "NoAdvection"};

typedef std::complex<double> NekComplexDouble;

struct WomersleyParams
{
    WomersleyParams(int dim)
    {
        m_axisnormal = Array<OneD, NekDouble>(dim, 0.0);
        m_axispoint  = Array<OneD, NekDouble>(dim, 0.0);
    };

    virtual ~WomersleyParams(){};

    // Real and imaginary velocity comp. of wom
    std::vector<NekComplexDouble> m_wom_vel;

    // Womersley  BC constants
    NekDouble m_radius;
    NekDouble m_period;
    Array<OneD, NekDouble> m_axisnormal;
    // currently this needs to be the point in the middle of the
    // axis but should be generalised to be any point on the axis
    Array<OneD, NekDouble> m_axispoint;

    // poiseuille flow and fourier coefficients
    Array<OneD, Array<OneD, NekDouble>> m_poiseuille;
    Array<OneD, Array<OneD, Array<OneD, NekComplexDouble>>> m_zvel;
};
typedef std::shared_ptr<WomersleyParams> WomersleyParamsSharedPtr;

/**
 * \brief This class is the base class for Navier Stokes problems
 *
 */
class IncNavierStokes : public SolverUtils::AdvectionSystem,
                        public SolverUtils::FluidInterface
{
public:
    // Destructor
    ~IncNavierStokes() override;

    void v_InitObject(bool DeclareField = true) override;

    int GetNConvectiveFields(void)
    {
        return m_nConvectiveFields;
    }

    void AddForcing(const SolverUtils::ForcingSharedPtr &pForce);

    void v_GetPressure(
        const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        Array<OneD, NekDouble> &pressure) override;

    void v_GetDensity(
        const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        Array<OneD, NekDouble> &density) override;

    bool v_HasConstantDensity() override
    {
        return true;
    }

    void v_GetVelocity(
        const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        Array<OneD, Array<OneD, NekDouble>> &velocity) override;

    void v_SetMovingFrameVelocities(const Array<OneD, NekDouble> &vFrameVels,
                                    const int step) override;
    bool v_GetMovingFrameVelocities(Array<OneD, NekDouble> &vFrameVels,
                                    const int step) override;
    void v_SetMovingFrameDisp(const Array<OneD, NekDouble> &vFrameDisp,
                              const int step) override;
    void v_SetMovingFramePivot(
        const Array<OneD, NekDouble> &vFramePivot) override;
    bool v_GetMovingFrameDisp(Array<OneD, NekDouble> &vFrameDisp,
                              const int step) override;
    void v_SetAeroForce(Array<OneD, NekDouble> forces) override;
    void v_GetAeroForce(Array<OneD, NekDouble> forces) override;

    bool DefinedForcing(const std::string &sForce);

protected:
    // pointer to the extrapolation class for sub-stepping and HOPBS

    ExtrapolateSharedPtr m_extrapolation;
    IncBoundaryConditionsSharedPtr m_IncNavierStokesBCs;

    /// modal energy file
    std::ofstream m_mdlFile;

    /// bool to identify if advection term smoothing is requested
    bool m_SmoothAdvection;

    /// Forcing terms
    std::vector<SolverUtils::ForcingSharedPtr> m_forcing;

    /// Number of fields to be convected;
    int m_nConvectiveFields;

    /// int which identifies which components of m_fields contains the
    /// velocity (u,v,w);
    Array<OneD, int> m_velocity;

    /// Pointer to field holding pressure field
    MultiRegions::ExpListSharedPtr m_pressure;
    /// Kinematic viscosity
    NekDouble m_kinvis;
    /// dump energy to file at steps time
    int m_energysteps;

    /// equation type;
    EquationType m_equationType;

    /// Mapping from BCs to Elmt IDs
    Array<OneD, Array<OneD, int>> m_fieldsBCToElmtID;
    /// Mapping from BCs to Elmt Edge IDs
    Array<OneD, Array<OneD, int>> m_fieldsBCToTraceID;
    /// RHS Factor for Radiation Condition
    Array<OneD, Array<OneD, NekDouble>> m_fieldsRadiationFactor;

    /// Number of time integration steps AND Order of extrapolation for
    /// pressure boundary conditions.
    int m_intSteps;

    /// pivot point for moving reference frame
    Array<OneD, NekDouble> m_pivotPoint;
    Array<OneD, NekDouble> m_aeroForces;

    /// Constructor.
    IncNavierStokes(const LibUtilities::SessionReaderSharedPtr &pSession,
                    const SpatialDomains::MeshGraphSharedPtr &pGraph);

    EquationType GetEquationType(void)
    {
        return m_equationType;
    }
    static std::string eqTypeLookupIds[];

    void EvaluateAdvectionTerms(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time);

    void WriteModalEnergy(void);

    /// time dependent boundary conditions updating
    void SetBoundaryConditions(NekDouble time);

    /// Set Radiation forcing term
    void SetRadiationBoundaryForcing(int fieldid);

    /// Set Normal Velocity Component to Zero
    void SetZeroNormalVelocity();

    /// Set Womersley Profile if specified
    void SetWomersleyBoundary(const int fldid, const int bndid);

    /// Set Up Womersley details
    void SetUpWomersley(const int fldid, const int bndid, std::string womstr);

    /// Womersley parameters if required
    std::map<int, std::map<int, WomersleyParamsSharedPtr>> m_womersleyParams;

    MultiRegions::ExpListSharedPtr v_GetPressure() override
    {
        return m_pressure;
    }

    void v_TransCoeffToPhys(void) override
    {
        ASSERTL0(false, "This method is not defined in this class");
    }

    void v_TransPhysToCoeff(void) override
    {
        ASSERTL0(false, "This method is not defined in this class");
    }

    virtual int v_GetForceDimension() = 0;

    Array<OneD, NekDouble> v_GetMaxStdVelocity(
        const NekDouble SpeedSoundFactor) override;

    bool v_PreIntegrate(int step) override;

private:
};

typedef std::shared_ptr<IncNavierStokes> IncNavierStokesSharedPtr;

} // namespace Nektar

#endif // NEKTAR_SOLVERS_INCNAVIERSTOKES_H
