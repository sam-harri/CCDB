///////////////////////////////////////////////////////////////////////////////
//
// File: CompressibleFlowSystemImplicit.h
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
// Description: Auxiliary functions for the incompressible
// compressible flow system
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_COMPRESSIBLEFLOWSYSTEMIMPLICIT_H
#define NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_COMPRESSIBLEFLOWSYSTEMIMPLICIT_H

#include <CompressibleFlowSolver/EquationSystems/CompressibleFlowSystem.h>
#include <CompressibleFlowSolver/Preconditioner/PreconCfs.h>
#include <CompressibleFlowSolver/Preconditioner/PreconCfsOp.h>
#include <LibUtilities/LinearAlgebra/NekNonlinSysIter.h>

namespace Nektar
{
/**
 *
 */
class CFSImplicit : virtual public CompressibleFlowSystem

{
public:
    CFSImplicit(const LibUtilities::SessionReaderSharedPtr &pSession,
                const SpatialDomains::MeshGraphSharedPtr &pGraph);

    ~CFSImplicit() override = default;

protected:
    bool m_viscousJacFlag;
    bool m_advectionJacFlag;
    bool m_flagImplicitItsStatistics;

    int m_nPadding = 1;

    // Flag to control the update of preconditioning matrix.
    int m_TotNewtonIts = 0;
    int m_TotLinIts    = 0;
    int m_TotImpStages = 0;

    /// Estimate the magnitude of each conserved varibles
    Array<OneD, NekDouble> m_magnitdEstimat;

    Array<OneD, Array<OneD, NekDouble>> m_solutionPhys;

    /// coefff of spacial derivatives(rhs or m_F in GLM) in calculating the
    /// residual of the whole equation(used in unsteady time integrations)
    NekDouble m_TimeIntegLambda = 0.0;
    NekDouble m_inArrayNorm     = -1.0;
    NekDouble m_jacobiFreeEps;

    TensorOfArray4D<NekSingle> m_stdSMatDataDBB;
    TensorOfArray5D<NekSingle> m_stdSMatDataDBDB;

    LibUtilities::NekNonlinSysIterSharedPtr m_nonlinsol;

    PreconCfsSharedPtr m_preconCfs;

    // flag to update shock capturing artificial viscosity. this flag is
    // switched on/off in DoImplicitSolve() to ensure that the AV is only
    // updated once every stage in a multi-stage time integration scheme
    bool m_updateShockCaptPhys{true};

    void v_InitObject(bool DeclareFields = true) override;

    void InitialiseNonlinSysSolver();

    void v_DoSolve() override;

    void v_PrintStatusInformation(const int step,
                                  const NekDouble cpuTime) override;

    void v_PrintSummaryStatistics(const NekDouble intTime) override;

    void v_ALEInitObject(
        int spaceDim,
        Array<OneD, MultiRegions::ExpListSharedPtr> &fields) override;

    void NonlinSysEvaluatorCoeff1D(const Array<OneD, const NekDouble> &inarray,
                                   Array<OneD, NekDouble> &out,
                                   const bool &flag);

    void NonlinSysEvaluatorCoeff(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &out, const bool &flag);

    void DoOdeImplicitRhs(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time);

    void DoOdeRhsCoeff(const Array<OneD, const Array<OneD, NekDouble>> &inarray,
                       Array<OneD, Array<OneD, NekDouble>> &outarray,
                       const NekDouble time);

    void DoAdvectionCoeff(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
        const Array<OneD, const Array<OneD, NekDouble>> &pFwd,
        const Array<OneD, const Array<OneD, NekDouble>> &pBwd);

    void DoDiffusionCoeff(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray,
        const Array<OneD, const Array<OneD, NekDouble>> &pFwd,
        const Array<OneD, const Array<OneD, NekDouble>> &pBwd);

    void DoImplicitSolve(
        const Array<OneD, const Array<OneD, NekDouble>> &inpnts,
        Array<OneD, Array<OneD, NekDouble>> &outpnt, const NekDouble time,
        const NekDouble lambda);

    void DoImplicitSolveCoeff(
        const Array<OneD, const Array<OneD, NekDouble>> &inpnts,
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &out, const NekDouble time,
        const NekDouble lambda);

    void MatrixMultiplyMatrixFreeCoeff(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &out, const bool &centralDifferenceFlag);

    void CalcRefValues(const Array<OneD, const NekDouble> &inarray);

    void PreconCoeff(const Array<OneD, NekDouble> &inarray,
                     Array<OneD, NekDouble> &outarray, const bool &flag);

    template <typename DataType, typename TypeNekBlkMatSharedPtr>
    void AddMatNSBlkDiagVol(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        const Array<OneD, const TensorOfArray2D<NekDouble>> &qfield,
        Array<OneD, Array<OneD, TypeNekBlkMatSharedPtr>> &gmtxarray,
        TensorOfArray4D<DataType> &StdMatDataDBB,
        TensorOfArray5D<DataType> &StdMatDataDBDB);

    template <typename DataType>
    void CalcVolJacStdMat(TensorOfArray4D<DataType> &StdMatDataDBB,
                          TensorOfArray5D<DataType> &StdMatDataDBDB);

    template <typename DataType, typename TypeNekBlkMatSharedPtr>
    void AddMatNSBlkDiagBnd(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        TensorOfArray3D<NekDouble> &qfield,
        TensorOfArray2D<TypeNekBlkMatSharedPtr> &gmtxarray,
        Array<OneD, TypeNekBlkMatSharedPtr> &TraceJac,
        Array<OneD, TypeNekBlkMatSharedPtr> &TraceJacDeriv,
        Array<OneD, Array<OneD, DataType>> &TraceJacDerivSign,
        TensorOfArray5D<DataType> &TraceIPSymJacArray);

    template <typename DataType, typename TypeNekBlkMatSharedPtr>
    void ElmtVarInvMtrx(
        Array<OneD, Array<OneD, TypeNekBlkMatSharedPtr>> &gmtxarray,
        TypeNekBlkMatSharedPtr &gmtVar, const DataType &tmpDatatype);

    template <typename DataType, typename TypeNekBlkMatSharedPtr>
    void GetTraceJac(const Array<OneD, const Array<OneD, NekDouble>> &inarray,
                     TensorOfArray3D<NekDouble> &qfield,
                     Array<OneD, TypeNekBlkMatSharedPtr> &TraceJac,
                     Array<OneD, TypeNekBlkMatSharedPtr> &TraceJacDeriv,
                     Array<OneD, Array<OneD, DataType>> &TraceJacDerivSign,
                     TensorOfArray5D<DataType> &TraceIPSymJacArray);

    template <typename DataType, typename TypeNekBlkMatSharedPtr>
    void NumCalcRiemFluxJac(
        const int nConvectiveFields,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, const Array<OneD, NekDouble>> &AdvVel,
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        TensorOfArray3D<NekDouble> &qfield, const NekDouble &time,
        const Array<OneD, const Array<OneD, NekDouble>> &Fwd,
        const Array<OneD, const Array<OneD, NekDouble>> &Bwd,
        TypeNekBlkMatSharedPtr &FJac, TypeNekBlkMatSharedPtr &BJac,
        TensorOfArray5D<DataType> &TraceIPSymJacArray);

    void PointFluxJacobianPoint(const Array<OneD, NekDouble> &Fwd,
                                const Array<OneD, NekDouble> &normals,
                                DNekMatSharedPtr &FJac, const NekDouble efix,
                                const NekDouble fsw);

    template <typename DataType, typename TypeNekBlkMatSharedPtr>
    void TranSamesizeBlkDiagMatIntoArray(const TypeNekBlkMatSharedPtr &BlkMat,
                                         TensorOfArray3D<DataType> &MatArray);

    template <typename DataType, typename TypeNekBlkMatSharedPtr>
    void TransTraceJacMatToArray(
        const Array<OneD, TypeNekBlkMatSharedPtr> &TraceJac,
        TensorOfArray4D<DataType> &TraceJacDerivArray);

    template <typename DataType, typename TypeNekBlkMatSharedPtr>
    void Fill2DArrayOfBlkDiagonalMat(
        Array<OneD, Array<OneD, TypeNekBlkMatSharedPtr>> &gmtxarray,
        const DataType valu);

    template <typename DataType, typename TypeNekBlkMatSharedPtr>
    void Fill1DArrayOfBlkDiagonalMat(
        Array<OneD, TypeNekBlkMatSharedPtr> &gmtxarray, const DataType valu);

    inline void AllocateNekBlkMatDig(SNekBlkMatSharedPtr &mat,
                                     const Array<OneD, unsigned int> nrow,
                                     const Array<OneD, unsigned int> ncol)
    {
        mat =
            MemoryManager<SNekBlkMat>::AllocateSharedPtr(nrow, ncol, eDIAGONAL);
        SNekMatSharedPtr loc_matNvar;
        for (int nelm = 0; nelm < nrow.size(); ++nelm)
        {
            int nrowsVars = nrow[nelm];
            int ncolsVars = ncol[nelm];

            loc_matNvar = MemoryManager<SNekMat>::AllocateSharedPtr(
                nrowsVars, ncolsVars, 0.0);
            mat->SetBlock(nelm, nelm, loc_matNvar);
        }
    }

    void CalcPreconMatBRJCoeff(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, SNekBlkMatSharedPtr>> &gmtxarray,
        SNekBlkMatSharedPtr &gmtVar, Array<OneD, SNekBlkMatSharedPtr> &TraceJac,
        Array<OneD, SNekBlkMatSharedPtr> &TraceJacDeriv,
        Array<OneD, Array<OneD, NekSingle>> &TraceJacDerivSign,
        TensorOfArray4D<NekSingle> &TraceJacArray,
        TensorOfArray4D<NekSingle> &TraceJacDerivArray,
        TensorOfArray5D<NekSingle> &TraceIPSymJacArray);

    template <typename DataType, typename TypeNekBlkMatSharedPtr>
    void MultiplyElmtInvMassPlusSource(
        Array<OneD, Array<OneD, TypeNekBlkMatSharedPtr>> &gmtxarray,
        const NekDouble dtlamda);

    void GetFluxVectorJacDirElmt(
        const int nConvectiveFields, const int nElmtPnt,
        const Array<OneD, const Array<OneD, NekDouble>> &locVars,
        const Array<OneD, NekDouble> &normals, DNekMatSharedPtr &wspMat,
        Array<OneD, Array<OneD, NekDouble>> &PntJacArray);

    void GetFluxVectorJacPoint(const int nConvectiveFields,
                               const Array<OneD, NekDouble> &conservVar,
                               const Array<OneD, NekDouble> &normals,
                               DNekMatSharedPtr &fluxJac);

    void CalcTraceNumericalFlux(
        const int nConvectiveFields, const int nDim, const int nPts,
        const int nTracePts, const NekDouble PenaltyFactor2,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, const Array<OneD, NekDouble>> &AdvVel,
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        const NekDouble time, TensorOfArray3D<NekDouble> &qfield,
        const Array<OneD, const Array<OneD, NekDouble>> &vFwd,
        const Array<OneD, const Array<OneD, NekDouble>> &vBwd,
        const Array<OneD, const TensorOfArray2D<NekDouble>> &qFwd,
        const Array<OneD, const TensorOfArray2D<NekDouble>> &qBwd,
        const Array<OneD, NekDouble> &MuVarTrace,
        Array<OneD, int> &nonZeroIndex,
        Array<OneD, Array<OneD, NekDouble>> &traceflux);

    void MinusDiffusionFluxJacPoint(
        const int nConvectiveFields, const int nElmtPnt,
        const Array<OneD, const Array<OneD, NekDouble>> &locVars,
        const TensorOfArray3D<NekDouble> &locDerv,
        const Array<OneD, NekDouble> &locmu,
        const Array<OneD, NekDouble> &locDmuDT,
        const Array<OneD, NekDouble> &normals, DNekMatSharedPtr &wspMat,
        Array<OneD, Array<OneD, NekDouble>> &PntJacArray)
    {
        v_MinusDiffusionFluxJacPoint(nConvectiveFields, nElmtPnt, locVars,
                                     locDerv, locmu, locDmuDT, normals, wspMat,
                                     PntJacArray);
    }

    void GetFluxDerivJacDirctn(
        const MultiRegions::ExpListSharedPtr &explist,
        const Array<OneD, const Array<OneD, NekDouble>> &normals,
        const int nDervDir,
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        TensorOfArray5D<NekDouble> &ElmtJacArray, const int nFluxDir)
    {
        v_GetFluxDerivJacDirctn(explist, normals, nDervDir, inarray,
                                ElmtJacArray, nFluxDir);
    }

    void GetFluxDerivJacDirctnElmt(
        const int nConvectiveFields, const int nElmtPnt, const int nDervDir,
        const Array<OneD, const Array<OneD, NekDouble>> &locVars,
        const Array<OneD, NekDouble> &locmu,
        const Array<OneD, const Array<OneD, NekDouble>> &locnormal,
        DNekMatSharedPtr &wspMat,
        Array<OneD, Array<OneD, NekDouble>> &PntJacArray)
    {
        v_GetFluxDerivJacDirctnElmt(nConvectiveFields, nElmtPnt, nDervDir,
                                    locVars, locmu, locnormal, wspMat,
                                    PntJacArray);
    }

    void GetFluxDerivJacDirctn(
        const MultiRegions::ExpListSharedPtr &explist,
        const Array<OneD, const Array<OneD, NekDouble>> &normals,
        const int nDervDir,
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, DNekMatSharedPtr>> &ElmtJac)
    {
        v_GetFluxDerivJacDirctn(explist, normals, nDervDir, inarray, ElmtJac);
    }

    void CalcPhysDeriv(const Array<OneD, const Array<OneD, NekDouble>> &inarray,
                       TensorOfArray3D<NekDouble> &qfield)
    {
        v_CalcPhysDeriv(inarray, qfield);
    }

    void CalcMuDmuDT(const Array<OneD, const Array<OneD, NekDouble>> &inarray,
                     Array<OneD, NekDouble> &mu, Array<OneD, NekDouble> &DmuDT)
    {
        v_CalcMuDmuDT(inarray, mu, DmuDT);
    }

    virtual void v_DoDiffusionCoeff(
        [[maybe_unused]] const Array<OneD, const Array<OneD, NekDouble>>
            &inarray,
        [[maybe_unused]] Array<OneD, Array<OneD, NekDouble>> &outarray,
        [[maybe_unused]] const Array<OneD, const Array<OneD, NekDouble>> &pFwd,
        [[maybe_unused]] const Array<OneD, const Array<OneD, NekDouble>> &pBwd)
    {
    }

    virtual void v_CalcMuDmuDT(
        [[maybe_unused]] const Array<OneD, const Array<OneD, NekDouble>>
            &inarray,
        [[maybe_unused]] Array<OneD, NekDouble> &mu,
        [[maybe_unused]] Array<OneD, NekDouble> &DmuDT)
    {
    }

    virtual void v_CalcPhysDeriv(
        [[maybe_unused]] const Array<OneD, const Array<OneD, NekDouble>>
            &inarray,
        [[maybe_unused]] TensorOfArray3D<NekDouble> &qfield)
    {
    }

    virtual void v_MinusDiffusionFluxJacPoint(
        const int nConvectiveFields, const int nElmtPnt,
        const Array<OneD, const Array<OneD, NekDouble>> &locVars,
        const TensorOfArray3D<NekDouble> &locDerv,
        const Array<OneD, NekDouble> &locmu,
        const Array<OneD, NekDouble> &locDmuDT,
        const Array<OneD, NekDouble> &normals, DNekMatSharedPtr &wspMat,
        Array<OneD, Array<OneD, NekDouble>> &PntJacArray);

    virtual void v_GetFluxDerivJacDirctn(
        const MultiRegions::ExpListSharedPtr &explist,
        const Array<OneD, const Array<OneD, NekDouble>> &normals,
        const int nDervDir,
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        TensorOfArray5D<NekDouble> &ElmtJacArray, const int nFluxDir);

    virtual void v_GetFluxDerivJacDirctnElmt(
        const int nConvectiveFields, const int nElmtPnt, const int nDervDir,
        const Array<OneD, const Array<OneD, NekDouble>> &locVars,
        const Array<OneD, NekDouble> &locmu,
        const Array<OneD, const Array<OneD, NekDouble>> &locnormal,
        DNekMatSharedPtr &wspMat,
        Array<OneD, Array<OneD, NekDouble>> &PntJacArray);

    virtual void v_GetFluxDerivJacDirctn(
        const MultiRegions::ExpListSharedPtr &explist,
        const Array<OneD, const Array<OneD, NekDouble>> &normals,
        const int nDervDir,
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, DNekMatSharedPtr>> &ElmtJac);

    bool v_UpdateTimeStepCheck() override;
};
} // namespace Nektar
#endif
