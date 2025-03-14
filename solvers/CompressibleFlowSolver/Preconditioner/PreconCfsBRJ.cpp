///////////////////////////////////////////////////////////////////////////////
//
// File:  PreconCfsBRJ.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
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
// Description:  PreconCfsBRJ definition
//
///////////////////////////////////////////////////////////////////////////////

#include <CompressibleFlowSolver/Preconditioner/PreconCfsBRJ.h>
#include <LibUtilities/BasicUtils/Timer.h>

using namespace std;

namespace Nektar
{
/**
 * @class  PreconCfsBRJ
 *
 * Solves a linear system using iterative methods.
 */
std::string PreconCfsBRJ::className =
    GetPreconCfsFactory().RegisterCreatorFunction(
        "PreconCfsBRJ", PreconCfsBRJ::create,
        "Block Relaxed Jacobi Preconditioner for CFS.");

PreconCfsBRJ::PreconCfsBRJ(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const LibUtilities::CommSharedPtr &vComm)
    : PreconCfs(pFields, pSession, vComm)
{
    pSession->LoadParameter("PreconItsStep", m_PreconItsStep, 7);
    pSession->LoadParameter("BRJRelaxParam", m_BRJRelaxParam, 1.0);

    size_t nvariables = pFields.size();

    m_PreconMatVarsSingle =
        Array<OneD, Array<OneD, SNekBlkMatSharedPtr>>(nvariables);
    for (size_t i = 0; i < nvariables; i++)
    {
        m_PreconMatVarsSingle[i] = Array<OneD, SNekBlkMatSharedPtr>(nvariables);
    }
    AllocatePreconBlkDiagCoeff(pFields, m_PreconMatVarsSingle);

    AllocateSIMDPreconBlkMatDiag(pFields);
}

/**
 *
 */
void PreconCfsBRJ::v_InitObject()
{
}

/**
 *
 */
void PreconCfsBRJ::v_DoPreconCfs(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    const Array<OneD, NekDouble> &inarray, Array<OneD, NekDouble> &outarray,
    [[maybe_unused]] const bool &flag)
{
    size_t nBRJIterTot = m_PreconItsStep;
    if (nBRJIterTot == 0)
    {
        DoNullPrecon(inarray, outarray, flag);
    }
    else
    {
        const NekDouble BRJParam   = m_BRJRelaxParam;
        const NekDouble OmBRJParam = 1.0 - BRJParam;

        size_t nvariables = pFields.size();
        size_t npoints    = pFields[0]->GetNcoeffs();
        size_t ntotpnt    = inarray.size();

        ASSERTL0(nvariables * npoints == ntotpnt,
                 "nvariables*npoints!=ntotpnt in PreconCoeff");

        Array<OneD, NekDouble> rhs(ntotpnt);

        Array<OneD, NekDouble> outN(ntotpnt);
        Array<OneD, NekDouble> outTmp(ntotpnt);
        Array<OneD, Array<OneD, NekDouble>> rhs2d(nvariables);
        Array<OneD, Array<OneD, NekDouble>> out_2d(nvariables);
        Array<OneD, Array<OneD, NekDouble>> outTmp_2d(nvariables);
        for (size_t m = 0; m < nvariables; m++)
        {
            size_t moffset = m * npoints;
            rhs2d[m]       = rhs + moffset;
            out_2d[m]      = outarray + moffset;
            outTmp_2d[m]   = outTmp + moffset;
            pFields[m]->MultiplyByMassMatrix(inarray + moffset, rhs2d[m]);
        }

        size_t nphysic   = pFields[0]->GetNpoints();
        size_t nTracePts = pFields[0]->GetTrace()->GetNpoints();
        TensorOfArray3D<NekDouble> qfield(m_spacedim);
        for (int i = 0; i < m_spacedim; i++)
        {
            qfield[i] = Array<OneD, Array<OneD, NekDouble>>(nvariables);
            for (size_t j = 0; j < nvariables; j++)
            {
                qfield[i][j] = Array<OneD, NekDouble>(nphysic);
            }
        }
        size_t ntmpTrace = 4 + 2 * m_spacedim;
        TensorOfArray3D<NekDouble> tmpTrace(ntmpTrace);
        for (size_t i = 0; i < ntmpTrace; i++)
        {
            tmpTrace[i] = Array<OneD, Array<OneD, NekDouble>>(nvariables);
            for (size_t j = 0; j < nvariables; j++)
            {
                tmpTrace[i][j] = Array<OneD, NekDouble>(nTracePts);
            }
        }
        Array<OneD, Array<OneD, NekDouble>> FwdFluxDeriv(nvariables);
        Array<OneD, Array<OneD, NekDouble>> BwdFluxDeriv(nvariables);
        for (size_t j = 0; j < nvariables; j++)
        {
            FwdFluxDeriv[j] = Array<OneD, NekDouble>(nTracePts);
            BwdFluxDeriv[j] = Array<OneD, NekDouble>(nTracePts);
        }

        const size_t nwspTraceDataType = nvariables + 1;
        Array<OneD, Array<OneD, NekSingle>> wspTraceDataType(nwspTraceDataType);
        for (size_t m = 0; m < nwspTraceDataType; m++)
        {
            wspTraceDataType[m] = Array<OneD, NekSingle>(nTracePts);
        }

        LibUtilities::Timer timer;
        timer.Start();
        PreconBlkDiag(pFields, rhs, outarray);
        timer.Stop();
        timer.AccumulateRegion("PreconCfsBRJ::PreconBlkDiag", 2);

        for (size_t nrelax = 0; nrelax < nBRJIterTot - 1; nrelax++)
        {
            Vmath::Smul(ntotpnt, OmBRJParam, outarray, 1, outN, 1);

            timer.Start();
            MinusOffDiag2Rhs(pFields, nvariables, npoints, rhs2d, out_2d,
                             tmpTrace, wspTraceDataType, m_TraceJacArraySingle);
            timer.Stop();
            timer.AccumulateRegion("PreconCfsBRJ::MinusOffDiag2Rhs", 2);

            timer.Start();
            PreconBlkDiag(pFields, outarray, outTmp);
            timer.Stop();
            timer.AccumulateRegion("PreconCfsBRJ::PreconBlkDiag", 2);

            Vmath::Svtvp(ntotpnt, BRJParam, outTmp, 1, outN, 1, outarray, 1);
        }
    }
}

/**
 *
 */
void PreconCfsBRJ::v_BuildPreconCfs(
    [[maybe_unused]] const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    const Array<OneD, const Array<OneD, NekDouble>> &intmp,
    [[maybe_unused]] const NekDouble time, const NekDouble lambda)
{
    if (m_PreconItsStep > 0)
    {
        SNekBlkMatSharedPtr PreconMatSingle;
        using vec_t    = simd<NekSingle>;
        int nvariables = pFields.size();
        int nelmts     = pFields[0]->GetNumElmts();
        Array<OneD, unsigned int> matdim(nelmts);
        for (int i = 0; i < nelmts; i++)
        {
            matdim[i] = pFields[0]->GetExp(i)->GetNcoeffs() * nvariables;
        }
        AllocateNekBlkMatDig(PreconMatSingle, matdim, matdim);

        m_operator.DoCalcPreconMatBRJCoeff(
            intmp, m_PreconMatVarsSingle, PreconMatSingle, m_TraceJacSingle,
            m_TraceJacDerivSingle, m_TraceJacDerivSignSingle,
            m_TraceJacArraySingle, m_TraceJacDerivArraySingle,
            m_TraceIPSymJacArraySingle);

        if (m_verbose && m_Comm->GetRank() == 0)
        {
            cout << "     ## CalcuPreconMat " << endl;
        }

        // copy matrix to simd layout
        // load matrix
        int cnt             = 0;
        const auto vecwidth = vec_t::width;

        alignas(vec_t::alignment) std::array<NekSingle, vec_t::width> tmp;

        for (int ne = 0; ne < nelmts; ne++)
        {
            const auto nElmtDof = matdim[ne];
            const auto nblocks  = nElmtDof / vecwidth;

            const NekSingle *mmat =
                PreconMatSingle->GetBlockPtr(ne, ne)->GetRawPtr();
            /// Copy array into column major blocks of vector width
            for (int i1 = 0; i1 < nblocks; ++i1)
            {
                for (int j = 0; j < nElmtDof; ++j)
                {
                    for (int i = 0; i < vecwidth; ++i)
                    {
                        tmp[i] = mmat[j + (i1 * vecwidth + i) * nElmtDof];
                    }
                    // store contiguous vec_t array.
                    m_sBlkDiagMat[cnt++].load(tmp.data());
                }
            }

            const auto endwidth = nElmtDof - nblocks * vecwidth;

            // end rows that do not fit into vector widths
            if (endwidth)
            {
                for (int j = 0; j < nElmtDof; ++j)
                {
                    for (int i = 0; i < endwidth; ++i)
                    {
                        tmp[i] = mmat[j + (nblocks * vecwidth + i) * nElmtDof];
                    }

                    for (int i = endwidth; i < vecwidth; ++i)
                    {
                        tmp[i] = 0.0;
                    }
                    m_sBlkDiagMat[cnt++].load(tmp.data());
                }
            }
        }
    }

    m_DtLambdaPreconMat  = lambda;
    m_CalcPreconMatFlag  = false;
    m_PreconTimesCounter = 1;
}

/**
 *
 */
bool PreconCfsBRJ::v_UpdatePreconMatCheck(
    [[maybe_unused]] const Array<OneD, const NekDouble> &res,
    const NekDouble dtLambda)
{
    bool flag = (m_CalcPreconMatFlag || m_DtLambdaPreconMat != dtLambda ||
                 m_PreconMatFreezNumb < m_PreconTimesCounter);
    m_CalcPreconMatFlag = flag;
    return flag;
}

/**
 *
 */
void PreconCfsBRJ::DoNullPrecon(const Array<OneD, NekDouble> &pInput,
                                Array<OneD, NekDouble> &pOutput,
                                [[maybe_unused]] const bool &flag)
{
    Vmath::Vcopy(pInput.size(), pInput, 1, pOutput, 1);
}

/**
 *
 */
void PreconCfsBRJ::PreconBlkDiag(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    const Array<OneD, NekDouble> &inarray, Array<OneD, NekDouble> &outarray)
{
    unsigned int nvariables = pFields.size();

    int nTotElmt = pFields[0]->GetNumElmts();

    using vec_t         = simd<NekSingle>;
    const auto vecwidth = vec_t::width;

    // vectorized matrix multiply
    std::vector<vec_t, tinysimd::allocator<vec_t>> Sinarray(m_max_nblocks);
    std::vector<vec_t, tinysimd::allocator<vec_t>> Soutarray(m_max_nElmtDof);

    alignas(vec_t::alignment) std::array<NekSingle, vec_t::width> tmp;

    for (int ne = 0, cnt = 0, icnt = 0, icnt1 = 0; ne < nTotElmt; ne++)
    {
        const auto nElmtCoef = pFields[0]->GetNcoeffs(ne);
        const auto nElmtDof  = nElmtCoef * nvariables;
        const auto nblocks   = (nElmtDof % vecwidth) ? nElmtDof / vecwidth + 1
                                                     : nElmtDof / vecwidth;

        // gather data into blocks - could probably be done with a
        // gather call? can be replaced with a gather op when working
        for (int j = 0; j < nblocks; ++j, icnt += vecwidth)
        {
            for (int i = 0; i < vecwidth; ++i)
            {
                tmp[i] = inarray[m_inputIdx[icnt + i]];
            }

            Sinarray[j].load(tmp.data());
        }

        // Do matrix multiply
        // first block just needs multiplying
        vec_t in = Sinarray[0];
        for (int i = 0; i < nElmtDof; ++i)
        {
            Soutarray[i] = m_sBlkDiagMat[cnt++] * in;
        }

        // rest of blocks are multiply add operations;
        for (int n = 1; n < nblocks; ++n)
        {
            in = Sinarray[n];
            for (int i = 0; i < nElmtDof; ++i)
            {
                Soutarray[i].fma(m_sBlkDiagMat[cnt++], in);
            }
        }

        // get block aligned index for this expansion
        NekSingle val;
        for (int i = 0; i < nElmtDof; ++i)
        {
            // Get hold of data
            Soutarray[i].store(tmp.data());

            // Sum vector width
            val = tmp[0];
            for (int j = 1; j < vecwidth; ++j)
            {
                val += tmp[j];
            }
            // put data into outarray
            outarray[m_inputIdx[icnt1 + i]] = NekDouble(val);
        }

        icnt1 += nblocks * vecwidth;
    }
}

/**
 *
 */
template <typename DataType>
void PreconCfsBRJ::MinusOffDiag2Rhs(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    const size_t nvariables, const size_t nCoeffs,
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray,
    TensorOfArray3D<NekDouble> &wspTrace,
    Array<OneD, Array<OneD, DataType>> &wspTraceDataType,
    const TensorOfArray4D<DataType> &TraceJacArray)
{
    size_t nTracePts = pFields[0]->GetTrace()->GetNpoints();
    size_t npoints   = pFields[0]->GetNpoints();
    size_t nDim      = m_spacedim;

    Array<OneD, Array<OneD, NekDouble>> outpnts(nvariables);
    for (size_t i = 0; i < nvariables; i++)
    {
        outpnts[i] = Array<OneD, NekDouble>(npoints, 0.0);
        pFields[i]->BwdTrans(outarray[i], outpnts[i]);
    }

    // Store forwards/backwards space along trace space
    Array<OneD, Array<OneD, NekDouble>> Fwd;
    Array<OneD, Array<OneD, NekDouble>> Bwd;
    Array<OneD, Array<OneD, NekDouble>> FwdFlux;
    Array<OneD, Array<OneD, NekDouble>> BwdFlux;
    TensorOfArray3D<NekDouble> numDerivBwd(nDim);
    TensorOfArray3D<NekDouble> numDerivFwd(nDim);
    size_t indexwspTrace = 0;
    Fwd                  = wspTrace[indexwspTrace], indexwspTrace++;
    Bwd                  = wspTrace[indexwspTrace], indexwspTrace++;
    FwdFlux              = wspTrace[indexwspTrace], indexwspTrace++;
    BwdFlux              = wspTrace[indexwspTrace], indexwspTrace++;

    LibUtilities::Timer timer;
    for (size_t i = 0; i < nvariables; ++i)
    {
        timer.Start();
        pFields[i]->GetFwdBwdTracePhys(outpnts[i], Fwd[i], Bwd[i]);
        timer.Stop();
        timer.AccumulateRegion("ExpList::GetFwdBwdTracePhys", 10);
    }

    size_t indexwspTraceDataType = 0;
    Array<OneD, Array<OneD, DataType>> Fwdarray(nvariables);
    for (size_t m = 0; m < nvariables; ++m)
    {
        Fwdarray[m] = wspTraceDataType[indexwspTraceDataType],
        indexwspTraceDataType++;
    }
    Array<OneD, DataType> Fwdreslt;
    Fwdreslt = wspTraceDataType[indexwspTraceDataType], indexwspTraceDataType++;

    for (size_t m = 0; m < nvariables; ++m)
    {
        for (size_t i = 0; i < nTracePts; ++i)
        {
            Fwdarray[m][i] = DataType(Fwd[m][i]);
        }
    }
    for (size_t m = 0; m < nvariables; ++m)
    {
        Vmath::Zero(nTracePts, Fwdreslt, 1);
        for (size_t n = 0; n < nvariables; ++n)
        {
            for (size_t p = 0; p < nTracePts; ++p)
            {
                Fwdreslt[p] += TraceJacArray[0][m][n][p] * Fwdarray[n][p];
            }
        }
        for (size_t i = 0; i < nTracePts; ++i)
        {
            FwdFlux[m][i] = NekDouble(Fwdreslt[i]);
        }
    }

    for (size_t m = 0; m < nvariables; ++m)
    {
        for (size_t i = 0; i < nTracePts; ++i)
        {
            Fwdarray[m][i] = DataType(Bwd[m][i]);
        }
    }
    for (size_t m = 0; m < nvariables; ++m)
    {
        Vmath::Zero(nTracePts, Fwdreslt, 1);
        for (size_t n = 0; n < nvariables; ++n)
        {
            for (size_t p = 0; p < nTracePts; ++p)
            {
                Fwdreslt[p] += TraceJacArray[1][m][n][p] * Fwdarray[n][p];
            }
        }
        for (size_t i = 0; i < nTracePts; ++i)
        {
            BwdFlux[m][i] = NekDouble(Fwdreslt[i]);
        }
    }

    for (size_t i = 0; i < nvariables; ++i)
    {
        Vmath::Fill(nCoeffs, 0.0, outarray[i], 1);
        timer.Start();
        pFields[i]->AddTraceIntegralToOffDiag(FwdFlux[i], BwdFlux[i],
                                              outarray[i]);
        timer.Stop();
        timer.AccumulateRegion("ExpList::AddTraceIntegralToOffDiag", 10);
    }

    for (size_t i = 0; i < nvariables; ++i)
    {
        for (size_t p = 0; p < nCoeffs; ++p)
        {
            outarray[i][p] =
                -m_DtLambdaPreconMat * outarray[i][p] + inarray[i][p];
        }
    }
}

/**
 *
 */
template <typename TypeNekBlkMatSharedPtr>
void PreconCfsBRJ::AllocatePreconBlkDiagCoeff(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    Array<OneD, Array<OneD, TypeNekBlkMatSharedPtr>> &gmtxarray,
    const int &nscale)
{

    size_t nvars  = pFields.size();
    size_t nelmts = pFields[0]->GetNumElmts();
    size_t nelmtcoef;
    Array<OneD, unsigned int> nelmtmatdim(nelmts);
    for (size_t i = 0; i < nelmts; i++)
    {
        nelmtcoef      = pFields[0]->GetExp(i)->GetNcoeffs();
        nelmtmatdim[i] = nelmtcoef * nscale;
    }

    for (size_t i = 0; i < nvars; i++)
    {
        for (size_t j = 0; j < nvars; j++)
        {
            AllocateNekBlkMatDig(gmtxarray[i][j], nelmtmatdim, nelmtmatdim);
        }
    }
}

} // namespace Nektar
