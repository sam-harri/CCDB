///////////////////////////////////////////////////////////////////////////////
//
// File: NekLinSysIterGMRESLoc.cpp
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
// Description: NekLinSysIterGMRESLoc definition
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/Timer.h>
#include <LibUtilities/LinearAlgebra/NekLinSysIterGMRESLoc.h>

using namespace std;

namespace Nektar::LibUtilities
{
/**
 * @class  NekLinSysIterGMRESLoc
 *
 * Solves a linear system using iterative methods using local storage rather
 * than global
 */
string NekLinSysIterGMRESLoc::className =
    LibUtilities::GetNekLinSysIterFactory().RegisterCreatorFunction(
        "GMRESLoc", NekLinSysIterGMRESLoc::create,
        "NekLinSysIterGMRES local storage solver.");

NekLinSysIterGMRESLoc::NekLinSysIterGMRESLoc(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const LibUtilities::CommSharedPtr &vRowComm, const int nDimen,
    const NekSysKey &pKey)
    : NekLinSysIter(pSession, vRowComm, nDimen, pKey)
{
    m_NekLinSysLeftPrecon  = pKey.m_NekLinSysLeftPrecon;
    m_NekLinSysRightPrecon = pKey.m_NekLinSysRightPrecon;

    m_KrylovMaxHessMatBand = pKey.m_KrylovMaxHessMatBand;

    m_maxrestart       = ceil(NekDouble(m_NekLinSysMaxIterations) /
                              NekDouble(pKey.m_LinSysMaxStorage));
    m_LinSysMaxStorage = min(m_NekLinSysMaxIterations, pKey.m_LinSysMaxStorage);

    m_GMRESCentralDifference = pKey.m_GMRESCentralDifference;

    m_isLocal = true;

    // Allocate array storage of coefficients
    // Hessenburg matrix
    m_hes = Array<OneD, Array<OneD, NekDouble>>(m_LinSysMaxStorage);
    for (size_t nd = 0; nd < m_LinSysMaxStorage; nd++)
    {
        m_hes[nd] = Array<OneD, NekDouble>(m_LinSysMaxStorage + 1, 0.0);
    }
    // Hesseburg matrix after rotation
    m_Upper = Array<OneD, Array<OneD, NekDouble>>(m_LinSysMaxStorage);
    for (size_t nd = 0; nd < m_LinSysMaxStorage; nd++)
    {
        m_Upper[nd] = Array<OneD, NekDouble>(m_LinSysMaxStorage + 1, 0.0);
    }
    // Total search directions
    m_V_total = Array<OneD, Array<OneD, NekDouble>>(m_LinSysMaxStorage + 1);
}

void NekLinSysIterGMRESLoc::v_InitObject()
{
    NekLinSysIter::v_InitObject();
}

/**
 *
 */
int NekLinSysIterGMRESLoc::v_SolveSystem(
    const int nLocal, const Array<OneD, const NekDouble> &pInput,
    Array<OneD, NekDouble> &pOutput, [[maybe_unused]] const int nDir)
{
    int niterations = DoGMRES(nLocal, pInput, pOutput);

    return niterations;
}

/**  
 * Solve a global linear system using the Gmres 
 * We solve only for the non-Dirichlet modes. The operator is evaluated  
 * using an auxiliary function v_DoMatrixMultiply defined by the  
 * specific solver. Distributed math routines are used to support  
 * parallel execution of the solver.  
 *  
 * @param       pInput      Input residual  in local format.
 * @param       pOutput     Solution vector in local format but with continuous
 * values
 */

int NekLinSysIterGMRESLoc::DoGMRES(const int nLocal,
                                   const Array<OneD, const NekDouble> &pInput,
                                   Array<OneD, NekDouble> &pOutput)
{
    m_prec_factor = NekConstants::kNekUnsetDouble;
    if (m_rhs_magnitude == NekConstants::kNekUnsetDouble)
    {
        Set_Rhs_Magnitude(pInput);
    }

    NekDouble eps = 0.0;
    Array<OneD, NekDouble> tmp;

    m_totalIterations = 0;
    m_converged       = false;

    bool restarted = false;
    bool truncted  = false;

    if (m_KrylovMaxHessMatBand > 0)
    {
        truncted = true;
    }

    for (int nrestart = 0; nrestart < m_maxrestart; ++nrestart)
    {
        eps = DoGmresRestart(restarted, truncted, nLocal, pInput, pOutput);

        if (m_converged)
        {
            break;
        }
        restarted = true;
    }

    if (m_verbose)
    {
        Array<OneD, NekDouble> r0(nLocal);
        Array<OneD, NekDouble> wk(nLocal);

        // calculate difference in residual of solution
        m_operator.DoNekSysLhsEval(pOutput, r0, m_GMRESCentralDifference);

        // Note this is finding the difference between the whole
        // residual not jsut the non-Dirichlet values.
        // Probably OK since just an monitoring output?
        Vmath::Vsub(nLocal, pInput, 1, r0, 1, r0, 1);

        m_operator.DoAssembleLoc(r0, wk, true);
        NekDouble vExchange = Vmath::Dot(nLocal, wk, r0);

        m_rowComm->AllReduce(vExchange, LibUtilities::ReduceSum);
        NekDouble eps1 = vExchange;

        if (m_root)
        {
            int nwidthcolm = 13;

            cout << std::scientific << std::setw(nwidthcolm)
                 << std::setprecision(nwidthcolm - 8)
                 << "       GMRES iterations made = " << m_totalIterations
                 << " using tolerance of " << m_NekLinSysTolerance
                 << " (error = " << sqrt(eps * m_prec_factor / m_rhs_magnitude)
                 << ")";

            cout << " WITH (GMRES eps = " << eps << " REAL eps= " << eps1
                 << ")";

            if (m_converged)
            {
                cout << " CONVERGED" << endl;
            }
            else
            {
                cout << " WARNING: Exceeded maxIt" << endl;
            }
        }
    }

    if (m_FlagWarnings)
    {
        WARNINGL1(m_converged, "GMRES did not converge.");
    }
    return m_totalIterations;
}

NekDouble NekLinSysIterGMRESLoc::DoGmresRestart(
    const bool restarted, const bool truncted, const int nLocal,
    const Array<OneD, const NekDouble> &pInput, Array<OneD, NekDouble> &pOutput)
{
    // Allocate array storage of coefficients
    // Residual
    Array<OneD, NekDouble> eta(m_LinSysMaxStorage + 1, 0.0);
    // Givens rotation c
    Array<OneD, NekDouble> cs(m_LinSysMaxStorage, 0.0);
    // Givens rotation s
    Array<OneD, NekDouble> sn(m_LinSysMaxStorage, 0.0);
    // Total coefficients, just for check
    Array<OneD, NekDouble> y_total(m_LinSysMaxStorage, 0.0);
    // Search direction order
    Array<OneD, int> id(m_LinSysMaxStorage, 0);
    Array<OneD, int> id_start(m_LinSysMaxStorage, 0);
    Array<OneD, int> id_end(m_LinSysMaxStorage, 0);
    // Temporary variables
    int idtem;
    int starttem;
    int endtem;

    NekDouble eps;
    NekDouble beta, alpha;
    NekDouble vExchange = 0;
    // Temporary Array
    Array<OneD, NekDouble> w(nLocal, 0.0);
    Array<OneD, NekDouble> wk(nLocal, 0.0);
    Array<OneD, NekDouble> r0(nLocal, 0.0);
    Array<OneD, NekDouble> V1;
    Array<OneD, NekDouble> V2;
    Array<OneD, NekDouble> h1;
    Array<OneD, NekDouble> h2;

    if (restarted)
    {
        // This is A*x
        m_operator.DoNekSysLhsEval(pOutput, r0, m_GMRESCentralDifference);

        // The first search direction
        beta = -1.0;

        // This is r0 = b-A*x
        Vmath::Svtvp(nLocal, beta, r0, 1, pInput, 1, r0, 1);
    }
    else
    {
        // If not restarted, x0 should be zero
        Vmath::Vcopy(nLocal, pInput, 1, r0, 1);
    }

    if (m_NekLinSysLeftPrecon)
    {
        m_operator.DoNekSysPrecon(r0, r0);
    }

    // Norm of (r0)
    // The m_map tells how to connect
    m_operator.DoAssembleLoc(r0, wk, true);
    vExchange = Vmath::Dot(nLocal, wk, r0);
    m_rowComm->AllReduce(vExchange, LibUtilities::ReduceSum);
    eps = vExchange;

    if (!restarted)
    {
        if (m_prec_factor == NekConstants::kNekUnsetDouble)
        {
            if (m_NekLinSysLeftPrecon)
            {
                // Evaluate initial residual error for exit check
                ASSERTL0(false, "Need to set up/debugging");

                m_operator.DoAssembleLoc(pInput, wk, true);
                vExchange = Vmath::Dot(nLocal, wk, pInput);
                m_rowComm->AllReduce(vExchange, LibUtilities::ReduceSum);
                m_prec_factor = vExchange / eps;
            }
            else
            {
                m_prec_factor = 1.0;
            }
        }
    }

    Vmath::Smul(nLocal, sqrt(m_prec_factor), r0, 1, r0, 1);
    eps    = eps * m_prec_factor;
    eta[0] = sqrt(eps);

    // Give an order for the entries in Hessenburg matrix
    for (int nd = 0; nd < m_LinSysMaxStorage; ++nd)
    {
        id[nd]     = nd;
        id_end[nd] = nd + 1;
        starttem   = id_end[nd] - m_KrylovMaxHessMatBand;
        if (truncted && (starttem) > 0)
        {
            id_start[nd] = starttem;
        }
        else
        {
            id_start[nd] = 0;
        }
    }

    // Normlized by r0 norm V(:,1)=r0/norm(r0)
    alpha = 1.0 / eta[0];

    // Scalar multiplication
    if (m_V_total[0].size() == 0)
    {
        m_V_total[0] = Array<OneD, NekDouble>(nLocal, 0.0);
    }
    Vmath::Smul(nLocal, alpha, r0, 1, m_V_total[0], 1);

    // restarted Gmres(m) process
    if (m_NekLinSysRightPrecon)
    {
        V1 = Array<OneD, NekDouble>(nLocal, 0.0);
    }

    int nswp = 0;
    for (int nd = 0; nd < m_LinSysMaxStorage; ++nd)
    {
        if (m_V_total[nd + 1].size() == 0)
        {
            m_V_total[nd + 1] = Array<OneD, NekDouble>(nLocal, 0.0);
        }
        Vmath::Zero(nLocal, m_V_total[nd + 1], 1);
        Vmath::Zero(m_LinSysMaxStorage + 1, m_hes[nd], 1);
        V2 = m_V_total[nd + 1];
        h1 = m_hes[nd];

        if (m_NekLinSysRightPrecon)
        {
            m_operator.DoNekSysPrecon(m_V_total[nd], V1, true);
        }
        else
        {
            V1 = m_V_total[nd];
        }

        // w here is no need to add nDir due to temporary Array
        idtem    = id[nd];
        starttem = id_start[idtem];
        endtem   = id_end[idtem];

        DoArnoldi(starttem, endtem, nLocal, w, wk, V1, V2, h1);

        if (starttem > 0)
        {
            starttem = starttem - 1;
        }

        h2 = m_Upper[nd];
        Vmath::Vcopy(m_LinSysMaxStorage + 1, &h1[0], 1, &h2[0], 1);
        DoGivensRotation(starttem, endtem, cs, sn, h2, eta);

        eps = eta[nd + 1] * eta[nd + 1];

        // This Gmres merge truncted Gmres to accelerate.
        // If truncted, cannot jump out because
        // the last term of eta is not residual
        if ((!truncted) || (nd < m_KrylovMaxHessMatBand))
        {
            if ((eps < m_NekLinSysTolerance * m_NekLinSysTolerance *
                           m_rhs_magnitude)) //&& nd > 0)
            {
                m_converged = true;
            }
        }
        nswp++;
        m_totalIterations++;

        if (m_converged)
        {
            break;
        }
    }

    DoBackward(nswp, m_Upper, eta, y_total);

    // calculate output y_total*V_total
    Array<OneD, NekDouble> solution(nLocal, 0.0);
    for (int i = 0; i < nswp; ++i)
    {
        beta = y_total[i];
        Vmath::Svtvp(nLocal, beta, m_V_total[i], 1, solution, 1, solution, 1);
    }

    if (m_NekLinSysRightPrecon)
    {
        m_operator.DoNekSysPrecon(solution, solution, true);
    }

    // Update output.
    Vmath::Vadd(nLocal, solution, 1, pOutput, 1, pOutput, 1);

    return eps;
}

// Arnoldi Subroutine
void NekLinSysIterGMRESLoc::DoArnoldi(const int starttem, const int endtem,
                                      const int nLocal,
                                      Array<OneD, NekDouble> &w,
                                      Array<OneD, NekDouble> &wk,
                                      Array<OneD, NekDouble> &V1,
                                      Array<OneD, NekDouble> &V2,
                                      Array<OneD, NekDouble> &h)
{
    NekDouble alpha, beta, vExchange = 0.0;
    LibUtilities::Timer timer;
    timer.Start();
    m_operator.DoNekSysLhsEval(V1, w, m_GMRESCentralDifference);
    timer.Stop();
    timer.AccumulateRegion("NekSysOperators::DoNekSysLhsEval", 10);

    if (m_NekLinSysLeftPrecon)
    {
        m_operator.DoNekSysPrecon(w, w);
    }

    Vmath::Smul(nLocal, sqrt(m_prec_factor), w, 1, w, 1);

    // Modified Gram-Schmidt
    for (int i = starttem; i < endtem; ++i)
    {
        // Do inner product on equivalent of global values excluding
        // Diriclet conditions. To do this need to assmble (and
        // scatter back vector and then zero Dirichlet conditions.
        m_operator.DoAssembleLoc(m_V_total[i], wk, true);
        vExchange = Vmath::Dot(nLocal, wk, w);
        m_rowComm->AllReduce(vExchange, LibUtilities::ReduceSum);

        h[i] = vExchange;

        beta = -1.0 * vExchange;
        Vmath::Svtvp(nLocal, beta, m_V_total[i], 1, w, 1, w, 1);
    }
    // end of Modified Gram-Schmidt

    // calculate the L2 norm and normalize
    m_operator.DoAssembleLoc(w, wk, true);
    vExchange = Vmath::Dot(nLocal, wk, w);
    m_rowComm->AllReduce(vExchange, LibUtilities::ReduceSum);

    h[endtem] = sqrt(vExchange);

    alpha = 1.0 / h[endtem];
    Vmath::Smul(nLocal, alpha, w, 1, V2, 1);
}

// QR factorization through Givens rotation -> Put into a helper class
void NekLinSysIterGMRESLoc::DoGivensRotation(const int starttem,
                                             const int endtem,
                                             Array<OneD, NekDouble> &c,
                                             Array<OneD, NekDouble> &s,
                                             Array<OneD, NekDouble> &h,
                                             Array<OneD, NekDouble> &eta)
{
    NekDouble temp_dbl;
    NekDouble dd;
    NekDouble hh;
    int idtem = endtem - 1;
    // The starttem and endtem are beginning and ending order of Givens rotation
    // They usually equal to the beginning position and ending position of
    // Hessenburg matrix But sometimes starttem will change, like if it is
    // initial 0 and becomes nonzero because previous Givens rotation See Yu
    // Pan's User Guide
    for (int i = starttem; i < idtem; ++i)
    {
        temp_dbl = c[i] * h[i] - s[i] * h[i + 1];
        h[i + 1] = s[i] * h[i] + c[i] * h[i + 1];
        h[i]     = temp_dbl;
    }
    dd = h[idtem];
    hh = h[endtem];
    if (hh == 0.0)
    {
        c[idtem] = 1.0;
        s[idtem] = 0.0;
    }
    else if (abs(hh) > abs(dd))
    {
        temp_dbl = -dd / hh;
        s[idtem] = 1.0 / sqrt(1.0 + temp_dbl * temp_dbl);
        c[idtem] = temp_dbl * s[idtem];
    }
    else
    {
        temp_dbl = -hh / dd;
        c[idtem] = 1.0 / sqrt(1.0 + temp_dbl * temp_dbl);
        s[idtem] = temp_dbl * c[idtem];
    }

    h[idtem]  = c[idtem] * h[idtem] - s[idtem] * h[endtem];
    h[endtem] = 0.0;

    temp_dbl    = c[idtem] * eta[idtem] - s[idtem] * eta[endtem];
    eta[endtem] = s[idtem] * eta[idtem] + c[idtem] * eta[endtem];
    eta[idtem]  = temp_dbl;
}

// Backward calculation
// To notice, Hesssenburg matrix's column
// and row changes due to use Array<OneD,Array<OneD,NekDouble>> --> Put into a
// helper class
void NekLinSysIterGMRESLoc::DoBackward(const int number,
                                       Array<OneD, Array<OneD, NekDouble>> &A,
                                       const Array<OneD, const NekDouble> &b,
                                       Array<OneD, NekDouble> &y)
{
    // Number is the entry number
    // but C++'s order need to be one smaller
    int maxid = number - 1;
    NekDouble sum;
    y[maxid] = b[maxid] / A[maxid][maxid];
    for (int i = maxid - 1; i > -1; --i)
    {
        sum = b[i];
        for (int j = i + 1; j < number; ++j)
        {
            // i and j changes due to use Array<OneD,Array<OneD,NekDouble>>
            sum -= y[j] * A[j][i];
        }
        y[i] = sum / A[i][i];
    }
}
} // namespace Nektar::LibUtilities
