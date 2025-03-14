///////////////////////////////////////////////////////////////////////////////
//
// File: MeshPartitionPtScotch.cpp
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
// Description: PtScotch partitioner interface
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Communication/CommMpi.h>
#include <SpatialDomains/MeshPartitionPtScotch.h>

#include <ptscotch.h>

#define SCOTCH_CALL(scotchFunc, args)                                          \
    {                                                                          \
        ASSERTL0(scotchFunc args == 0,                                         \
                 std::string("Error in Scotch calling function ") +            \
                     std::string(#scotchFunc));                                \
    }

namespace Nektar::SpatialDomains
{

std::string MeshPartitionPtScotch::className =
    GetMeshPartitionFactory().RegisterCreatorFunction(
        "PtScotch", MeshPartitionPtScotch::create,
        "Parallel partitioning using the PtScotch library.");

std::string MeshPartitionPtScotch::cmdSwitch =
    LibUtilities::SessionReader::RegisterCmdLineFlag(
        "use-ptscotch", "", "Use PtScotch for parallel mesh partitioning.");

MeshPartitionPtScotch::MeshPartitionPtScotch(
    const LibUtilities::SessionReaderSharedPtr session,
    LibUtilities::CommSharedPtr comm, int meshDim,
    std::map<int, MeshEntity> element, CompositeDescriptor compMap)
    : MeshPartition(session, comm, meshDim, element, compMap)
{
    m_parallel = true;
}

MeshPartitionPtScotch::~MeshPartitionPtScotch()
{
}

void MeshPartitionPtScotch::v_PartitionGraphImpl(
    int &nVerts, [[maybe_unused]] int &nVertConds, Array<OneD, int> &xadj,
    Array<OneD, int> &adjcy, Array<OneD, int> &vertWgt,
    [[maybe_unused]] Array<OneD, int> &vertSize,
    [[maybe_unused]] Array<OneD, int> &edgeWgt, int &nparts,
    [[maybe_unused]] int &volume, Array<OneD, int> &part)
{
    LibUtilities::CommMpiSharedPtr mpiComm =
        std::dynamic_pointer_cast<LibUtilities::CommMpi>(
            m_comm->GetSpaceComm());

    ASSERTL0(mpiComm, "PtScotch not supported in serial execution.");

    SCOTCH_Dgraph scGraph;
    SCOTCH_CALL(SCOTCH_dgraphInit, (&scGraph, mpiComm->GetComm()));
    SCOTCH_CALL(SCOTCH_dgraphBuild,
                (&scGraph, 0, nVerts, nVerts, &xadj[0], &xadj[1], &vertWgt[0],
                 nullptr, adjcy.size(), adjcy.size(), &adjcy[0], nullptr,
                 nullptr));
    SCOTCH_CALL(SCOTCH_dgraphCheck, (&scGraph));

    SCOTCH_Strat strat;
    SCOTCH_CALL(SCOTCH_stratInit, (&strat));
    SCOTCH_CALL(SCOTCH_stratDgraphMapBuild,
                (&strat, SCOTCH_STRATQUALITY, nparts, nparts, 0.05));

    SCOTCH_CALL(SCOTCH_dgraphPart, (&scGraph, nparts, &strat, &part[0]));
}

} // namespace Nektar::SpatialDomains
