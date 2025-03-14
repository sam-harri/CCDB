///////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessRevolve.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description: Revolve a two-dimensional mesh around the y-axis to form a
//  three-dimensional mesh.  Arguments are "angle" and "layers"
//
////////////////////////////////////////////////////////////////////////////////

#include "ProcessRevolve.h"
#include <NekMesh/MeshElements/Element.h>

using namespace std;
using namespace Nektar::NekMesh;

namespace Nektar::NekMesh
{
ModuleKey ProcessRevolve::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "revolve"), ProcessRevolve::create);

ProcessRevolve::ProcessRevolve(MeshSharedPtr m) : ProcessModule(m)
{
    m_config["layers"] = ConfigOption(false, "24", "Number of layers");
    m_config["angle"] =
        ConfigOption(false, "full", "Angle of revolution (in radians)");
}

ProcessRevolve::~ProcessRevolve()
{
}

void ProcessRevolve::Process()
{
    m_log(VERBOSE) << "Revolving grid." << endl;

    if (m_mesh->m_spaceDim != 2)
    {
        m_log(FATAL) << "Revolve should only be called for a two dimensional "
                     << "mesh" << endl;
    }

    int nLayers = m_config["layers"].as<int>();
    // Special case needed for full revolutions
    bool full       = m_config["angle"].as<std::string>() == "full";
    NekDouble angle = full ? 2 * M_PI : m_config["angle"].as<NekDouble>();
    if (angle < 0 || angle > 2 * M_PI)
    {
        m_log(FATAL) << "Angle should be between 0 and 2pi" << endl;
    }
    NekDouble dphi = angle / nLayers;
    // Increment space and expansion dimensions.
    m_mesh->m_spaceDim++;
    m_mesh->m_expDim++;

    // Grab a copy of the existing two-dimensional elements.
    vector<ElementSharedPtr> el = m_mesh->m_element[2];

    // Grab a copy of existing composites.
    CompositeMap oldComp = m_mesh->m_composite;
    m_log(VERBOSE) << "Boundary composites" << endl;
    for (auto &it : oldComp)
    {
        if (it.second->m_tag != "E")
        {
            continue;
        }
        m_log(VERBOSE) << it.first << "\t" << it.second->m_tag;
        for (int i = 0; i < it.second->m_items.size(); ++i)
        {
            m_log(VERBOSE) << "\t" << it.second->m_items[i]->GetId() << " ("
                           << it.second->m_items[i]->GetVertex(0) << ", "
                           << it.second->m_items[i]->GetVertex(1) << ")";
            vector<NodeSharedPtr> vv = it.second->m_items[i]->GetVertexList();
            m_log(VERBOSE) << "\t(" << vv[0]->GetID() << ", " << vv[1]->GetID()
                           << ")";
        }
        m_log(VERBOSE) << endl;
    }

    // Reset mesh.
    for (int d = 0; d <= 3; ++d)
    {
        m_mesh->m_element[d].clear();
    }

    NodeSet nodes = m_mesh->m_vertexSet;

    map<int, NodeSharedPtr> id2node;

    for (auto &n : nodes)
    {
        id2node[n->m_id] = n;
    }
    // Save z plane coordinate
    NekDouble phi0 = 0;

    // Create vertices for subsequent layers.
    for (int i = 1; i < nLayers + 1; ++i)
    {
        if (full && i == nLayers)
        {
            // For last layer we will connect up to first layer
            break;
        }
        for (auto &n : nodes)
        {
            if (n->m_x <= 0)
            {
                // Disallow x=0 for now
                m_log(FATAL)
                    << "All vertices require positive (non-zero) x coordinates"
                    << endl;
            }
            NekDouble r0 = sqrt(n->m_x * n->m_x + n->m_z * n->m_z);
            NodeSharedPtr newNode(new Node(i * nodes.size() + n->m_id,
                                           r0 * cos(i * dphi), n->m_y,
                                           r0 * sin(i * dphi)));

            m_mesh->m_vertexSet.insert(newNode);
            id2node[i * nodes.size() + n->m_id] = newNode;
        }
    }

    EdgeSet esOld = m_mesh->m_edgeSet; // copy edges for curvature

    for (int j = 0; j < nLayers; ++j)
    {
        int next = full ? (j + 1) % nLayers : j + 1;
        for (int i = 0; i < el.size(); ++i)
        {
            vector<NodeSharedPtr> verts = el[i]->GetVertexList();
            if (verts.size() == 4)
            {
                vector<NodeSharedPtr> nodeList(8);

                nodeList[0] = id2node[verts[0]->m_id + j * nodes.size()];
                nodeList[1] = id2node[verts[1]->m_id + j * nodes.size()];
                nodeList[2] = id2node[verts[2]->m_id + j * nodes.size()];
                nodeList[3] = id2node[verts[3]->m_id + j * nodes.size()];
                nodeList[4] = id2node[verts[0]->m_id + next * nodes.size()];
                nodeList[5] = id2node[verts[1]->m_id + next * nodes.size()];
                nodeList[6] = id2node[verts[2]->m_id + next * nodes.size()];
                nodeList[7] = id2node[verts[3]->m_id + next * nodes.size()];

                vector<int> tags(1);
                tags[0] = 0;

                ElmtConfig conf(LibUtilities::eHexahedron, 1, false, false,
                                false);
                ElementSharedPtr E = GetElementFactory().CreateInstance(
                    LibUtilities::eHexahedron, conf, nodeList, tags);

                m_mesh->m_element[3].push_back(E);
            }
            else
            {
                vector<NodeSharedPtr> nodeList(6);
                nodeList[0] = id2node[verts[0]->m_id + next * nodes.size()];
                nodeList[1] = id2node[verts[1]->m_id + next * nodes.size()];
                nodeList[2] = id2node[verts[1]->m_id + j * nodes.size()];
                nodeList[3] = id2node[verts[0]->m_id + j * nodes.size()];
                nodeList[4] = id2node[verts[2]->m_id + next * nodes.size()];
                nodeList[5] = id2node[verts[2]->m_id + j * nodes.size()];

                vector<int> tags(1);
                tags[0] = 1;

                ElmtConfig conf(LibUtilities::ePrism, 1, false, false, false);
                ElementSharedPtr E = GetElementFactory().CreateInstance(
                    LibUtilities::ePrism, conf, nodeList, tags);

                m_mesh->m_element[3].push_back(E);
            }
        }
    }

    ProcessEdges();
    ProcessFaces();
    ProcessElements();
    ProcessComposites();

    // Copy edge information
    for (auto &edge : esOld)
    {
        if (edge->m_edgeNodes.size() > 0)
        {
            for (int j = 0; j < nLayers + 1; ++j)
            {
                if (full && j == nLayers)
                {
                    // For last layer we will connect up to first layer
                    break;
                }
                vector<NodeSharedPtr> ns(edge->m_edgeNodes.size());
                for (int i = 0; i < ns.size(); i++)
                {
                    NodeSharedPtr n = edge->m_edgeNodes[i];
                    NekDouble r0    = sqrt(n->m_x * n->m_x + n->m_z * n->m_z);

                    ns[i] = std::shared_ptr<Node>(new Node(
                        0, r0 * cos(j * dphi), n->m_y, r0 * sin(j * dphi)));
                }

                EdgeSharedPtr e = std::shared_ptr<Edge>(
                    new Edge(id2node[edge->m_n1->m_id + j * nodes.size()],
                             id2node[edge->m_n2->m_id + j * nodes.size()]));

                auto f = m_mesh->m_edgeSet.find(e);
                ASSERTL1(f != m_mesh->m_edgeSet.end(), "could not find edge");

                // Copy edge type
                (*f)->m_curveType = edge->m_curveType;
                // Copy points
                if ((*f)->m_n1 == e->m_n1)
                {
                    (*f)->m_edgeNodes = ns;
                }
                else
                {
                    reverse(ns.begin(), ns.end());
                    (*f)->m_edgeNodes = ns;
                }
            }
        }
    }

    // Get composites max id
    unsigned int maxCompId = 0;
    for (auto &it : oldComp)
    {
        if (it.second->m_id >= maxCompId)
        {
            maxCompId = it.second->m_id;
        }
    }

    // First rename surface to volume composites to out of range
    int outCompId = maxCompId + 1;
    std::vector<int> toErase;
    for (auto &it2 : m_mesh->m_composite)
    {
        if (it2.second->m_id > maxCompId)
        {
            // done!
            break;
        }
        if (it2.second->m_tag == "H" || it2.second->m_tag == "R")
        {
            it2.second->m_id = outCompId;
            m_mesh->m_composite.insert(std::make_pair(outCompId, it2.second));
            toErase.push_back(it2.first);
            outCompId += 1;
        }
    }

    for (auto &e : toErase)
    {
        m_mesh->m_composite.erase(e);
    }

    toErase.clear();

    // Then copy surface to volume composites names
    for (auto &it2 : m_mesh->m_composite)
    {
        if (it2.second->m_tag == "H" || it2.second->m_tag == "R")
        {
            for (auto &it1 : oldComp)
            {
                if (it2.second->m_tag == "H" && it1.second->m_tag == "Q")
                {
                    it2.second->m_id = it1.second->m_id;
                    m_mesh->m_composite.insert(
                        std::make_pair(it1.second->m_id, it2.second));
                    toErase.push_back(it2.first);
                    oldComp.erase(it1.first);
                    break;
                }
                else if (it2.second->m_tag == "R" && it1.second->m_tag == "T")
                {
                    it2.second->m_id = it1.second->m_id;
                    m_mesh->m_composite.insert(
                        std::make_pair(it1.second->m_id, it2.second));
                    toErase.push_back(it2.first);
                    oldComp.erase(it1.first);
                    break;
                }
            }
        }
    }

    for (auto &e : toErase)
    {
        m_mesh->m_composite.erase(e);
    }

    // Add new composite to be filled with all boundary faces
    CompositeSharedPtr comp(new Composite());
    comp->m_id = ++maxCompId;
    unsigned int compAllFaceId =
        maxCompId; // save it so we can remove it later on
    comp->m_tag = "F";
    m_mesh->m_composite.insert(std::make_pair(maxCompId, comp));

    // Add all boundary faces to the composite
    auto allFaceC = m_mesh->m_composite.find(maxCompId);
    m_log(VERBOSE) << "Faces boundary list" << endl;

    for (auto &it : m_mesh->m_faceSet)
    {
        // Add to composite if boundary face
        if (it->m_elLink.size() < 2)
        {
            if (it->m_vertexList.size() == 3)
            {
                // Triangle
                ElmtConfig conf(LibUtilities::eTriangle, 1, false, false);
                vector<int> tags(1);
                tags[0]            = 1;
                ElementSharedPtr E = GetElementFactory().CreateInstance(
                    LibUtilities::eTriangle, conf, it->m_vertexList, tags);
                E->SetId(it->m_id);
                allFaceC->second->m_items.push_back(E);
            }
            else if (it->m_vertexList.size() == 4)
            {
                // Quad
                ElmtConfig conf(LibUtilities::eQuadrilateral, 1, false, false);
                vector<int> tags(1);
                tags[0]            = 0;
                ElementSharedPtr E = GetElementFactory().CreateInstance(
                    LibUtilities::eQuadrilateral, conf, it->m_vertexList, tags);
                E->SetId(it->m_id);
                allFaceC->second->m_items.push_back(E);
            }
        }
    }

    // Create boundary composites
    for (auto &itOc : oldComp)
    {
        CompositeSharedPtr comp(new Composite());
        comp->m_id  = itOc.second->m_id;
        comp->m_tag = "F";
        m_mesh->m_composite.insert(std::make_pair(itOc.second->m_id, comp));
    }
    // Create periodic composites
    if (!full)
    {
        for (int i = 0; i < 2; i++)
        {
            CompositeSharedPtr comp(new Composite());
            comp->m_id  = ++maxCompId;
            comp->m_tag = "F";
            m_mesh->m_composite.insert(std::make_pair(maxCompId, comp));
        }
    }

    // Populates boundary composites
    for (auto &itQ : allFaceC->second->m_items)
    {
        // Check if this quad belongs to previous boundary
        for (auto &itOc : oldComp)
        {
            for (int iEd = 0; iEd < itOc.second->m_items.size(); ++iEd)
            {
                int inCommon = 0;
                for (int iV = 0; iV < itQ->GetVertexList().size(); iV++)
                {
                    for (int j = 0; j < 2; j++)
                    {
                        NekDouble oldr2 =
                            itOc.second->m_items[iEd]->GetVertex(j)->m_x *
                                itOc.second->m_items[iEd]->GetVertex(j)->m_x +
                            itOc.second->m_items[iEd]->GetVertex(j)->m_z *
                                itOc.second->m_items[iEd]->GetVertex(j)->m_z;
                        NekDouble newr2 =
                            itQ->GetVertex(iV)->m_x * itQ->GetVertex(iV)->m_x +
                            itQ->GetVertex(iV)->m_z * itQ->GetVertex(iV)->m_z;

                        if (LibUtilities::IsRealEqual(oldr2, newr2) &&
                            LibUtilities::IsRealEqual(
                                itQ->GetVertex(iV)->m_y,
                                itOc.second->m_items[iEd]->GetVertex(j)->m_y))
                        {
                            ++inCommon;
                        }
                    }
                }
                // If the face contains 4 xy pairs in common with 1 edge it
                // must be an extruded edge and it should be added to the
                // corresponding composite
                if (inCommon == 4)
                {
                    auto newC = m_mesh->m_composite.find(itOc.second->m_id);
                    // Quad
                    ElmtConfig conf(LibUtilities::eQuadrilateral, 1, false,
                                    false);
                    vector<int> tags(1);
                    tags[0]            = 0;
                    ElementSharedPtr E = GetElementFactory().CreateInstance(
                        LibUtilities::eQuadrilateral, conf,
                        itQ->GetVertexList(), tags);
                    E->SetId(itQ->GetId());
                    newC->second->m_items.push_back(E);
                }
            }
        }

        // Populates periodic composites
        if (!full)
        {
            NekDouble phidist = 0.0;
            for (int iV = 0; iV < itQ->GetVertexList().size(); iV++)
            {
                phidist +=
                    (atan2(-itQ->GetVertex(iV)->m_z, -itQ->GetVertex(iV)->m_x) +
                     M_PI);
            }
            phidist                = phidist / itQ->GetVertexList().size();
            unsigned int compPerId = 0;
            if (LibUtilities::IsRealEqual(phidist, phi0))
            {
                compPerId = maxCompId - 1;
            }
            else if (LibUtilities::IsRealEqual(phidist - phi0, angle))
            {
                compPerId = maxCompId;
            }
            if (compPerId > 0 && itQ->GetVertexList().size() == 3)
            {
                // Triangle
                auto perC = m_mesh->m_composite.find(compPerId);
                ElmtConfig conf(LibUtilities::eTriangle, 1, false, false);
                vector<int> tags(1);
                tags[0]            = 1;
                ElementSharedPtr E = GetElementFactory().CreateInstance(
                    LibUtilities::eTriangle, conf, itQ->GetVertexList(), tags);
                E->SetId(itQ->GetId());
                perC->second->m_items.push_back(E);
            }
            else if (compPerId > 0 && itQ->GetVertexList().size() == 4)
            {
                // Quad
                auto perC = m_mesh->m_composite.find(compPerId);
                ElmtConfig conf(LibUtilities::eQuadrilateral, 1, false, false);
                vector<int> tags(1);
                tags[0]            = 0;
                ElementSharedPtr E = GetElementFactory().CreateInstance(
                    LibUtilities::eQuadrilateral, conf, itQ->GetVertexList(),
                    tags);
                E->SetId(itQ->GetId());
                perC->second->m_items.push_back(E);
            }
        }
    }
    // Remove all faces composite
    m_mesh->m_composite.erase(compAllFaceId);
}
} // namespace Nektar::NekMesh
