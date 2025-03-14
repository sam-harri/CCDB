////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessPerAlign.cpp
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
//  Description: Reorder composites to align periodic boundaries.
//
////////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/CheckedCast.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include <LocalRegions/NodalTriExp.h>
#include <LocalRegions/QuadExp.h>
#include <LocalRegions/SegExp.h>
#include <LocalRegions/TriExp.h>

#include <NekMesh/MeshElements/Element.h>

#include "ProcessPerAlign.h"

#include <LibUtilities/Interpreter/Interpreter.h>
#include <boost/algorithm/string.hpp>

#include <LibUtilities/BasicUtils/ParseUtils.h>

using namespace std;
using namespace Nektar::NekMesh;

namespace Nektar::NekMesh
{

ModuleKey ProcessPerAlign::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "peralign"), ProcessPerAlign::create,
        "Align periodic boundaries");

/**
 * @class ProcessPerAlign
 */

/**
 * @brief Default constructor.
 */
ProcessPerAlign::ProcessPerAlign(MeshSharedPtr m) : ProcessModule(m)
{
    m_config["surf1"] =
        ConfigOption(false, "-1", "Tag identifying first surface.");
    m_config["surf2"] =
        ConfigOption(false, "-1", "Tag identifying first surface.");
    m_config["dir"] = ConfigOption(
        false, "",
        "Direction in which to align (either x, y, or z, "
        "or vector with components separated by a comma). "
        "If rot is specified this is interpreted as the axis of rotation");
    m_config["rot"] = ConfigOption(
        false, "", "Rotation to align composites in radians, e.g. PI/20");
    m_config["orient"] =
        ConfigOption(true, "0", "Attempt to reorient tets and prisms");
    m_config["tolfac"] = ConfigOption(
        false, "4",
        "Tolerance factor to which to check planes "
        "are the same after rotation/translation (default tolfac=4)");
    m_config["abstol"] = ConfigOption(
        false, "0",
        "Absolute tolerance to check if planes "
        "are the same after rotation/translation (default abstol=0)");
}

/**
 * @brief Destructor.
 */
ProcessPerAlign::~ProcessPerAlign()
{
}

void ProcessPerAlign::Process()
{
    m_log(VERBOSE)
        << "This module provides certain capability allowing "
           "orderings elements and vertices for multiple periodicities. ! "
        << endl;
    m_log(VERBOSE) << "Limitations : " << endl;
    m_log(VERBOSE) << "- Hexahedral meshes - none (run without 'orient' on the "
                      "final mesh)."
                   << endl;
    m_log(VERBOSE)
        << "- Tetrahedral meshes - none (run on the linear mesh with 'orient') "
           "use NekMesh for HO mesh."
        << endl;
    m_log(VERBOSE) << "- Hybrid Tet-Prism meshes - it retains the curvature of "
                      "the Prisms(only edgeNodes), but not on the Tet (run "
                      "after the BL splitting with 'orient') "
                   << endl;

    string surf1IO       = m_config["surf1"].as<string>();
    string surf2IO       = m_config["surf2"].as<string>();
    string dirIO         = m_config["dir"].as<string>();
    string rotIO         = m_config["rot"].as<string>();
    bool orient          = m_config["orient"].as<bool>();
    string toleranceFact = m_config["tolfac"].as<string>();
    string absoluteTol   = m_config["abstol"].as<string>();

    // create a periodicity couples
    vector<unsigned int> surf1;
    if (surf1IO.size() > 0)
    {
        ParseUtils::GenerateSeqVector(surf1IO, surf1);
    }
    else
    {
        m_log(FATAL) << "Define surf1." << endl;
    }

    vector<unsigned int> surf2;
    if (surf2IO.size() > 0)
    {
        ParseUtils::GenerateSeqVector(surf2IO, surf2);
    }
    else
    {
        m_log(FATAL) << "Define surf2." << endl;
    }

    vector<string> dirV;
    if (dirIO.size())
    {
        boost::split(dirV, dirIO, boost::is_any_of(";"));
    }

    vector<string> rotV;
    if (rotIO.size())
    {
        boost::split(rotV, rotIO, boost::is_any_of(";"));
    }

    if (surf1.size() != surf2.size())
    {
        m_log(FATAL)
            << " There is an unequal number of surface couples surf1, surf2."
            << endl;
    }

    int cnt = 0;
    map<int, pair<FaceSharedPtr, vector<int>>> perFaces;
    for (int pIt = 0; pIt < surf1.size(); pIt++)
    {
        string dir = "false";
        if (dirV.size() == surf1.size())
        {
            dir = dirV[pIt];
        }
        else
        {
            m_log(VERBOSE)
                << "Directions not selected or inconsistent with surface pairs"
                << endl;
        }

        string rot = "false";
        if (rotV.size() == surf1.size())
        {
            rot = rotV[pIt];
        }
        else
        {
            m_log(VERBOSE)
                << "Rotations not selected or inconsistent with surface pairs"
                << endl;
        }

        m_log(VERBOSE) << "surf1= " << surf1[pIt] << endl;
        m_log(VERBOSE) << "surf2= " << surf2[pIt] << endl;

        vector<string> tmp1;
        // boost::split() will return a vector of size 1 if input is an empty
        // string and a delimiter is specified
        if (dir.size() && dir != "false")
        {
            m_log(VERBOSE) << "dir = " << dir << endl;
            boost::split(tmp1, dir, boost::is_any_of(","));
        }

        vector<string> tmp2;
        if (rot.size() && rot != "false")
        {
            m_log(VERBOSE) << "rot = " << dir << endl;
            boost::split(tmp2, rot, boost::is_any_of(","));
        }
        bool rotalign = false;

        NekDouble alignDir[3] = {0.0, 0.0, 0.0};
        NekDouble rotangle    = 0.0;

        if (rot.size() && rot != "false")
        {
            // Rotationally periodic boundary
            rotalign = true;

            // Evaluate expression since may be give as function of PI
            LibUtilities::Interpreter strEval;
            int ExprId = strEval.DefineFunction(" ", tmp2[0]);
            rotangle   = strEval.Evaluate(ExprId);

            // Negate angle since we want to rotate second composite back
            // to this one.
            rotangle *= -1.0;
            // Check that user specified axis of rotation
            if (tmp1.size() != 1)
            {
                m_log(WARNING) << "dir must be specified to x, y, or z when "
                                  "\"rot\" is used. Skipping periodic alignment"
                               << endl;
                return;
            }
        }
        if (!dir.size() && m_mesh->m_spaceDim == 2 && m_mesh->m_cad)
        {
            // if the direction is not specified and its a 2D mesh and there is
            // CAD it can figure out the dir on its own
            auto T = m_mesh->m_cad->GetPeriodicTranslationVector(surf1[pIt],
                                                                 surf2[pIt]);
            NekDouble mag = sqrt(T[0] * T[0] + T[1] * T[1]);

            alignDir[0] = T[0] / mag;
            alignDir[1] = T[1] / mag;
            alignDir[2] = T[2] / mag;
        }
        else if (tmp1.size() == 1)
        {
            if (dir != "x" && dir != "y" && dir != "z")
            {
                m_log(WARNING) << "dir must be set to either x, y or z. "
                               << "Skipping periodic alignment." << endl;
                return;
            }

            alignDir[0] = (dir == "x") ? 1.0 : 0.0;
            alignDir[1] = (dir == "y") ? 1.0 : 0.0;
            alignDir[2] = (dir == "z") ? 1.0 : 0.0;
        }
        else if (tmp1.size() == 3)
        {
            alignDir[0] = std::stod(tmp1[0]);
            alignDir[1] = std::stod(tmp1[1]);
            alignDir[2] = std::stod(tmp1[2]);
        }
        else
        {
            m_log(FATAL) << "Expected three domponents or letter for option "
                            "\"dir\". Skipping periodic alignment,"
                         << endl;
            return;
        }

        auto it1 = m_mesh->m_composite.find(surf1[pIt]);
        auto it2 = m_mesh->m_composite.find(surf2[pIt]);

        if (it1 == m_mesh->m_composite.end())
        {
            m_log(WARNING) << "Couldn't find surface " << surf1[pIt]
                           << ". Skipping periodic alignment." << endl;
            return;
        }

        if (it2 == m_mesh->m_composite.end())
        {
            m_log(WARNING) << "Couldn't find surface " << surf2[pIt] << ", "
                           << "skipping periodic alignment." << endl;
            return;
        }

        CompositeSharedPtr c1 = it1->second;
        CompositeSharedPtr c2 = it2->second;

        if (c1->m_items.size() != c2->m_items.size())
        {
            m_log(WARNING) << "Surfaces " << surf1[pIt] << " and " << surf2[pIt]
                           << " have different numbers of elements. Skipping"
                           << " periodic alignment." << endl;
            return;
        }

        c1->m_reorder = false;
        c2->m_reorder = false;

        // Loop over elements, calculate centroids of elements in c2.
        map<int, Node> centroidMap;
        for (size_t i = 0; i < c2->m_items.size(); ++i)
        {
            Node centroid;
            for (size_t j = 0; j < c2->m_items[i]->GetVertexCount(); ++j)
            {
                centroid += *(c2->m_items[i]->GetVertex(j));
            }
            centroid /= (NekDouble)c2->m_items[i]->GetVertexCount();

            if (rotalign) // rotate centroid
            {
                centroid.Rotate(dir, rotangle);
            }

            centroidMap[i] = centroid;
        }

        std::unordered_set<int> elmtDone;
        map<int, int> elmtPairs;
        map<int, int> vertCheck;

        for (size_t i = 0; i < c1->m_items.size(); ++i)
        {
            Node centroid;
            for (size_t j = 0; j < c1->m_items[i]->GetVertexCount(); ++j)
            {
                centroid += *(c1->m_items[i]->GetVertex(j));
            }
            centroid /= (NekDouble)c1->m_items[i]->GetVertexCount();

            bool found           = false;
            unsigned int tolFact = LibUtilities::checked_cast<unsigned int>(
                std::stod(toleranceFact));
            NekDouble absTol = std::stod(absoluteTol);
            for (auto &it : centroidMap)
            {
                if (elmtDone.count(it.first) > 0)
                {
                    continue;
                }

                Node dx = it.second - centroid;
                bool match;
                if (rotalign)
                {
                    // match = it.second == centroid;
                    match = IsNodeEqual(it.second, centroid, tolFact) ||
                            IsNodeClose(it.second, centroid, absTol);
                }
                else
                {
                    // Check normalized inner product
                    NekDouble normInnProd =
                        fabs(dx.m_x * alignDir[0] + dx.m_y * alignDir[1] +
                             dx.m_z * alignDir[2]) /
                        sqrt(dx.abs2());
                    match =
                        LibUtilities::IsRealEqual(normInnProd, 1.0, tolFact) ||
                        LibUtilities::IsRealClose(normInnProd, 1.0, absTol);
                }

                if (match)
                {
                    // Found match
                    int id1, id2;

                    if (c1->m_items[i]->GetConf().m_e == LibUtilities::eSegment)
                    {
                        id1 = c1->m_items[i]->GetEdgeLink()->m_id;
                        id2 = c2->m_items[it.first]->GetEdgeLink()->m_id;
                    }
                    else
                    {
                        id1 = c1->m_items[i]->GetFaceLink()->m_id;
                        id2 = c2->m_items[it.first]->GetFaceLink()->m_id;
                    }

                    elmtDone.insert(it.first);
                    elmtPairs[i] = it.first;

                    // Identify periodic vertices
                    size_t nVerts = c1->m_items[i]->GetVertexCount();
                    vector<int> perVerts(nVerts, 0), perVertsInv(nVerts, 0);

                    if (orient)
                    {
                        for (size_t k = 0; k < nVerts; ++k)
                        {
                            NodeSharedPtr n1 =
                                c1->m_items[i]->GetFaceLink()->m_vertexList[k];
                            size_t l;
                            NekDouble mindn = 1000;

                            for (l = 0; l < nVerts; ++l)
                            {
                                NodeSharedPtr n2 = c2->m_items[it.first]
                                                       ->GetFaceLink()
                                                       ->m_vertexList[l];

                                if (rotalign) // rotate n2
                                {
                                    Node n2tmp = *n2;
                                    n2tmp.Rotate(dir, rotangle);
                                    // Check if same node
                                    // match = n2tmp == *n1;
                                    match = IsNodeEqual(n2tmp, *n1, tolFact) ||
                                            IsNodeClose(n2tmp, *n1, absTol);
                                    // Compute distance
                                    Node dn         = n2tmp - *n1;
                                    NekDouble dnabs = sqrt(dn.abs2());
                                    mindn = (dnabs < mindn) ? dnabs : mindn;
                                }
                                else
                                {
                                    Node dn = *n2 - *n1;

                                    // Check normalized inner product
                                    NekDouble dnabs =
                                        fabs(dn.m_x * alignDir[0] +
                                             dn.m_y * alignDir[1] +
                                             dn.m_z * alignDir[2]) /
                                        sqrt(dn.abs2());
                                    match = LibUtilities::IsRealEqual(
                                                dnabs, 1.0, tolFact) ||
                                            LibUtilities::IsRealClose(
                                                dnabs, 1.0, absTol);
                                    mindn = (dnabs < mindn) ? dnabs : mindn;
                                }

                                if (match)
                                {
                                    perVerts[k]    = l;
                                    perVertsInv[l] = k;

                                    int id1 = n1->m_id;
                                    int id2 = n2->m_id;
                                    if (vertCheck.count(id1) == 0)
                                    {
                                        vertCheck[id1] = id2;
                                    }
                                    else
                                    {
                                        ASSERTL0(vertCheck[id1] == id2,
                                                 "Periodic vertex already "
                                                 "identified!");
                                    }
                                    break;
                                }
                            }
                            ASSERTL1(l < nVerts,
                                     "Could not identify periodic vertices, "
                                     "nearest distance was " +
                                         boost::lexical_cast<string>(mindn));
                        }

                        int tot1 = 0, tot2 = 0;
                        for (size_t k = 0; k < nVerts; ++k)
                        {
                            tot1 += perVerts[k];
                            tot2 += perVertsInv[k];
                        }
                        ASSERTL0(tot1 == nVerts * (nVerts - 1) / 2 &&
                                     tot2 == nVerts * (nVerts - 1) / 2,
                                 "Error identifying periodic vertices");
                    }

                    if (c2->m_items[i]->GetConf().m_e != LibUtilities::eSegment)
                    {
                        // manually adding the flag for the
                        perVerts.push_back(cnt);
                        perVertsInv.push_back(cnt);
                        // construct perFaces
                        perFaces[id1] = make_pair(
                            c2->m_items[it.first]->GetFaceLink(), perVerts);
                        perFaces[id2] = make_pair(c1->m_items[i]->GetFaceLink(),
                                                  perVertsInv);
                    }
                    found = true;
                    break;
                }
            }

            if (!found)
            {
                m_log(WARNING) << "Could not find matching edge for surface "
                               << "element " << c1->m_items[i]->GetId() << ". "
                               << "Rot = " << rot << "."
                               << "Skipping periodic alignment." << endl;
                return;
            }
        }

        // Reorder vectors.
        vector<ElementSharedPtr> tmp = c2->m_items;

        for (size_t i = 0; i < tmp.size(); ++i)
        {
            c2->m_items[i] = tmp[elmtPairs[i]];
        }

        // Important the pair flag about for the perFace
        cnt++;
    }

    if (orient)
    {
        ReorderPrisms(perFaces);
    }

    m_log(VERBOSE) << " Orient(necessary for Tet/Prisms) = " << orient << endl;
}

} // namespace Nektar::NekMesh
