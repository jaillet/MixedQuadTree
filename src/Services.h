/*
 <Mix-mesher: region type. This program generates a mixed-elements 2D mesh>

 Copyright (C) <2013,2020>  <Claudio Lobos> All rights reserved.

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/lgpl.txt>
 */
/**
* @file Services.h
* @author Claudio Lobos, Fabrice Jaillet
* @version 0.1
* @brief
**/

#ifndef Services_h
#define Services_h 1

#include "Polyline.h"
#include "FEMesh.h"
#include "RefinementCubeRegion.h"
#include "RefinementSurfaceRegion.h"
#include "RefinementInputSurfaceRegion.h"
#include "RefinementAllRegion.h"
#include "MeshPoint.h"
#include "Quadrant.h"
#include "QuadEdge.h"
#include <stdlib.h>
#include <memory>

using Clobscode::Point3D;
using Clobscode::Polyline;
using Clobscode::MeshPoint;
using Clobscode::Quadrant;
using Clobscode::QuadEdge;
using std::vector;

namespace Clobscode
{

class Services {

public:


    //-------------------------------------------------------------------
    //-------------------------------------------------------------------

    static unsigned short readRefinementRegions(string name,
                                                list<RefinementRegion *> &regions);

    //-------------------------------------------------------------------
    //-------------------------------------------------------------------
    static bool readSurfaceRefinementRegion(string name,
                                            list<RefinementRegion *> &regions,
                                            const unsigned short &rrl);

    //-------------------------------------------------------------------
    //-------------------------------------------------------------------
    static bool ReadOffMesh(string name,
                            vector<Clobscode::Polyline> &clobs_inputs);


    //-------------------------------------------------------------------
    //-------------------------------------------------------------------
    static bool ReadMdlMesh(std::string name,
                            vector<Clobscode::Polyline> &clobs_inputs);

    //-------------------------------------------------------------------
    //-------------------------------------------------------------------
    //        First line: <# of vertices> <dimension (must be 2)> <# of attributes> <# of boundary markers (0 or 1)>
    //        Following lines: <vertex #> <x> <y> [attributes] [boundary marker]
    //        One line: <# of segments> <# of boundary markers (0 or 1)>
    //        Following lines: <segment #> <endpoint> <endpoint> [boundary marker]
    //        One line: <# of holes>
    //        Following lines: <hole #> <x> <y>
    //        Optional line: <# of regional attributes and/or area constraints>
    //        Optional following lines: <region #> <x> <y> <attribute> <maximum area>

    static bool ReadPolyFile(std::string name,
                             vector<Clobscode::Polyline> &clobs_inputs);


    //-------------------------------------------------------------------
    //-------------------------------------------------------------------
    static bool ReadQuadrantList(std::string name, list<unsigned int> &olist,
                                 vector<unsigned int> &ele_quadt_ref) ;

    //-------------------------------------------------------------------
    //-------------------------------------------------------------------
    static bool ReadQuadMesh(std::string name, vector<MeshPoint> &points,
                             vector<Quadrant> &Quadrants,
                             set<QuadEdge> &edges,
                             vector<unsigned int> &ele_oct_ref,
                             GeometricTransform &gt,
                             unsigned short &minrl,
                             unsigned short &maxrl);

    //-------------------------------------------------------------------
    //-------------------------------------------------------------------
    static bool WriteQuadtreeMesh(std::string name, const vector<MeshPoint> &points,
                                  const vector<Quadrant> &Quadrants,
                                  const set<QuadEdge> &edges,
                                  unsigned int nels,
                                  const GeometricTransform &gt);

    //-------------------------------------------------------------------
    //-------------------------------------------------------------------
    static bool WriteGMSH(std::string name, const shared_ptr<FEMesh> &output);

    //-------------------------------------------------------------------
    //-------------------------------------------------------------------
    static bool WriteVTK(std::string name, const shared_ptr<FEMesh> &output);

    //-------------------------------------------------------------------
    //-------------------------------------------------------------------
    static bool WritePolyVTK(std::string name, vector<Polyline> inputs);

    //-------------------------------------------------------------------
    //-------------------------------------------------------------------
    static bool WritePolyFile(std::string name, vector<Polyline> inputs);

    //-------------------------------------------------------------------
    //-------------------------------------------------------------------
    static bool WriteOFF(std::string name, const shared_ptr<FEMesh> &output);

    //-------------------------------------------------------------------
    //-------------------------------------------------------------------
    static bool WriteMixedVolumeMesh(std::string name, const shared_ptr<FEMesh> &output);

    //-------------------------------------------------------------------
    //-------------------------------------------------------------------
    static bool WriteOutputMesh(std::string name, const shared_ptr<FEMesh> &output);

    //-------------------------------------------------------------------
    //-------------------------------------------------------------------
    static bool WriteMeshGetfem(std::string name, const std::shared_ptr<FEMesh> &output);

    //-------------------------------------------------------------------
    //-------------------------------------------------------------------
    static bool WriteHistogram(std::string name, const std::shared_ptr<FEMesh> &output);
};

}
#endif
