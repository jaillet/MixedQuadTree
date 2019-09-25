/*
 <Mix-mesher: region type. This program generates a mixed-elements mesh>
 
 Copyright (C) <2013,2017>  <Claudio Lobos>
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/gpl.txt>
 */

#ifndef SurfaceTemplatesVisitor_h
#define SurfaceTemplatesVisitor_h 1

#include "../QuadEdge.h"
#include "../MeshPoint.h"
#include "../Polyline.h"
#include "../Quadrant.h"
#include "../QuadSurfTemplate.h"

#include "Visitor.h"

#include <set>
#include <vector>

using Clobscode::MeshPoint;
using Clobscode::QuadEdge;
using std::set;
using std::vector;

namespace Clobscode
{
    class SurfaceTemplatesVisitor : public Visitor {
    public:
        SurfaceTemplatesVisitor();

        bool visit(Quadrant *o) override;

        void setPoints(const vector<MeshPoint> &meshpts);
        
        void setPolyline(const Polyline &polyline);

    private:
        
        const vector<MeshPoint> *meshpts;
        const Polyline *polyline;

    };
    
    inline void SurfaceTemplatesVisitor::setPoints(const vector<MeshPoint> &meshpts) {
        this->meshpts = &meshpts;
    }
    
    inline void SurfaceTemplatesVisitor::setPolyline(const Polyline &polyline) {
        this->polyline = &polyline;
    }
}

#endif //MESHER_ROI_SURFACETEMPLATESVISITOR_H
