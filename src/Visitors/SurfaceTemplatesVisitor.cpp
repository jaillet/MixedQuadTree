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

#include "SurfaceTemplatesVisitor.h"
#include "../Quadrant.h"

namespace Clobscode
{

    /*Note: newpts will probably needed later when analysing
     quadrants that presents both, surface and transition patterns.
     */
    
    SurfaceTemplatesVisitor::SurfaceTemplatesVisitor() {
        meshpts = NULL;
        //newpts = NULL;
        input = NULL;
        e_idx = NULL;
    }


    void SurfaceTemplatesVisitor::setPoints(vector<MeshPoint> &meshpts) {
        this->meshpts = &meshpts;
    }

    //This will probably needed later
    /*void SurfaceTemplatesVisitor::setNewPoints(list<MeshPoint> &newpts) {
        this->newpts = &newpts;
    }*/

    void SurfaceTemplatesVisitor::setInput(Polyline &input) {
        this->input = &input;
    }

    void SurfaceTemplatesVisitor::setIdx(unsigned int &e_idx) {
        this->e_idx = &e_idx;
    }

    bool SurfaceTemplatesVisitor::visit(Quadrant *o) {
        
        const vector<unsigned int> &pointindex = o->pointindex;
        if (pointindex.size()>4) {
            return true;
        }
        
        //cout << "SurfaceTemplates" << endl;
        //A surface template should be applyed only over elements
        //that intersect one surface or all of them. In both cases
        //at least one element node should be outside the sum of
        //input surfaces
        vector<vector<unsigned int> > &sub_elements = o->sub_elements;
        vector<bool> in(pointindex.size(),true);
        unsigned int nin = 0;

        for (unsigned int i=0; i<pointindex.size(); i++) {

            if (meshpts->at(pointindex[i]).isOutside()) {
                in[i] = false;
            }
            else {
                nin++;
            }
        }
        
        QuadSurfTemplate st;
        
        switch (nin) {
            case 0:
                return true;
            case 1:
                return st.one(pointindex,in,sub_elements);
            case 2:
                return st.two(pointindex,in,sub_elements);
            case 3:
                return st.three(pointindex,in,sub_elements);
            default:
                return true;
        }

        /*if (pointindex.size()==4 && sub_elements.size()==1) {

            //if we are here, the "inside" projected nodes will count as
            //outside nodes.
            return applyHexSurfaceTemplates(o,in);
        }*/

        return false;

    }
}
