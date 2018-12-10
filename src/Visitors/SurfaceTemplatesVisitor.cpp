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
    
    SurfaceTemplatesVisitor::SurfaceTemplatesVisitor():meshpts(NULL) {
    }

    bool SurfaceTemplatesVisitor::visit(Quadrant *o) {
        
        const vector<unsigned int> &pointindex = o->pointindex;
        if (pointindex.size()!=4 /*|| o->getSubElements().size()>1*/) {
            return true;
        }
        
        //cout << "SurfaceTemplates" << endl;
        //A surface template should be applyed only over elements
        //that intersect one surface or all of them. In both cases
        //at least one element node should be outside the sum of
        //input surfaces
        vector<vector<unsigned int> > newsubelems; // = o->sub_elements;

        //----------------------------------------------------------------------
        // not a transition quadrant
        //----------------------------------------------------------------------
        if (o->getSubElements().size()==1) {

            vector<bool> in(pointindex.size(),false);
            unsigned int nin = 0, onei=0, oneo=0;

            for (unsigned int i=0; i<pointindex.size(); i++) {
                if (!meshpts->at(pointindex[i]).isOutside()){
                    in[i] = true;
                    onei = i;
                    nin++;
                }
                else {
                    oneo = i;
                }
            }

            if (nin==0) {
                //All the nodes of the octant were projected onto the
                //boundary however it is not outside.
                //Every node with angle below 150 will be considered
                //as inside, while the rest as outside. The same logic
                //as the rest of the octant will be employed after this
                //step.
                bool tri = false;


                for (unsigned int i=0; i<pointindex.size(); i++) {
                    double angle = o->getAngle(i,*meshpts);
                    
                    //this is to treat when two next nodes are projected
                    //onto the same place (candidate: become a triangle)
                    if (isnan(angle) && !tri) {
                        tri = true;
                        in[i] = true;
                        nin++;
                        continue;
                    }
                    
                    if (angle<=150.0) {
                        in[i]=true;
                        nin++;
                    }
                    else {
                        //if the angle is not good, this is the candidate (index)
                        //to be used for splitting the Quad into two triangles.
                        //this will be performed by the final switch of this method
                        //in the case of 3 nodes in.
                        oneo=i;
                        if (angle>175. && angle<185.) {
                            tri = true;
                        }
                    }
                }

                //NOTE: nin was updated with badAngles.
                if (nin==4) {
                    //nothing to do, the element will
                    return true;
                }
                if (nin==3 && tri) {
                    //it should be replaced by 1 or 2 triangles.
                    vector<unsigned int> t;
                    t.reserve(3);
                    for (unsigned int i=0; i<4; i++) {
                        if (in[i]) {
                            t.push_back(pointindex[i]);
                        }
                    }
                    newsubelems.push_back(t);
                    o->updateSubElements(newsubelems);
                    return true;
                }
            }
            else {
                if (nin==1) {
                    unsigned int op = (onei+2)%4;
                    bool good_angle = false;
                    //test if opposed node has acceptable angle.
                    if (!o->badAngle(op,*meshpts)) {
                        op = (onei+1)%4;
                        //if we are here we must also test the direct
                        //neighbors of the inside node. They may cause
                        //concave elements due to feature projection.
                        if (!o->badAngle(op,*meshpts)) {
                            op = (onei+3)%4;
                            if (!o->badAngle(op,*meshpts)) {
                                good_angle = true;
                            }
                        }
                    }
                    if (good_angle) {
                        //If the angles are acceptable, test if the center
                        //of the triangle formed by projected nodes (not the
                        //inside one) is outside the domain. In this case,
                        //we must replace the quadrant with a triangle.
                        const Point3D &P0 = meshpts->at(pointindex[(onei+1)%4]).getPoint();
                        const Point3D &P1 = meshpts->at(pointindex[(onei+2)%4]).getPoint();
                        const Point3D &P2 = meshpts->at(pointindex[(onei+3)%4]).getPoint();
                        Point3D centroid = (P0+P1+P2)/3.;

                        if (polyline->pointIsInMesh(centroid,o->getIntersectedEdges())) {
                            return true;
                        }
                    }
                }
            }

            QuadSurfTemplate st;
            bool res;
            switch (nin) {
                case 0:
                    return true;
                case 1:
                    res = st.one(pointindex,in,newsubelems);
                    o->updateSubElements(newsubelems);
                    return res;
                case 2:
                    res = st.two(pointindex,in,newsubelems);
                    o->updateSubElements(newsubelems);
                    return res;
                case 3:
                    //if the angle at the outside node is not close to 180, it is not necessary
                    //to split the Quad into two triangles.
                    if (!o->badAngle(oneo,*meshpts)) {
                        return true;
                    }
                    res = st.three(pointindex,in,newsubelems);
                    o->updateSubElements(newsubelems);
                    return res;
                default:
                    return true;
            }
        } else
        //----------------------------------------------------------------------
        // Transition quadrant, must handle sub-elements....
        //----------------------------------------------------------------------
        {
//            return true;

            const vector <vector <unsigned int> > &subelems=o->getSubElements();
            //First attempt:
            //remove element if all nodes outside (or projected) && gravity outside
            for (unsigned isub=0; isub< subelems.size(); ++isub) {

                const vector <unsigned int> &selem=subelems[isub];

                vector<bool> in(selem.size(),false);
                unsigned int nin = 0, onei=0, oneo=0;
                Point3D centroid;

                for (unsigned int i=0; i<selem.size(); i++) {

                    if (meshpts->at(selem[i]).isInside() /*|| !meshpts->at(selem[i]).wasProjected()*/ ){
                        in[i] = true;
                        onei = i;
                        nin++;
                    }
                    else {
                        oneo = i;
                    }
                    centroid += meshpts->at(selem[i]).getPoint();
                }
                centroid /=selem.size();

                if (nin==0) {
                    if (polyline->pointIsInMesh(centroid,o->getIntersectedEdges())) {

                        //keep the element
                        newsubelems.push_back(subelems[isub]);
                        //return true;
                        continue;
                    }
                } else {
                    newsubelems.push_back(subelems[isub]);
                }
            }
            o->updateSubElements(newsubelems);
            return true;
        }
        //----------------------------------------------------------------------

        return false;
        
    }
}
