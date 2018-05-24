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

    void SurfaceTemplatesVisitor::setInput(TriMesh &input) {
        this->input = &input;
    }

    void SurfaceTemplatesVisitor::setIdx(unsigned int &e_idx) {
        this->e_idx = &e_idx;
    }

    bool SurfaceTemplatesVisitor::visit(Quadrant *o) {
        //cout << "SurfaceTemplates" << endl;
        //A surface template should be applyed only over elements
        //that intersect one surface or all of them. In both cases
        //at least one element node should be outside the sum of
        //input surfaces
        list<unsigned int>::iterator piter;
        vector<unsigned int> &pointindex = o->pointindex;
        vector<vector<unsigned int> > &sub_elements = o->sub_elements;
        vector<bool> in(pointindex.size(),true);
        unsigned int nin = 0;

        for (unsigned int i=0; i<pointindex.size(); i++) {

            if (meshpts->at(pointindex[i]).isOutside()) {
                in[i] = false;
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
        
        /*list<vector<unsigned int> >::iterator ne_iter;
        sub_elements.clear();
        sub_elements.reserve(new_eles_lst.size());
        for (ne_iter=new_eles_lst.begin(); ne_iter!=new_eles_lst.end(); ne_iter++) {
            sub_elements.push_back(*ne_iter);
        }*/

        return false;

    }

    /*bool SurfaceTemplatesVisitor::applyHexSurfaceTemplates(Quadrant *o, vector<bool> &inpts) {
        
        vector<unsigned int> &pointindex = o->pointindex;
        vector<vector<unsigned int>> &sub_elements = o->sub_elements;

        //if the Quadrant has more the 8 nodes, it was split by surface or
        //transition patterns. In either case, this method is not able
        //to split it again.
        if (pointindex.size()!=8) {
            return true;
        }

        //The sub-elements of this Quadrant will be now the ones
        //given by the surface pattern that is employed.
        sub_elements.clear();

        //select the patter to apply
        switch (inpts.size()) {
            case 0:{
                //If at this point, the element has 0 node inside,
                 it might be tangencial to input mesh, in which case
                 it should be removed, or represent a feature of
                 the domain (e.g. all nodes outside, but there is
                 something like a pipeline crossing it). This algorithm
                 isn't yet "future sensitive", therefore the element
                 is simply removed.
     
                //sub_elements.push_back(pointindex);
                //return false;
                sub_elements.clear();
                return true;
            }
            case 1: {
                SurfTemplate1 surf_t1;
                return surf_t1.getSubelements(pointindex,inpts,sub_elements);
            }
            case 2: {
                SurfTemplate2 surf_t2;
                return surf_t2.getSubelements(pointindex,inpts,sub_elements);
            }
            case 3: {
                SurfTemplate3 surf_t3;
                return surf_t3.getSubelements(pointindex,inpts,sub_elements);
            }
            case 4: {
                SurfTemplate4 surf_t4;
                return surf_t4.getSubelements(pointindex,inpts,sub_elements);
            }
            case 5: {
                SurfTemplate5 surf_t5;
                unsigned int old_size = newpts->size();
                bool succeed = surf_t5.getSubelements(pointindex,inpts,*meshpts,
                                                      *newpts,sub_elements,*input,
                                                      o->intersected_faces,*e_idx);
                if (old_size != newpts->size()) {
                    pointindex.push_back((meshpts->size()+newpts->size())-1);
                }
                return succeed;
            }
            case 6: {
                SurfTemplate6 surf_t6;
                return surf_t6.getSubelements(pointindex,outpts,sub_elements);
            }
            case 7: {
                SurfTemplate7 surf_t7;
                return surf_t7.getSubelements(pointindex,outpts,sub_elements);
            }
            case 8: {
                //If this happens the element is inside the overall
                //geometry, but intersects inner features.
                sub_elements.push_back(pointindex);
                return true;
            }
            default: {
                cerr << " Error at EnhancedElement::applySurfacePatterns\n";
                cerr << " Number of inside nodes: " << inpts.size() << "\n";
                cerr << " Surface Patterns must be applied over elements";
                cerr << " with 1 to 7 inside nodes.\n";
                return false;
            }
        }
    }*/
}