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

#include "QuadTransition.h"

namespace patterns {
    
    QuadTransition::QuadTransition(vector<unsigned int> &nodes) {
        this->nodes = nodes;
    }
	
    bool QuadTransition::getNewElements(vector<vector<unsigned int> > &sub_elements,
                                        const unsigned int &nedges) {
        switch (nedges) {
            case 1:
                return oneEdge(sub_elements);
            case 2:
                return twoEdges(sub_elements);
            case 3:
                return threeEdges(sub_elements);
            case 4:
                return fourEdges(sub_elements);
            default:
                return false;
        }
        return false;
    }
    
    bool QuadTransition::oneEdge(vector<vector<unsigned int> > &sub_elements) {
    
        unsigned int rotation = 0;
        for (unsigned int i=0; i<4; i++) {
            if (nodes[i+4]!=0) {
                break;
            }
            rotation++;
        }
        if (rotation!=0) {
            nodes = rotate(nodes,rotation);
        }
        
        sub_elements.clear();
        sub_elements.reserve(3);
        
        vector<unsigned int> t1(3,0);
        t1[0]=nodes[0];
        t1[1]=nodes[4];
        t1[2]=nodes[3];
        
        vector<unsigned int> t2(3,0);
        t2[0]=nodes[3];
        t2[1]=nodes[4];
        t2[2]=nodes[2];
        
        vector<unsigned int> t3(3,0);
        t3[0]=nodes[2];
        t3[1]=nodes[4];
        t3[2]=nodes[1];
        
        sub_elements.push_back(t1);
        sub_elements.push_back(t2);
        sub_elements.push_back(t3);
        
        return true;
    }
    
    bool QuadTransition::twoEdges(vector<vector<unsigned int> > &sub_elements) {
        
        unsigned int rotation = 0;
        for (unsigned int i=0; i<4; i++) {
            if (nodes[i+4]!=0) {
                break;
            }
            rotation++;
        }
        
        if (rotation!=0) {
            nodes = rotate(nodes,rotation);
        }
        
        sub_elements.clear();
        
        if (nodes[6]!=0) {
            sub_elements.reserve(2);
            vector<unsigned int> q1(4,0);
            q1[0] = nodes[0];
            q1[1] = nodes[4];
            q1[2] = nodes[6];
            q1[3] = nodes[3];
            
            vector<unsigned int> q2(4,0);
            q2[0] = nodes[4];
            q2[1] = nodes[1];
            q2[2] = nodes[2];
            q2[3] = nodes[6];
            
            sub_elements.push_back(q1);
            sub_elements.push_back(q2);
            
            return true;
        }
        
        if (nodes[7]!=0) {
            nodes = inverseRotate(nodes,1);
        }
        
        sub_elements.reserve(4);
        
        vector<unsigned int> t1(3,0);
        t1[0]=nodes[0];
        t1[1]=nodes[4];
        t1[2]=nodes[3];
        
        vector<unsigned int> t2(3,0);
        t2[0]=nodes[4];
        t2[1]=nodes[5];
        t2[2]=nodes[3];
        
        vector<unsigned int> t3(3,0);
        t3[0]=nodes[5];
        t3[1]=nodes[2];
        t3[2]=nodes[3];
        
        vector<unsigned int> t4(3,0);
        t4[0]=nodes[4];
        t4[1]=nodes[1];
        t4[2]=nodes[5];
        
        sub_elements.push_back(t1);
        sub_elements.push_back(t2);
        sub_elements.push_back(t3);
        sub_elements.push_back(t4);
        
        return true;
    }
    
    bool QuadTransition::threeEdges(vector<vector<unsigned int> > &sub_elements) {
        
        if (nodes[4]==0) {
            nodes = rotate(nodes,2);
        }
        else {
            if (nodes[5]==0) {
                nodes = inverseRotate(nodes,1);
            }
            else {
                if (nodes[7]==0) {
                    nodes = rotate(nodes,1);
                }
            }
        }
        
        sub_elements.clear();
        sub_elements.reserve(4);
        
        vector<unsigned int> t1(3,0);
        t1[0]=nodes[0];
        t1[1]=nodes[4];
        t1[2]=nodes[7];
        
        vector<unsigned int> t2(3,0);
        t2[0]=nodes[4];
        t2[1]=nodes[5];
        t2[2]=nodes[7];
        
        vector<unsigned int> t3(3,0);
        t3[0]=nodes[4];
        t3[1]=nodes[1];
        t3[2]=nodes[5];
        
        vector<unsigned int> q1(4,0);
        q1[0] = nodes[7];
        q1[1] = nodes[5];
        q1[2] = nodes[2];
        q1[3] = nodes[3];
        
        sub_elements.push_back(t1);
        sub_elements.push_back(t2);
        sub_elements.push_back(t3);
        sub_elements.push_back(q1);
        
        return true;
    }
    
    bool QuadTransition::fourEdges(vector<vector<unsigned int> > &sub_elements) {
        
        sub_elements.clear();
        sub_elements.reserve(5);
        
        vector<unsigned int> t1(3,0);
        t1[0]=nodes[0];
        t1[1]=nodes[4];
        t1[2]=nodes[7];
        
        vector<unsigned int> t2(3,0);
        t2[0]=nodes[4];
        t2[1]=nodes[1];
        t2[2]=nodes[5];
        
        vector<unsigned int> t3(3,0);
        t3[0]=nodes[5];
        t3[1]=nodes[2];
        t3[2]=nodes[6];
        
        vector<unsigned int> t4(3,0);
        t3[0]=nodes[6];
        t3[1]=nodes[3];
        t3[2]=nodes[7];
        
        vector<unsigned int> q1(4,0);
        q1[0] = nodes[4];
        q1[1] = nodes[5];
        q1[2] = nodes[6];
        q1[3] = nodes[1];
        
        sub_elements.push_back(t1);
        sub_elements.push_back(t2);
        sub_elements.push_back(t3);
        sub_elements.push_back(t4);
        sub_elements.push_back(q1);
        
        return true;
    }
    
    vector<unsigned int> QuadTransition::rotate(const vector<unsigned int> &nodes,
                                                       const unsigned int &times) {
        
        /*
         0 1 2 3 4 5 6 7
         | | | | | | | | >>1
         1 2 3 0 5 6 7 4
         */
        
        vector<unsigned int> res(8,0);
        for (unsigned int i=0; i<4; i++) {
            unsigned int pos = (i+times)%4;
            res[i] = nodes[pos];
            res[i+4] = nodes[pos+4];
        }
        return res;
    }
    
    vector<unsigned int> QuadTransition::inverseRotate(const vector<unsigned int> &nodes,
                                                              const unsigned int &times) {
        /*
         0 1 2 3 4 5 6 7
         | | | | | | | | >>1
         3 0 1 2 7 4 5 6
         */
        
        vector<unsigned int> res(8,0);
        for (unsigned int i=0; i<4; i++) {
            unsigned int pos = (4+i-times)%4;
            res[i] = nodes[pos];
            res[i+4] = nodes[pos+4];
        }
        return res;
    
    }
	
}
