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

#include "QuadSurfTemplate.h"

namespace Clobscode
{
    QuadSurfTemplate::QuadSurfTemplate(){
        
    }
    
    QuadSurfTemplate::~QuadSurfTemplate(){
        
    }
    
    const bool QuadSurfTemplate::one(const vector<unsigned int> &nodes, vector<bool> &in,
                                     vector<vector<unsigned int> > &newsubs) const {
        
        vector<unsigned int> rotated=nodes;
        
        unsigned int rotation = 0;
        for (unsigned int i=0; i<4; i++) {
            if (in[i]) {
                break;
            }
            rotation++;
        }
        
        if (rotation!=0) {
            std::rotate(rotated.begin(),rotated.begin()+rotation,rotated.end());
        }
        
        newsubs.resize(1);
        newsubs[0]={rotated[0],rotated[1],rotated[3]};
        
        return true;
    }
    
    const bool QuadSurfTemplate::two(const vector<unsigned int> &nodes, vector<bool> &in,
                                     vector<vector<unsigned int> > &newsubs) const {
        
        vector<unsigned int> rotated=nodes;
        
        unsigned int rotation = 0;
        for (unsigned int i=0; i<4; i++) {
            if (in[i]) {
                break;
            }
            rotation++;
        }
        
        //case: two consecutive nodes in: we return the same element.
        if (in[(rotation+1)%4] || in[(rotation+3)%4]) {
            
            //new subs is not the original Quadrant sub elements.
            newsubs.push_back(nodes);
            return true;
        }
        
        if (rotation!=0) {
            //          rotated = rotate(nodes,rotation);
            std::rotate(rotated.begin(),rotated.begin()+rotation,rotated.end());
        }
        
        newsubs.resize(2);
        newsubs[0]= {rotated[0],rotated[1],rotated[3]} ;
        newsubs[1]= {rotated[1],rotated[2],rotated[3]} ;
        
        return true;
    }
    
    const bool QuadSurfTemplate::three(const vector<unsigned int> &nodes, vector<bool> &in,
                                       vector<vector<unsigned int> > &newsubs) const {
        
        vector<unsigned int> rotated= nodes;
        
        unsigned int rotation = 0;
        for (unsigned int i=0; i<4; i++) {
            if (!in[(i+2)%4]) {
                break;
            }
            rotation++;
        }
        
        
        if (rotation!=0) {
            std::rotate(rotated.begin(),rotated.begin()+rotation,rotated.end());
        }
        
        newsubs.resize(2);
        newsubs[0]= {rotated[0],rotated[2],rotated[3]} ;
        newsubs[1]= {rotated[0],rotated[1],rotated[2]} ;
        
        return true;
        
    }
    
    vector<unsigned int> QuadSurfTemplate::rotated(const vector<unsigned int> &nodes,
                                                   const unsigned int &times) {
        
        std::cerr << "Obsolete: replaced by std::rotate(),"
        <<   " Since no need to copy the vector in QuadSurfTemplate::one(), two(), three() ! \n";
        
        /*
         0 1 2 3
         | | | | >>1
         1 2 3 0
         */
        
        vector<unsigned int> res(4,0);
        for (unsigned int i=0; i<4; i++) {
            unsigned int pos = (i+times)%4;
            res[i] = nodes[pos];
        }
        return res;
    }
}
