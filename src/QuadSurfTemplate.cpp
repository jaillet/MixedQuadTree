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
	
	bool QuadSurfTemplate::one(const vector<unsigned int> &nodes, vector<bool> &in,
                               vector<vector<unsigned int> > &newsubs){
		
        //cout << "\n\n QuadSurfTemplate::one\n";
        
        vector<unsigned int> rotated = nodes;
        
        unsigned int rotation = 0;
        for (unsigned int i=0; i<4; i++) {
            if (in[i]) {
                break;
            }
            rotation++;
        }
        
        if (rotation!=0) {
            rotated = rotate(nodes,rotation);
        }
        
        newsubs.clear();
        newsubs.reserve(1);
        vector<unsigned int> t(3,0);
		
		t[0] = rotated[0];
		t[1] = rotated[1];
		t[2] = rotated[3];
		
		newsubs.push_back(t);
		
		return true;
	}
    
    bool QuadSurfTemplate::two(const vector<unsigned int> &nodes, vector<bool> &in,
                               vector<vector<unsigned int> > &newsubs){
        
        //cout << "\n\n QuadSurfTemplate::two\n";
        vector<unsigned int> rotated = nodes;
        
        unsigned int rotation = 0;
        for (unsigned int i=0; i<4; i++) {
            if (in[i]) {
                break;
            }
            rotation++;
        }
        
        //case: two consecutive nodes in: we return the same element.
        if (in[(rotation+1)%4] || in[(rotation+3)%4]) {
            newsubs.push_back(nodes);
            return true;
        }
        
        if (rotation!=0) {
            rotated = rotate(nodes,rotation);
        }
        
        newsubs.clear();
        newsubs.reserve(2);
        vector<unsigned int> t1(3,0), t2(3,0);
        
        t1[0] = rotated[0];
        t1[1] = rotated[1];
        t1[2] = rotated[3];
        
        t2[0] = rotated[1];
        t2[1] = rotated[2];
        t2[2] = rotated[3];
        
        newsubs.push_back(t1);
        newsubs.push_back(t2);
        
        return true;
    }
    
    bool QuadSurfTemplate::three(const vector<unsigned int> &nodes, vector<bool> &in,
                                 vector<vector<unsigned int> > &newsubs){
        
        //cout << "\n\n QuadSurfTemplate::three\n";
        vector<unsigned int> rotated = nodes;
        
        unsigned int rotation = 0;
        for (unsigned int i=0; i<4; i++) {
            if (!in[(i+2)%4]) {
                break;
            }
            rotation++;
        }
        
        if (rotation!=0) {
            rotated = rotate(nodes,rotation);
        }
        
        newsubs.clear();
        newsubs.reserve(2);

        vector<unsigned int> t1(3,0), t2(3,0);
        
        t1[0] = rotated[0];
        t1[1] = rotated[2];
        t1[2] = rotated[3];
        
        t2[0] = rotated[0];
        t2[1] = rotated[1];
        t2[2] = rotated[2];
        
        newsubs.push_back(t1);
        newsubs.push_back(t2);
        
        return true;

    }
    
    vector<unsigned int> QuadSurfTemplate::rotate(const vector<unsigned int> &nodes,
                                                  const unsigned int &times) {
        
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
