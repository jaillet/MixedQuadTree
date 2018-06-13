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

#ifndef QuadSurfTemplate_h
#define QuadSurfTemplate_h 1

#include <vector>
#include <iostream>
#include "HexRotation.h"

using std::vector;
using std::cout;

namespace Clobscode
{
class QuadSurfTemplate {
	
  public:
    
    QuadSurfTemplate();

    virtual ~QuadSurfTemplate();

    virtual bool one(const vector<unsigned int> &all, vector<bool> &in,
                     vector<vector<unsigned int> > &newsubs);

    virtual bool two(const vector<unsigned int> &all, vector<bool> &in,
                     vector<vector<unsigned int> > &newsubs);
    
    virtual bool three(const vector<unsigned int> &all, vector<bool> &in,
                       vector<vector<unsigned int> > &newsubs);
    
protected:
    
    virtual vector<unsigned int> rotate(const vector<unsigned int> &nodes,
                                        const unsigned int &times);

};
}
#endif
