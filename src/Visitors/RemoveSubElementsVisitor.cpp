/*
 <Mix-mesher: region type. This program generates a mixed-elements 2D mesh>

 Copyright (C) <2013,2018>  <Claudio Lobos> All rights reserved.

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
* @file RemoveSubElementsVisitor.cpp
* @author Claudio Lobos, Fabrice Jaillet
* @version 0.1
* @brief
**/

#include "RemoveSubElementsVisitor.h"
#include "../Quadrant.h"

namespace Clobscode
{

    RemoveSubElementsVisitor::RemoveSubElementsVisitor() :points(NULL) {
    }

    void RemoveSubElementsVisitor::setPoints(vector<MeshPoint> &points) {
        this->points = &points;
    }

    bool RemoveSubElementsVisitor::visit(Quadrant *q) {

        vector<vector<unsigned int>> &sub_elements = q->sub_elements;

        list<vector<unsigned int> > still_in;
        list<vector<unsigned int> >::const_iterator iter;

        for (unsigned int i=0; i<sub_elements.size(); i++) {

            bool onein = false;
            const vector<unsigned int> &e_pts = sub_elements[i];
            // for all vertices in subelements
            for (unsigned int j=0; j<e_pts.size(); j++) {
                // check if one at least is inside
                if (points->at(e_pts[j]).isInside()) {
                    onein = true;
                    break;
                }
            }
            if (onein) {
                still_in.push_back(sub_elements[i]);
            }
        }

        // at least one vertex per subelem is in
        if (still_in.size()==sub_elements.size()) {
            return false;
        }
        // or none is inside
        if (still_in.empty()) {
            return true;
        }

        sub_elements.clear();
        sub_elements.insert(sub_elements.begin(),still_in.begin(),still_in.end());
//        sub_elements.reserve(still_in.size());
//        for (iter=still_in.begin(); iter!=still_in.end(); iter++) {
//            sub_elements.push_back(*iter);
//        }
        return false;
    }
}
