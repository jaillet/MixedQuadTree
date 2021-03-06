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
* @file Visitor.h
* @author Claudio Lobos, Fabrice Jaillet
* @version 0.1
* @brief
**/

#ifndef Visitor_h
#define Visitor_h 1

/*
 * Parent definition of visitors
 * Visitors shall extend this main class to get forward declarations,
 * and also for operations that extends to more than one class...
 */

namespace Clobscode {

    //forward declarations
    class Quadrant;

    class Visitor {

    public:
        Visitor() { }

        virtual bool visit(Quadrant *q) { return false; }
    };


}

#endif //Visitor_h
