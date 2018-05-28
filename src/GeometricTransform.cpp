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
* @file GeometricTransform.cpp
* @author Claudio Lobos, Fabrice Jaillet
* @version 0.1
* @brief
**/

#include "GeometricTransform.h"
#include <algorithm>

namespace Clobscode
{
    
    //--------------------------------------------------------------------------------
	//--------------------------------------------------------------------------------
    GeometricTransform::GeometricTransform()
        :x(0.0),y(0.0),z(0.0) {	} //centroid initialized by Point3D default constructor
    
    //--------------------------------------------------------------------------------
	//--------------------------------------------------------------------------------
    GeometricTransform::GeometricTransform(Point3D &p, double &x, double &y, double &z)
     :centroid(p),x(x),y(y),z(z) {	}
	
    //--------------------------------------------------------------------------------
	//--------------------------------------------------------------------------------
	GeometricTransform::~GeometricTransform(){
		
	}

    //--------------------------------------------------------------------------------
	//--------------------------------------------------------------------------------
    void GeometricTransform::rotatePolyline(Polyline &ply) {
        //calculateAnglesAndCentroid(ply);
		
        for (unsigned int i=0; i<ply.getPoints().size(); i++){
            apply(ply.getPoints()[i]);
        }
        ply.update();

//FJA rq: previous version simpler....
//        // edges are the same, no need to copy
//        vector<vector<unsigned int> > edges (ply.getEdges().size(), vector<unsigned int>(2, 0));
//        for (unsigned int i=0; i<ply.getEdges().size(); i++) {
//            edges[i] = ply.getEdges()[i].getPoints();
//        }
//        //creating a new Polyline with the new vertices rotated
//        Polyline ply2 (ply.getPoints(), edges);
//        ply = ply2;
        
    }
    
    
    //--------------------------------------------------------------------------------
	//--------------------------------------------------------------------------------
    void GeometricTransform::rotatePolylineInverse(Polyline &ply) {
        //to implement here
        std::cerr << "warning at GeometricTransform::rotateSurfaceMeshInverse(Polyline &ply)\n";
        std::cerr << "not implemented yet...\n";
    }
    
    //--------------------------------------------------------------------------------
	//--------------------------------------------------------------------------------
    void GeometricTransform::calculateAnglesAndCentroid(Polyline &input) {

        std::cerr << "GeometricTransform::calculateAnglesAndCentroid(Polyline &ply)\n";
        std::cerr << "!!!!!!!!!!!!!!! seems to be different in 3D... to be checked\n";

//        const vector<Point3D> &vertices= tm.getPoints();
        const vector<Point3D> &normals= input.getNormals();
        vector<Point3D> maxNormals;
        
        double alpha = 0.0;//77/57.2957795;
        double beta = 0.0;//35/57.2957795;
        double gamma = 0.0;//56/57.2957795;
        double epsilon = toRadians(5.), angle;
        double maxNorm = 0.0;
        int maxIndex = 0;
        
        maxNormals.push_back(normals[0]);
        //getting a distribution of mean maximum normals
        for (unsigned int i=1; i<normals.size(); i++){
            for (unsigned int j=0; j<maxNormals.size(); j++){
                angle = acos((normals[i]*maxNormals[j])/(normals[i].Norm()*maxNormals[j].Norm()));
                if(angle<epsilon){
                    maxNormals[j] += normals[i];
                    break;
                }
                else if(angle>toRadians(175.)){
                    maxNormals[j] += normals[i]*(-1);
                    break;
                }
                else if(j==maxNormals.size()-1){
                    maxNormals.push_back(normals[i]);
                    break;
                }
            }
        }
          
        //sorting the maxnormals from the greatest norm to lower one
        std::sort (maxNormals.begin(), maxNormals.end(), largestNormal);
        
        Point3D xAxis (1., 0., 0.);
        Point3D zAxis (0., 0., 1.);
        Point3D yAxis (0., 1., 0.);
        Point3D aux =  maxNormals[0];
        //finding a orthogonal maxnormal for the greatest maxnormal
        for (unsigned int j=1; j<maxNormals.size(); j++){
            if(maxNormals[j].Norm()/maxNormals[0].Norm() > 0.2){
                maxNorm = acos((maxNormals[0]*maxNormals[j])/(maxNormals[0].Norm()*maxNormals[j].Norm()));
                if(maxNorm > toRadians(80.) && maxNorm < toRadians(100.)){
                    maxIndex = j;
                    break;
                }
            }
        }


        //        //finding the rotation angles
        //        aux.X() = 0;
        //        alpha = acos(zAxis*aux/(zAxis.Norm()*aux.Norm()));
        //        if(aux.Y() < 0) alpha*=-1;

        //        maxNormals[0].xAxisRotation(alpha);

        //        aux =  maxNormals[0];
        //        aux.Y() = 0;
        //        beta = acos(zAxis*aux/(zAxis.Norm()*aux.Norm()));
        //        if(aux.X() > 0) beta*=-1;
        //        maxNormals[maxIndex].xAxisRotation(alpha);
        //        maxNormals[maxIndex].yAxisRotation(beta);
        //        aux = maxNormals[maxIndex];
        //        aux.Z() = 0;
        //        if(maxIndex==0 || acos(xAxis*aux/(xAxis.Norm()*aux.Norm()))==acos(-1)) {
        //            gamma=0;
        //        }
        //        else{
        //            gamma = acos(xAxis*aux/(xAxis.Norm()*aux.Norm()));
        //            if(aux.Y() > 0) {
        //                gamma*=-1;
        //            }
        //        }

        //finding the rotation angles
        // align maxNormal[0] on X, wrt Z
        aux.Z() = 0;
        gamma = acos(xAxis*aux/(xAxis.Norm()*aux.Norm()));
        if(aux.Y() > 0) gamma*=-1;

        maxNormals[0].zAxisRotation(gamma);

        // align maxNormal[0] on X, wrt Y
        aux =  maxNormals[0];
        aux.Y() = 0;
        beta = acos(xAxis*aux/(xAxis.Norm()*aux.Norm()));
        if(aux.X() < 0) beta*=-1;
        maxNormals[maxIndex].zAxisRotation(gamma);
        maxNormals[maxIndex].yAxisRotation(beta);
        // align maxNormals[maxIndex] on Y, wrt X
        aux = maxNormals[maxIndex];
        aux.X() = 0;
        if(maxIndex==0 || acos(yAxis*aux/(yAxis.Norm()*aux.Norm()))==acos(-1)) {
            alpha=0;
        }
        else{
            alpha = acos(yAxis*aux/(yAxis.Norm()*aux.Norm()));
            if(aux.Z() > 0) {
                alpha*=-1;
            }
        }

		this->x=alpha;
		this->y=beta;
		this->z=gamma;
        this->centroid=input.getCentroid();
	}
    //--------------------------------------------------------------------------------
	//--------------------------------------------------------------------------------
    void GeometricTransform::apply(Point3D &p) {
        p -= centroid;
        p.xAxisRotation(x);
        p.yAxisRotation(y);
        p.zAxisRotation(z);
    }
    
    
    //--------------------------------------------------------------------------------
	//--------------------------------------------------------------------------------
    void GeometricTransform::applyInverse(Point3D &p) {
       
        p.zAxisRotation(-z);
        p.yAxisRotation(-y);
        p.xAxisRotation(-x);
        p += centroid;
    }
}
