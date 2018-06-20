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
* @file Mesher.cpp
* @author Claudio Lobos, Fabrice Jaillet
* @version 0.1
* @brief
**/

#include "Mesher.h"
#include <math.h>

namespace Clobscode
{
//vector<MeshPoint> points;
//vector<Quadrant> Quadrants;
//set<QuadEdge> QuadEdges;
//list<RefinementRegion *> regions;

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    Mesher::Mesher(){

    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------

    Mesher::~Mesher(){

    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------

    //Create a grid mesh regarding the Bounding Box of input mesh.
    //This will produce several cubes as roots of an octree structure.
    //Then split each initial element 8^rl times (where rl stands
    //for Refinement Level).
    FEMesh Mesher::refineMesh(Polyline &input, const unsigned short &rl,
                              const string &name, list<unsigned int> &roctli,
                              list<RefinementRegion *> &all_reg,
                              GeometricTransform &gt, const unsigned short &minrl,
                              const unsigned short &omaxrl, bool decoration){

        //Note: rotation are not enabled when refining an already produced mesh.
        bool rotated = !gt.Default();
        if(rotated) {
            /*cout << "rotating input surface mesh\n";
            cout << gt.getCentroid() << "\n";
            cout << gt.getXAxis() << " " << gt.getYAxis();
            cout << " " << gt.getZAxis() << "\n";*/
            gt.rotatePolyline(input);
        }

        //split Quadrants until the refinement level (rl) is achieved.
        //The output will be a one-irregular mesh.
        splitQuadrants(rl,input,roctli,all_reg,name,minrl,omaxrl);

        //The points of the Quadrant mesh must be saved at this point, otherwise node are
        //projected onto the surface and causes further problems with knowing if nodes
        //are inside, outside or projected.
        const vector<MeshPoint> &oct_points = points;

        //link element and node info for code optimization.
        /*
        detectFeatureQuadrants(input);
        linkElementsToNodes();
        detectInsideNodes(input);

        projectCloseToBoundaryNodes(input);
        removeOnSurface(input);

        //apply the surface Patterns
        applySurfacePatterns(input);
        removeOnSurface(input);*/

        //Now that we have all the elements, we can save the Quadrant mesh.
        unsigned int nels = Quadrants.size();
        Services::WriteQuadtreeMesh(name,oct_points,Quadrants,QuadEdges,nels,gt);

        /*detectInsideNodes(input);

        //update element and node info.
        linkElementsToNodes();

        //shrink outside nodes to the input domain boundary
        shrinkToBoundary(input);*/

        if (rotated) {
            for (unsigned int i=0; i<points.size(); i++) {
                gt.applyInverse(points[i].getPoint());
            }
        }

        //the almighty output mesh
        FEMesh mesh;

        //save the data of the mesh
        saveOutputMesh(mesh,decoration);

        return mesh;
    }


    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------

    //Create a grid mesh regarding the Bounding Box of input mesh.
    //This will produce several cubes as roots of an octree structure.
    //Then split each initial element 8^rl times (where rl stands
    //for Refinement Level).
    FEMesh Mesher::generateMesh(Polyline &input, const unsigned short &rl,
                                const string &name, list<RefinementRegion *> &all_reg, bool decoration){

        //ATTENTION: geometric transform causes invalid input rotation when the
        //input is a cube.
        GeometricTransform gt;

        //rotate: This method is written below and its mostly commented because
        //it causes conflicts when the input is a cube. Must be checked.
        bool rotated = rotateGridMesh(input, all_reg, gt);

        //generate root Quadrants
        generateGridMesh(input);

        //split Quadrants until the refinement level (rl) is achieved.
        //The output will be a one-irregular mesh.
        generateQuadtreeMesh(rl,input,all_reg,name);

        //The points of the Quadrant mesh must be saved at this point, otherwise node are
        //projected onto the surface and causes further problems with knowing if nodes
        //are inside, outside or projected.
        vector<MeshPoint> oct_points = points;

        //link element and node info for code optimization, also
        //detect Quadrants with features.
        detectFeatureQuadrants(input);
        linkElementsToNodes();
        detectInsideNodes(input);

        projectCloseToBoundaryNodes(input);
        //removeOnSurface(input);
        
        //Now that we have all the elements, we can save the Quadrant mesh.
        unsigned int nels = Quadrants.size();
        Services::WriteQuadtreeMesh(name,points,Quadrants,QuadEdges,nels,gt);

        //update element and node info.
        linkElementsToNodes();
        detectInsideNodes(input);

        //shrink outside nodes to the input domain boundary
        shrinkToBoundary(input);
        
        //apply the surface Patterns
        applySurfacePatterns(input);
        removeOnSurface(input);

        if (rotated) {
            // rotate the mesh
            for (unsigned int i=0; i<points.size(); i++) {
                gt.applyInverse(points[i].getPoint());
            }
            // rotate back the polyline as well
            for (unsigned int i=0; i<input.getPoints().size(); i++) {
                gt.applyInverse(input.getPoints()[i]);
            }
        }

        //the almighty output mesh
        FEMesh mesh;

        //save the data of the mesh in its final state
        saveOutputMesh(mesh,decoration);
        
        
        //Debugging:
        /*vector<unsigned int> color (Quadrants.size(),0);
        for (unsigned int i=0; i<Quadrants.size(); i++) {
            
            if (Quadrants[i].isSurface()) {
                
                if (input.hasFeature()) {
                    color[i]=1;
                }
            }
        }
        mesh.setColoredCells(color);*/

        return mesh;
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    bool Mesher::rotateGridMesh(Polyline &input, list<RefinementRegion *> &all_reg,
                                GeometricTransform &gt){

        list<RefinementRegion *>::const_iterator it, rrot;
        bool inputHasbeenRotated = false;

        for (it = all_reg.begin(); it!=all_reg.end(); it++) {
            //in case of input roi
            if((*it)->needsInputRotation()){
                if(!inputHasbeenRotated){
                    gt = (*it)->rotateWithinYou(input);
                    inputHasbeenRotated = true;
                    rrot=it;
                    break;
                }
            }
        }

        if (inputHasbeenRotated) {
            for (it = all_reg.begin(); it!=all_reg.end(); it++) {
                if (it!=rrot) {
                    if ((*it)->needsLocalRotation()) {
                        (*it)->rotate(gt);
                    }
                }
            }
        }

        return inputHasbeenRotated;
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------

    void Mesher::generateGridMesh(Polyline &input){
        //vectors with each coordinate per axis
        vector<double> all_x, all_y;
        vector<vector<unsigned int> > elements;

        GridMesher gm;
        gm.generatePoints(input.getBounds(),all_x,all_y);
        gm.generateMesh(all_x,all_y,points,elements);

        Quadrants.reserve(elements.size());

        //create the root Quadrants
        for (unsigned int i=0; i<elements.size(); i++) {
            Quadrant o (elements[i], 0);
            //Only when the Quadrant intersects the input
            //add it to the list of current Quadrants. As
            //This is the first time Quadrants are checked
            //for intersections they must be made w.r.t.
            //all input edges.
            IntersectionsVisitor iv(false);
            //if (o.checkIntersections(input,points)) {
            iv.setPolyline(input);
            iv.setPoints(points);
            if (o.accept(&iv)) {
                EdgeVisitor::insertEdges(&o, QuadEdges);
                Quadrants.push_back(o);
            }
        }
        /* debug log
        std::cerr << "QuadEdges cont:";
        for (auto it=QuadEdges.begin(); it!=QuadEdges.end(); ++it)
            std::cerr << " -" << *it;
        std::cerr << '\n'<< std::flush;
        std::cerr << "Quadrants contains:";
        for (auto it=Quadrants.begin(); it!=Quadrants.end(); ++it)
            std::cerr << " -" << *it;
        std::cerr << '\n'<< std::flush; */

    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------

    void Mesher::splitQuadrants(const unsigned short &rl, Polyline &input,
                                list<unsigned int> &roctli,
                                list<RefinementRegion *> &all_reg, const string &name,
                                const unsigned short &minrl, const unsigned short &omaxrl){

        //to save m3d files per stage
        Clobscode::Services io;

        //list of temp Quadrants
        list<Quadrant> tmp_Quadrants, new_Quadrants;
        //list of the points added at this refinement iteration:
        list<Point3D> new_pts;

        list<Quadrant>::iterator iter;

        for (unsigned int i=0; i<Quadrants.size(); i++) {
            tmp_Quadrants.push_back(Quadrants[i]);
        }

        //create visitors and give them variables
        SplitVisitor sv;
        sv.setPoints(points);
        sv.setEdges(QuadEdges);
        sv.setNewPts(new_pts);

        //----------------------------------------------------------
        //refine once each Quadrant in the list
        //----------------------------------------------------------
        unsigned int cindex = 0;
        
        if (!roctli.empty()) {
            list<unsigned int>::iterator octidx = roctli.begin();
            
            for (iter=tmp_Quadrants.begin(); iter!=tmp_Quadrants.end(); ++iter) {
                if (octidx!=roctli.end() && cindex==(*octidx)) {
                    
                    octidx++;
                    unsigned short orl = (*iter).getRefinementLevel();
                    
                    list<unsigned int> inter_faces = iter->getIntersectedEdges();
                    
                    vector<vector<Point3D> > clipping_coords;
                    sv.setClipping(clipping_coords);
                    
                    vector<vector<unsigned int> > split_elements;
                    sv.setNewEles(split_elements);
                    
                    iter->accept(&sv);
                    
                    if (inter_faces.empty()) {
                        for (unsigned int j=0; j<split_elements.size(); j++) {
                            
                            Quadrant o (split_elements[j],orl+1);
                            new_Quadrants.push_back(o);
                        }
                    }
                    else {
                        for (unsigned int j=0; j<split_elements.size(); j++) {
                            Quadrant o (split_elements[j],orl+1);
                            
                            //the new points are inserted in bash at the end of this
                            //iteration. For this reason, the coordinates must be passed
                            //"manually" at this point (clipping_coords).
                            IntersectionsVisitor iv(true);
                            //if (o.checkIntersections(input,inter_faces,clipping_coords[j]))
                            iv.setPolyline(input);
                            iv.setEdges(inter_faces);
                            iv.setCoords(clipping_coords[j]);
                            
                            if (o.accept(&iv)) {
                                new_Quadrants.push_back(o);
                            }
                            else {
                                //The element doesn't intersect any input face.
                                //It must be checked if it's inside or outside.
                                //Only in the first case add it to new_Quadrants.
                                //Test this with parent Quadrant faces only.
                                
                                //Comment the following lines of this 'else' if
                                //only intersecting Quadrants are meant to be
                                //displayed.
                                
                                //note: inter_faces is quite enough to check if
                                //element is inside input, no Quadrant needed,
                                //so i moved the method to mesher  --setriva
                                
                                if (isItIn(input,inter_faces,clipping_coords[j])) {
                                    new_Quadrants.push_back(o);
                                }
                            }
                        }
                    }
                }
                else {
                    new_Quadrants.push_back(*iter);
                }
                cindex++;
            }
            
            
            tmp_Quadrants.clear();
            tmp_Quadrants = new_Quadrants;
            new_Quadrants.clear();
            
            //add the new points to the vector
            list<Point3D>::iterator piter;
            points.reserve(points.size() + new_pts.size());
            for (piter=new_pts.begin(); piter!=new_pts.end(); ++piter) {
                points.push_back(MeshPoint (*piter));
            }
        }

        //----------------------------------------------------------
        //refine each Quadrant until the Refinement Level is reached
        //----------------------------------------------------------
        
        for (unsigned short i=minrl; i<rl; i++) {

            bool changed = false;
            unsigned int refoct_int = 0, refoct_ext = 0;

            //the new_pts is a list that holds the coordinates of
            //new points inserted at this iteration. At the end of
            //this bucle, they are inserted in the point vector
            new_pts.clear();

            list<RefinementRegion *>::iterator reg_iter;

            //split the Quadrants as needed
            for (iter=tmp_Quadrants.begin(); iter!=tmp_Quadrants.end(); ++iter) {

                bool to_refine = false;

                for (reg_iter=all_reg.begin(); reg_iter!=all_reg.end(); ++reg_iter) {

                    unsigned short region_rl = (*reg_iter)->getRefinementLevel();
                    if (region_rl<i) {
                        continue;
                    }

                    //If the Quadrant has a greater RL than the region needs, continue
                    if (region_rl<=(*iter).getRefinementLevel()) {
                        continue;
                    }

                    //Get the two extreme nodes of the Quadrant to test intersection with
                    //this RefinementRegion. If not, conserve it as it is.
                    //unsigned int n_idx1 = (*iter).getPoints()[0];
                    //unsigned int n_idx2 = (*iter).getPoints()[6];

                    if ((*reg_iter)->intersectsQuadrant(points,*iter)) {
                        to_refine = true;
                        changed = true;
                    }
                }

                //now if refinement is not needed, we add the Quadrant as it was.
                if (!to_refine) {
                    new_Quadrants.push_back(*iter);
                    continue;
                }
                else {
                    unsigned short orl = (*iter).getRefinementLevel();

                    list<unsigned int> inter_faces = iter->getIntersectedEdges();

                    vector<vector<Point3D> > clipping_coords;
                    sv.setClipping(clipping_coords);

                    vector<vector<unsigned int> > split_elements;
                    sv.setNewEles(split_elements);

                    iter->accept(&sv);

                    if (inter_faces.empty()) {

                        refoct_int++;


                        for (unsigned int j=0; j<split_elements.size(); j++) {

                            Quadrant o (split_elements[j],orl+1);
                            new_Quadrants.push_back(o);
                        }
                        //break;
                    }
                    else {

                        refoct_ext++;

                        for (unsigned int j=0; j<split_elements.size(); j++) {
                            Quadrant o (split_elements[j],orl+1);
                            //the new points are inserted in bash at the end of this
                            //iteration. For this reason, the coordinates must be passed
                            //"manually" at this point (clipping_coords).
                            IntersectionsVisitor iv(true);
                            //if (o.checkIntersections(input,inter_faces,clipping_coords[j]))
                            iv.setPolyline(input);
                            iv.setEdges(inter_faces);
                            iv.setCoords(clipping_coords[j]);

                            if (o.accept(&iv)) {
                                new_Quadrants.push_back(o);
                            }
                            else {
                                //The element doesn't intersect any input face.
                                //It must be checked if it's inside or outside.
                                //Only in the first case add it to new_Quadrants.
                                //Test this with parent Quadrant faces only.

                                //Comment the following lines of this 'else' if
                                //only intersecting Quadrants are meant to be
                                //displayed.

                                //note: inter_faces is quite enough to check if
                                //element is inside input, no Quadrant needed,
                                //so i moved the method to mesher  --setriva

                                if (isItIn(input,inter_faces,clipping_coords[j])) {
                                    new_Quadrants.push_back(o);
                                }
                            }
                        }
                    }
                }
            }
            
            if (!changed) {
                continue;
            }

            //remove the old Quadrants
            tmp_Quadrants.clear();
            tmp_Quadrants = new_Quadrants;
            new_Quadrants.clear();

            //add the new points to the vector
            list<Point3D>::iterator piter;
            points.reserve(points.size() + new_pts.size());
            for (piter=new_pts.begin(); piter!=new_pts.end(); ++piter) {
                points.push_back(MeshPoint (*piter));
            }
        }

        //----------------------------------------------------------
        //produce a one-irregular mesh
        //----------------------------------------------------------

        bool one_irregular = false;
        new_Quadrants.clear();

        //visitante oneIrregular
        OneIrregularVisitor oiv;
        oiv.setEdges(QuadEdges);
        oiv.setMaxRefLevel(omaxrl);

        unsigned int irroct = 0;
        //cout << "number of regular Quadrants ";


        while (!one_irregular) {
            one_irregular = true;
            new_pts.clear();
            //refine until the mesh is one-irregular
            for (iter=tmp_Quadrants.begin(); iter!=tmp_Quadrants.end(); ++iter) {
                if (!(*iter).accept(&oiv)) {
                    //split this Quadrant
                    vector<vector<Point3D> > clipping_coords;
                    sv.setClipping(clipping_coords);
                    vector< vector <unsigned int> > split_elements;
                    sv.setNewEles(split_elements);

                    list<unsigned int> inter_faces = (*iter).getIntersectedEdges();

                    //(*iter).split(points,new_pts,QuadEdges,split_elements,clipping_coords);
                    //split the Quadrant
                    iter->accept(&sv);


                    unsigned short prl = (*iter).getRefinementLevel();
                    //insert the new elements
                    for (unsigned int j=0; j<split_elements.size(); j++) {
                        Quadrant o (split_elements[j],prl+1);

                        if (inter_faces.empty()) {
                            new_Quadrants.push_back(o);
                            continue;
                        }

                        //select_faces = true
                        IntersectionsVisitor iv(true);
                        //if (o.checkIntersections(input,inter_faces,clipping_coords[j]))
                        iv.setPolyline(input);
                        iv.setEdges(inter_faces);
                        iv.setCoords(clipping_coords[j]);

                        if (o.accept(&iv)) {
                            new_Quadrants.push_back(o);
                        }
                        else {
                            if (isItIn(input,inter_faces,clipping_coords[j])) {
                                new_Quadrants.push_back(o);
                            }
                        }
                    }
                    one_irregular = false;
                }
                else {
                    irroct++;
                    new_Quadrants.push_back(*iter);
                }
            }

            if (one_irregular) {
                break;
            }

            //remove the old Quadrants
            tmp_Quadrants.clear();
            tmp_Quadrants = new_Quadrants;
            new_Quadrants.clear();

            //if no points were added at this iteration, it is no longer
            //necessary to continue the refinement.
            if (new_pts.empty()) {
                break;
            }

            //add the new points to the vector
            list<Point3D>::iterator piter;
            points.reserve(points.size() + new_pts.size());
            for (piter=new_pts.begin(); piter!=new_pts.end(); ++piter) {
                points.push_back(MeshPoint (*piter));
            }
        }

        //----------------------------------------------------------
        // apply transition patterns
        //----------------------------------------------------------

        //clean tmp point list
        new_pts.clear();

        //TransitionPatternVisitor section 
        TransitionPatternVisitor tpv;
        tpv.setPoints(points);
        tpv.setEdges(QuadEdges);
        tpv.setMaxRefLevel(omaxrl);
        
        unsigned int cl3=0;

        for (iter = tmp_Quadrants.begin(); iter!=tmp_Quadrants.end(); ++iter) {

            if (!(*iter).accept(&tpv)) {
                std::cerr << "Error at Mesher::generateOctreeMesh";
                std::cerr << " Transition Pattern not found\n";
            }
        }

        //if no points were added at this iteration, it is no longer
        //necessary to continue the refinement.
        if (!new_pts.empty()) {
            //add the new points to the vector
            list<Point3D>::iterator piter;
            points.reserve(points.size() + new_pts.size());
            for (piter=new_pts.begin(); piter!=new_pts.end(); ++piter) {
                points.push_back(MeshPoint (*piter));
            }
        }

        /*{
            //save pure octree mesh
            FEMesh pure_octree;
            saveOutputMesh(pure_octree,points,tmp_Quadrants);
            string tmp_name = name + "_alloct";
            Services::WriteVTK(tmp_name,pure_octree);
            Services::WriteMixedVolumeMesh(tmp_name,pure_octree);
        }*/

        //put the Quadrants in a vector
        Quadrants.clear();
        Quadrants.reserve(tmp_Quadrants.size());
        for (iter=tmp_Quadrants.begin(); iter!=tmp_Quadrants.end(); ++iter) {
            Quadrants.push_back(*iter);
        }
    }


    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------

    void Mesher::generateQuadtreeMesh(const unsigned short &rl, Polyline &input,
                                    const list<RefinementRegion *> &all_reg,
                                    const string &name){

        //to save m3d files per stage
        Clobscode::Services io;

        //list of temp Quadrants
        list<Quadrant> tmp_Quadrants, new_Quadrants;
        //list of the points added at this refinement iteration:
        list<Point3D> new_pts;
        list<Quadrant>::iterator iter;

        //initialize list with vector, FJA: sure we want a list?
        tmp_Quadrants.assign(Quadrants.begin(),Quadrants.end());
//        for (unsigned int i=0; i<Quadrants.size(); i++) {
//            tmp_Quadrants.push_back(Quadrants[i]); //push_back is slow...
//        }

        //create visitors and give them variables
        SplitVisitor sv;
        sv.setPoints(points);
        sv.setEdges(QuadEdges);
        sv.setNewPts(new_pts);

        //----------------------------------------------------------
        //refine each Quadrant until the Refinement Level is reached
        //----------------------------------------------------------

        for (unsigned short i=0; i<rl; i++) {

            //the new_pts is a list that holds the coordinates of
            //new points inserted at this iteration. At the end of
            //this bucle, they are inserted in the point vector
            new_pts.clear();

            list<RefinementRegion *>::const_iterator reg_iter;

            //split the Quadrants as needed
            for (iter=tmp_Quadrants.begin(); iter!=tmp_Quadrants.end(); ++iter) {

                bool to_refine = false;

                for (reg_iter=all_reg.begin(); reg_iter!=all_reg.end(); ++reg_iter) {

                    unsigned short region_rl = (*reg_iter)->getRefinementLevel();
                    if (region_rl<i) {
                        continue;
                    }

                    //If the Quadrant has a greater RL than the region needs, continue
                    if (region_rl<=(*iter).getRefinementLevel()) {
                        continue;
                    }

                    //Get the two extreme nodes of the Quadrant to test intersection with
                    //this RefinementRegion. If not, conserve it as it is.
                    //unsigned int n_idx1 = (*iter).getPoints()[0];
                    //unsigned int n_idx2 = (*iter).getPoints()[2];

                    if ((*reg_iter)->intersectsQuadrant(points,*iter)) {
                        to_refine = true;
                    }
                }

                //now if refinement is not needed, we add the Quadrant as it was.
                if (!to_refine) {
                    new_Quadrants.push_back(*iter);
                    continue;
                }
                else {
                    list<unsigned int> inter_edges = iter->getIntersectedEdges();

                    vector<vector<Point3D> > clipping_coords;
                    sv.setClipping(clipping_coords);

                    vector<vector<unsigned int> > split_elements;
                    sv.setNewEles(split_elements);

                    //iter->split(points,new_pts,QuadEdges,split_elements,clipping_coords);
                    //cout << "Accept" << endl;
                    iter->accept(&sv);

                    if (inter_edges.empty()) {
                        for (unsigned int j=0; j<split_elements.size(); j++) {
                            Quadrant o (split_elements[j],i+1);
                            new_Quadrants.push_back(o);
                        }
                    }
                    else {
                        for (unsigned int j=0; j<split_elements.size(); j++) {
                            Quadrant o (split_elements[j],i+1);
                            //the new points are inserted in bash at the end of this
                            //iteration. For this reason, the coordinates must be passed
                            //"manually" at this point (clipping_coords).

                            //select_faces = true
                            IntersectionsVisitor iv(true);
                            //if (o.checkIntersections(input,inter_faces,clipping_coords[j]))
                            iv.setPolyline(input);
                            iv.setEdges(inter_edges);
                            iv.setCoords(clipping_coords[j]);

                            if (o.accept(&iv)) {
                                new_Quadrants.push_back(o);
                            }
                            else {
                                //The element doesn't intersect any input face.
                                //It must be checked if it's inside or outside.
                                //Only in the first case add it to new_Quadrants.
                                //Test this with parent Quadrant faces only.

                                //Comment the following lines of this 'else' if
                                //only intersecting Quadrants are meant to be
                                //displayed.

                                //note: inter_edges is quite enough to check if
                                //element is inside input, no Quadrant needed,
                                //so i moved the method to mesher  --setriva

                                if (isItIn(input,inter_edges,clipping_coords[j])) {
                                    new_Quadrants.push_back(o);
                                }
                            }
                        }
                    }
                }
            }

            //remove the old Quadrants
            tmp_Quadrants.clear();
            tmp_Quadrants = new_Quadrants;
            new_Quadrants.clear();

            //if no points were added at this iteration, it is no longer
            //necessary to continue the refinement.
            if (new_pts.empty()) {
                cout << "warning at Mesher::generateOctreeMesh no new points!!!\n";
                break;
            }
            //add the new points to the vector
            list<Point3D>::const_iterator piter;
            points.reserve(points.size() + new_pts.size());
            for (piter=new_pts.begin(); piter!=new_pts.end(); ++piter) {
                points.push_back(MeshPoint (*piter));
            }
        }

        //----------------------------------------------------------
        //produce a one-irregular mesh
        //----------------------------------------------------------

        //cout << "  > producing one-irregular mesh ...";
        //cout.flush();

        bool one_irregular = false;
        new_Quadrants.clear();

        //visitante oneIrregular
        OneIrregularVisitor oiv;
        oiv.setEdges(QuadEdges);
        oiv.setMaxRefLevel(rl);

        /* //if pointMoved is ever needed, uncomment this:
         * PointMovedVisitor pmv;
         * pmv.setPoints(points);
         * pmv.setEdges(QuadEdges);
         * pmv.setMaxRefLevel(rl);
         * //haven't tested it, but if it's like the others,
         * //it should work right out of the box
         */
        while (!one_irregular) {
            one_irregular = true;
            new_pts.clear();
            //refine until the mesh is one-irregular
            for (iter=tmp_Quadrants.begin(); iter!=tmp_Quadrants.end(); ++iter) {
                if (!(*iter).accept(&oiv)){
                    //split this Quadrant
                    vector<vector<Point3D> > clipping_coords;
                    sv.setClipping(clipping_coords);
                    vector< vector <unsigned int> > split_elements;
                    sv.setNewEles(split_elements);

                    list<unsigned int> inter_faces = (*iter).getIntersectedEdges();

                    //(*iter).split(points,new_pts,QuadEdges,split_elements,clipping_coords);
                    //split the Quadrant
                    iter->accept(&sv);


                    unsigned short prl = (*iter).getRefinementLevel();
                    //insert the new elements
                    for (unsigned int j=0; j<split_elements.size(); j++) {
                        Quadrant o (split_elements[j],prl+1);

                        if (inter_faces.empty()) {
                            new_Quadrants.push_back(o);
                            continue;
                        }

                        //select_faces = true
                        IntersectionsVisitor iv(true);
                        //if (o.checkIntersections(input,inter_faces,clipping_coords[j]))
                        iv.setPolyline(input);
                        iv.setEdges(inter_faces);
                        iv.setCoords(clipping_coords[j]);

                        if (o.accept(&iv)) {
                            new_Quadrants.push_back(o);
                        }
                        else {
                            if (isItIn(input,inter_faces,clipping_coords[j])) {
                                new_Quadrants.push_back(o);
                            }
                        }
                    }
                    one_irregular = false;
                }
                else {
                    new_Quadrants.push_back(*iter);
                }
            }

            if (one_irregular) {
                break;
            }

            //remove the old Quadrants
            tmp_Quadrants.clear();
            tmp_Quadrants = new_Quadrants;
            new_Quadrants.clear();

            //if no points were added at this iteration, it is no longer
            //necessary to continue the refinement.
            if (new_pts.empty()) {
                break;
            }
            //add the new points to the vector
            points.reserve(points.size() + new_pts.size());
            // append new_points to points
            points.insert(points.end(),new_pts.begin(),new_pts.end());
        }

         //----------------------------------------------------------
         // apply transition patterns
         //----------------------------------------------------------

         //visitante TransitionPattern in check mode
         //(we'll reuse it in apply mode later)
         TransitionPatternVisitor tpv;
         tpv.setPoints(points);
         tpv.setEdges(QuadEdges);
         tpv.setMaxRefLevel(rl);
         new_pts.clear();

        for (iter = tmp_Quadrants.begin(); iter!=tmp_Quadrants.end(); ++iter) {
            if (!(*iter).accept(&tpv)) {
                std::cerr << "Error at Mesher::generateOctreeMesh";
                std::cerr << " Transition Pattern not found\n";
            }
        }

        //if no points were added at this iteration, it is no longer
        //necessary to continue the refinement.
        if (!new_pts.empty()) {
            //add the new points to the vector
            points.insert(points.end(),new_pts.begin(),new_pts.end());
        }

        //put the Quadrants in a vector
        Quadrants.assign(tmp_Quadrants.begin(),tmp_Quadrants.end());
    }
    
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------

    bool Mesher::isItIn(const Polyline &mesh, const list<unsigned int> &faces, const vector<Point3D> &coords) const {
        //this method is meant to be used by Quadrants that don't
        //intersect input domains. If they are inside of at least
        //one input mesh, then they must remain in the output mesh.
        bool first = mesh.pointIsInMesh(coords[0],faces);
        bool second = mesh.pointIsInMesh(coords[1],faces);
        if (first==second) {
            return first;
        }

        //cout << "one inconsistency detected -> hard test\n";
        //return mesh.pointIsInMesh(coords[0],faces);
        return mesh.pointIsInMesh(coords[0]);
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------

    unsigned int Mesher::saveOutputMesh(FEMesh &mesh,bool decoration){

        vector<vector<unsigned int> > out_els;
        vector<unsigned short > out_els_ref_level;
        vector<double > out_els_min_angle;

        //new_idxs will hold the index of used nodes in the outside vector for points.
        //If the a node is not used by any element, its index will be 0 in this vector,
        //therefore the actual index is shiffted in 1. In other words, node 0 is node 1,
        //and node n is node n+1.
        vector<unsigned int> new_idxs (points.size(),0);
        unsigned int out_node_count = 0;
        vector<Point3D> out_pts;

        //recompute node indexes and update elements with them.
        for (unsigned int i=0; i<Quadrants.size(); i++) {
            vector<vector<unsigned int> > sub_els= Quadrants[i].getSubElements();
            for (unsigned int j=0; j<sub_els.size(); j++) {

                vector<unsigned int> sub_ele_new_idxs = sub_els[j];
                for (unsigned int k=0; k<sub_ele_new_idxs.size();k++) {

                    unsigned int p_idx = sub_ele_new_idxs[k];

                    if (new_idxs[p_idx]==0) {
                        sub_ele_new_idxs[k] = out_node_count++;
                        new_idxs[p_idx]=out_node_count;
                        out_pts.push_back(points[p_idx].getPoint());
                    }
                    else {
                        sub_ele_new_idxs[k] = new_idxs[p_idx]-1;
                    }
                }
                if (decoration) {
                    //refinment level herited from quad
                    out_els_ref_level.push_back(Quadrants[i].getRefinementLevel());
                    //compute minAngle
                    int np=sub_ele_new_idxs.size(); //nb points of the element
                    double minAngle=std::numeric_limits<double>::infinity();
                    for (int i=0; i<np; ++i) {

                        const Point3D &P0 = out_pts[sub_ele_new_idxs[(i-1+np)%np]];
                        const Point3D &P1 = out_pts[sub_ele_new_idxs[i]];
                        const Point3D &P2 = out_pts[sub_ele_new_idxs[(i+1)%np]];

                        minAngle=std::min(minAngle, P1.angle3Points(P0,P2));
                    }
                    out_els_min_angle.push_back(minAngle);
                }
                out_els.push_back(sub_ele_new_idxs);
            }
        }

        mesh.setPoints(out_pts);
        mesh.setElements(out_els);
        mesh.setRefLevels(out_els_ref_level);
        mesh.setMinAngles(out_els_min_angle);
        return out_els.size();
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------

    unsigned int Mesher::saveOutputMesh(FEMesh &mesh,vector<MeshPoint> &tmp_points,
                                list<Quadrant> &tmp_Quadrants){

        vector<Point3D> out_pts;
        list<vector<unsigned int> > tmp_elements;
        vector<vector<unsigned int> > out_els;

        unsigned int n = tmp_points.size();
        out_pts.reserve(n);
        for (unsigned int i=0; i<n; i++) {
            out_pts.push_back(points[i].getPoint());
        }

        list<Quadrant>::iterator o_iter;

        for (o_iter=tmp_Quadrants.begin(); o_iter!=tmp_Quadrants.end(); ++o_iter) {

            vector<vector<unsigned int> > sub_els= o_iter->getSubElements();
            for (unsigned int j=0; j<sub_els.size(); j++) {
                tmp_elements.push_back(sub_els[j]);
            }
        }

        out_els.reserve(tmp_elements.size());
        list<vector<unsigned int> >::iterator e_iter;

        for (e_iter=tmp_elements.begin(); e_iter!=tmp_elements.end(); ++e_iter) {
            out_els.push_back(*e_iter);
        }

        mesh.setPoints(out_pts);
        mesh.setElements(out_els);
        return out_els.size();
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    void Mesher::detectFeatureQuadrants(Polyline &input) {
        
        unsigned int featCount = 0;
        for (unsigned int i=0; i<Quadrants.size(); i++) {
            if (input.hasFeature(Quadrants[i],points)) {
                Quadrants[i].setFeature();
                featCount++;
            }
        }
        cout << "number of quadrants with features: " << featCount << "\n";
    }
    
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    void Mesher::linkElementsToNodes(){
        //clear previous information
        for (unsigned int i=0; i<points.size(); i++) {
            points[i].clearElements();
        }

        //link element info to nodes
        for (unsigned int i=0; i<Quadrants.size(); i++) {
            
            const vector <unsigned int> &q_indpts = Quadrants[i].getPointIndex();

            for (unsigned int j=0; j<q_indpts.size(); j++) {
                points.at(q_indpts[j]).addElement(i);
            }
        }
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------

    void Mesher::detectInsideNodes(Polyline &input){
        for (unsigned int i=0; i<points.size(); i++) {
            if (points[i].wasOutsideChecked()) {
                continue;
            }

            list<unsigned int> p_eles = points[i].getElements(), p_edges;
            points[i].outsideChecked();
            if (p_eles.empty()) {
                continue;
            }
            list<unsigned int>::const_iterator iter;
            for (iter=p_eles.begin(); iter!=p_eles.end(); ++iter) {
                const list<unsigned int> &qedges= Quadrants[*iter].getIntersectedEdges();
                list<unsigned int>::const_iterator qe_iter;
                if (qedges.empty()) {
                    continue;
                }
                // append qedges to p_edges
                p_edges.insert(p_edges.end(),qedges.begin(),qedges.end());
            }

            p_edges.sort();
            p_edges.unique();

            // p_edges: edges intersected
            if (p_edges.empty() || input.pointIsInMesh(points[i].getPoint(),p_edges)) {
                points[i].setInside();
            }
        }
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------

    void Mesher::removeOnSurface(Polyline &input){

        list<Quadrant> newele,removed;
        RemoveSubElementsVisitor rsv;
        rsv.setPoints(points);
        //remove elements without an inside node.
        for (unsigned int i=0; i<Quadrants.size(); i++) {
            if (Quadrants[i].isInside()) {
                newele.push_back(Quadrants[i]);
                continue;
            }
        
            if (Quadrants[i].hasFeature()) {
                if (Quadrants[i].accept(&rsv)) {
                    //Quadrant with feature to be removed
                    //only if average node is outside the
                    //the edeges it intersects.
                    Point3D avg;
                    for (auto quaNoIdx:Quadrants[i].getPointIndex()) {
                        avg+=points[quaNoIdx].getPoint();
                    }
                    avg/=Quadrants[i].getPointIndex().size();
                    if (input.pointIsInMesh(avg,Quadrants[i].getIntersectedEdges())) {
                        newele.push_back(Quadrants[i]);
                    }
                    else {
                        removed.push_back(Quadrants[i]);
                    }
                }
                else {
                    newele.push_back(Quadrants[i]);
                }
            }
            else { //FJA add a "else" here as some quadrants are inserted twice

                //if (Quadrants[i].removeOutsideSubElements(points)) {
                if (Quadrants[i].accept(&rsv)) {
                    removed.push_back(Quadrants[i]);
                }
                else {
                    newele.push_back(Quadrants[i]);
                }
            }

            /*bool onein = false;
             vector<unsigned int> epts = Quadrants[i].getPoints();

             for (unsigned int j=0; j< epts.size(); j++) {
             if (points.at(epts[j]).isInside()) {// &&
             //!points.at(epts[j]).wasProjected()) {
             onein = true;
             break;
             }
             }

             if (onein) {
             newele.push_back(Quadrants[i]);
             }
             else {
             removed.push_back(Quadrants[i]);
             }*/
        }

        if (removed.empty()) {
            return;
        }

        //clear removed elements
        removed.clear();
        //now element std::list from Vomule mesh can be cleared, as all remaining
        //elements are still in use and attached to newele std::list.
        Quadrants.clear();
        for (auto eiter = newele.begin(); eiter!=newele.end(); ++eiter) {
            Quadrants.push_back(*eiter);
        }
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------

    void Mesher::applySurfacePatterns(Polyline &input){
        
        //apply patterns to avoid flat, invalid and
        //poor quality elements.
        SurfaceTemplatesVisitor stv;
        stv.setPoints(points);
        stv.setInput(input);

        for (auto q:Quadrants) {
            
            if (q.getPointIndex().size()!=4) {
                continue;
            }
            
            bool feature = false;
            for (auto octNo:q.getPointIndex()) {
                if (points[octNo].isFeature()) {
                    feature = true;
                    break;
                }
            }
            if (!feature) {
                if (q.isSurface()) {
                    if (!q.accept(&stv)) {
                        cout << "Error in Mesher::applySurfacePatterns: coultd't apply";
                        cout << " a surface pattern\n";
                        cout << q << "\n";
                        continue;
                    }
                }
            }
            else {
                unsigned int inNod = 0, oneN=0, oneF=0;
                
                vector<unsigned int> octIndx = q.getPointIndex();
                
                cout << "Feature Octant: ";
                for (unsigned int j=0; j<4; j++) {
                    if (points[octIndx[j]].isInside()) {
                        cout << "I";
                        inNod++;
                        oneN = j;
                    }
                    else {
                        cout << "O";
                    }
                    if (points[octIndx[j]].isFeature()) {
                        cout << "F";
                        oneF = j;
                    }
                    cout << " ";
                }
                if (inNod==1) {
                    cout << "1I";
                    if (!points[octIndx[(oneN+2)%4]].isFeature()) {
                        cout << " surface Pat";
                        q.accept(&stv);
                    }
                }
                else if (inNod==3) {
                    if (q.badAngle(oneF,points)) {
                        cout << " surface 3In";
                        q.accept(&stv);
                    }
                }
                
                cout << "\n";
            }
        }
    }


    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    //shrink elements intersecting the envelope defined by all
    //input surfaces

    void Mesher::shrinkToBoundary(Polyline &input){

        //Slow element removed (but works): from elements intersecting the
        //input domain, detect inner nodes. Project this nodes onto the
        //surface. If after all is done, if an element counts only with "on
        //surface" and "outside" nodes, remove it.
        list<unsigned int> out_nodes;
        list<Quadrant>::iterator oiter;
        
        
        //Manage Quadrants with Features first
        for (auto q:Quadrants) {
            if (!q.hasFeature()) {
                continue;
            }
            
            list<Point3D> fs = input.getFeatureProjection(q,points);
            
            if (fs.empty()) {
                cout << "wtf!!\n";
                continue;
            }
            
            const vector<unsigned int> &epts = q.getPointIndex();
            
            unsigned int fsNum = fs.size(), outNo = 0;
            for (auto pIdx:epts) {
                if (points[pIdx].isOutside()) {
                    outNo++;
                }
            }
            
            list<Point3D>::const_iterator iter;
            
            for (iter=fs.begin(); iter!=fs.end(); ++iter) {
                
                if (outNo==0) {
                    break;
                }
                
                double best = std::numeric_limits<double>::infinity();
                Point3D projected = *iter;
                unsigned int pos = 0;
                bool push = false;
                
                for (auto pIdx:epts) {
                    if (points[pIdx].isOutside()) {
                        
                        const Point3D &current = points[pIdx].getPoint();
                        double dis = (current - projected).Norm();
                        
                        if(best>dis){
                            best = dis;
                            pos = pIdx;
                            push = true;
                        }
                    }
                }
                    
                if (push) {
                    points[pos].setProjected();
                    points[pos].setPoint(projected);
                    points[pos].featureProjected();
                    for (auto pe:points[pos].getElements()) {
                        
                        //this should be studied further.
                        if (Quadrants.at(pe).intersectsSurface()) {
                            Quadrants[pe].setSurface();
                        }
                    }
                    outNo--;
                }
            }
        }
    
        //Manage non Feature Quadrants.
        for (auto q:Quadrants) {
            if (q.isInside()) {
                continue;
            }

            //Put in a std::list inside nodes of boundary elements that
            //may be projected to the input domain.
            const vector<unsigned int> &epts = q.getPointIndex();
            for (auto pIdx:q.getPointIndex()) {
                
                if (points[pIdx].wasProjected()) {
                    continue;
                }

                if (!points[pIdx].wasOutsideChecked()) {
                    cout << "error!!! in Mesher::shrinkToBoundary\n";
                    cout << "  some nodes were not outside checked (";
                    cout << q.getPointIndex().size() << "): " << pIdx << "\n";

                    points[pIdx].outsideChecked();
                    Point3D oct_p = points[pIdx].getPoint();
                    if (input.pointIsInMesh(oct_p,q.getIntersectedEdges())) {
                        points[pIdx].setInside();
                    }
                }

                if (points[pIdx].isOutside()) {
                    out_nodes.push_back(pIdx);
                }
            }
        }

        out_nodes.sort();
        out_nodes.unique();

        for (auto p:out_nodes) {

            //get the faces of Quadrants sharing this node
            list<unsigned int> p_qInterEdges;
            
            if (points.at(p).wasProjected()) {
                continue;
            }

            for (auto pe:points[p].getElements()) { //elements containing p
                //append this to the list of edges
                p_qInterEdges.insert(p_qInterEdges.end(),
                                     Quadrants[pe].getIntersectedEdges().begin(),
                                     Quadrants[pe].getIntersectedEdges().end());
            }

            p_qInterEdges.sort();
            p_qInterEdges.unique();

            if (p_qInterEdges.empty()) {
                cout << "\nWarning at Mesher::shrinkToBoundary";
                cout << " no faces to project an outside node\n";
                cout << p << " n_els " << points.at(p).getElements().size() << ":";
                for (auto pe:points.at(p).getElements()) {
                    cout << " " << pe;
                }
                cout << "\n";
                continue;
            }

            const Point3D &current = points[p].getPoint();
            Point3D projected = input.getProjection(current,p_qInterEdges);

            points[p].setPoint(projected);
            points[p].setProjected();
        }
    }


    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------

    void Mesher::projectCloseToBoundaryNodes(Polyline &input){

        //Slow element removed (but works): from elements intersecting the
        //input domain, detect inner nodes. Project this nodes onto the
        //surface. If after all is done, if an element counts only with "on
        //surface" and "outside" nodes, remove it.
        list<unsigned int> in_nodes;

        for (auto q:Quadrants) {
            if (!q.isSurface()) {
                continue;
            }
            q.computeMaxDistance(points);
        }
        
        //Manage Quadrants with Features first
        for (auto q:Quadrants) {
            if (!q.hasFeature()) {
                continue;
            }
            
            list<Point3D> fs = input.getFeatureProjection(q,points);
            const vector<unsigned int> &epts = q.getPointIndex();
            double md = q.getMaxDistance();
            
            for (auto pIdx:epts) {
                points.at(pIdx).setMaxDistance(md);
            }
            
            unsigned int fsNum = fs.size();
            
            for (auto pIdx:epts) {
                if (fsNum==0) {
                    break;
                }
                if (points[pIdx].isInside()) {
                
                    const Point3D &current = points[pIdx].getPoint();
                    double best = std::numeric_limits<double>::infinity();
                    Point3D candidate;
                    
                    //we use a list and interator to erase the projected node
                    //from the list of features.
                    list<Point3D>::iterator iter, bIter;
                    bool push = false;
                    for (iter=fs.begin(); iter!=fs.end(); ++iter) {
                        Point3D projected = *iter;
                        double dis = (current - projected).Norm();
                        
                        if(points[pIdx].getMaxDistance()>dis && best>dis){
                            best = dis;
                            candidate = projected;
                            bIter = iter;
                            push = true;
                        }
                    }
                    
                    if (push) {
                        fs.erase(bIter);
                        points[pIdx].setProjected();
                        points[pIdx].setPoint(candidate);
                        points[pIdx].featureProjected();
                        for (auto pe:points[pIdx].getElements()) {
                            //this should be studied further.
                            if (Quadrants.at(pe).intersectsSurface()) {
                                Quadrants[pe].setSurface();
                            }
                        }
                        fsNum--;
                    }
                }
            }
        }

        //Manage non Feature Quadrants.
        for (unsigned int i=0; i<Quadrants.size(); i++) {
            if (Quadrants[i].isInside() || Quadrants[i].hasFeature()) {
                continue;
            }

            //Put in a std::list inside nodes of boundary elements that
            //may be projected to the input domain.
            const vector<unsigned int> &epts = Quadrants.at(i).getPointIndex();
            double md = Quadrants[i].getMaxDistance();
            
            for (unsigned int j=0; j < epts.size(); j++) {

                if (!points[epts[j]].wasOutsideChecked()) {
                    cerr << "Warning in Mesher::projectCloseToBoundaryNodes\n";
                    cerr << "  point wasn't outside checked >>> Managed\n";

                    points[epts[j]].outsideChecked();
                    Point3D oct_p = points.at(epts[j]).getPoint();

                    if (input.pointIsInMesh(oct_p,Quadrants[i].getIntersectedEdges())) {
                        points[epts[j]].setInside();
                    }
                }
                if (points[epts[j]].isInside()) {
                    in_nodes.push_back(epts[j]);
                    /*cerr << "warning!! in Mesher::projectCloseToBoundaryNodes\n";
                    cerr << "  hard coded test for Octree j>7 ???\n";
                    if (j>7) {
                        md*=0.5;
                    }*/
                    points[epts[j]].setMaxDistance(md);
                }
            }

        }

        in_nodes.sort();
        in_nodes.unique(); //need to be sorted to call this function

        //move (when possible) all inner points to surface
        for (auto p:in_nodes) {

            
            if (points[p].wasProjected()) {
                //projected due to a feature => already managed.
                continue;
            }
            
            //if this node is attached to a Quadrant which was split in
            //mixed-elements due to transition patterns, avoid the
            //displacement.

            //get the faces of Quadrants sharing this node
            list<unsigned int> p_qInterEdges;

            //elements containing p
            for (auto pe:points.at(p).getElements()) {
                //append this to the list of edges
                p_qInterEdges.insert(p_qInterEdges.end(),
                                     Quadrants[pe].getIntersectedEdges().begin(),
                                     Quadrants[pe].getIntersectedEdges().end());
            }

            p_qInterEdges.sort();
            p_qInterEdges.unique();

            const Point3D &current = points[p].getPoint();
            Point3D projected = input.getProjection(current,p_qInterEdges);
            double dis = (current - projected).Norm();

            if(dis<points[p].getMaxDistance()){
                //this node has been moved to boundary, thus every element
                //sharing this node must be set as a border element in order
                //to avoid topological problems. And if this node belonged
                //to a feature Quadrants, all the quadrants sharing the node
                //will labeled as Feature Quadrants too.
                //points.at(*piter).setOutside();
                points[p].setProjected();
                points[p].setPoint(projected);
                for (auto pe:points[p].getElements()) {
                    
                    //intersectsSurface() returns true when the Quad previously
                    //had non empty list of intersected edges.
                    if (Quadrants[pe].intersectsSurface()) {
                        //setSurface() will change the flags of "inside" wich
                        Quadrants[pe].setSurface();
                    }
                }
            }
        }
    }
}
