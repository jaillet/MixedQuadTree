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
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    Mesher::Mesher() {}

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    Mesher::~Mesher() {}

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------

    //Create a grid mesh regarding the Bounding Box of input mesh.
    //This will produce several cubes as roots of an octree structure.
    //Then split each initial element 4^rl times (where rl stands
    //for Refinement Level).
    std::shared_ptr<FEMesh> Mesher::refineMesh(Polyline &input, const unsigned short &rl,
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

        //link element and node info for code optimization.
        //detect Quadrants with features.
        detectFeatureQuadrants(input);
        linkElementsToNodes();
        detectInsideNodes(input);

        //Now that we have all the elements, we can save the Quadrant mesh.
        unsigned int nels = Quadrants.size();
        Services::WriteQuadtreeMesh(name,points,Quadrants,QuadEdges,nels,gt);

//FJAXAV        projectCloseToBoundaryNodes(input);
//FJAXAV        removeOnSurfaceSafe(input);
        projectCloseToBoundaryNodes(input);
        removeOnSurfaceSafe(input);

        //update element and node info.
        linkElementsToNodes();

        //shrink outside nodes to the input domain boundary
        shrinkToBoundary(input);

        //apply the surface Patterns
        //FJAXAV        applySurfacePatterns(input);
        applySurfacePatterns(input);
        //removeOnSurface(input);

        if (rotated) {
            for (unsigned int i=0; i<points.size(); i++) {
                gt.applyInverse(points[i].getPoint());
            }
        }

        //the almighty output mesh
        std::shared_ptr<FEMesh> mesh = std::make_shared<FEMesh> ();

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
    std::shared_ptr<FEMesh> Mesher::generateMesh(Polyline &input, const unsigned short &rl,
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

        //link element and node info for code optimization, also
        //detect Quadrants with features.
        detectFeatureQuadrants(input);
        linkElementsToNodes();
        detectInsideNodes(input);

        //Now that we have all the elements, we can save the Quadrant mesh.
        unsigned int nels = Quadrants.size();
        Services::WriteQuadtreeMesh(name,points,Quadrants,QuadEdges,nels,gt);
        
        //Debbuging
        {
            //save pure octree mesh
            std::shared_ptr<FEMesh> pure_octree=make_shared<FEMesh>();
            saveOutputMesh(pure_octree,points,Quadrants);
            string tmp_name = name + "_remSur";
            Services::WriteVTK(tmp_name,pure_octree);
        }

//FJAXAV        projectCloseToBoundaryNodes(input);
//FJAXAV        removeOnSurfaceSafe(input);
        projectCloseToBoundaryNodes(input);
        
        removeOnSurfaceSafe(input);

        //update element and node info.
        linkElementsToNodes();

        //shrink outside nodes to the input domain boundary
        shrinkToBoundary(input);
        
        //apply the surface Patterns
        applySurfacePatterns(input);
        //removeOnSurface(input);

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
        std::shared_ptr<FEMesh> mesh = std::make_shared<FEMesh>();

        //save the data of the mesh in its final state
        saveOutputMesh(mesh,decoration);
        
        
        //Debugging:
        /*vector<unsigned int> color (Quadrants.size(),0);*/
        /*for (auto q:Quadrants) {
            if (q.isDebugging()) {
                cout << "For Quadrant " << q << " elements (" << q.getSubElements().size() << ")\n";
                for (auto se:q.getSubElements()) {
                    for (auto nIdx:se) {
                        cout << " " << nIdx;
                    }
                    cout << "\n";
                }
            }
        }*/

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

        auto start_time = chrono::high_resolution_clock::now();

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

        auto end_time = chrono::high_resolution_clock::now();
        cout << "    * generateGridMesh in "
             << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count();
        cout << " ms"<< endl;
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------

    void Mesher::splitQuadrants(const unsigned short &rl, Polyline &input,
                                list<unsigned int> &roctli,
                                list<RefinementRegion *> &all_reg, const string &name,
                                const unsigned short &minrl, const unsigned short &omaxrl){
        auto start_time = chrono::high_resolution_clock::now();

        //list of temp Quadrants
        list<Quadrant> tmp_Quadrants, new_Quadrants;
        //list of the points added at this refinement iteration:
        list<Point3D> new_pts;

        list<Quadrant>::iterator iter;

        //initialising, moving quadrants to save memory
        tmp_Quadrants.assign(make_move_iterator(Quadrants.begin()),make_move_iterator(Quadrants.end()));
        Quadrants.clear();

        //create visitors and give them variables
        SplitVisitor sv;
        sv.setPoints(points);
        sv.setEdges(QuadEdges);
        sv.setNewPts(new_pts);

        //----------------------------------------------------------
        //refine once each Quadrant in the list
        //----------------------------------------------------------
        unsigned int cindex = 0;

        auto start_refine_once_quad_time = chrono::high_resolution_clock::now();

        if (!roctli.empty()) {
            list<unsigned int>::iterator octidx = roctli.begin();
            
            while (!tmp_Quadrants.empty()) {
                iter=tmp_Quadrants.begin();
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
                // remove yet processed Quad
                tmp_Quadrants.pop_front();
                cindex++;
            }
            
            // dont' forget to update list
            std::swap(tmp_Quadrants,new_Quadrants);
            
            //add the new points to the vector
            points.reserve(points.size() + new_pts.size());
            points.insert(points.end(),new_pts.begin(),new_pts.end());
        }

        auto end_refine_once_quad_time = chrono::high_resolution_clock::now();
        cout << "       * Refine Once Quad in "
             << std::chrono::duration_cast<chrono::milliseconds>(end_refine_once_quad_time-start_refine_once_quad_time).count();
        cout << " ms"<< endl;

        //----------------------------------------------------------
        //refine each Quadrant until the Refinement Level is reached
        //----------------------------------------------------------
        auto start_refine_quad_time = chrono::high_resolution_clock::now();

        for (unsigned short i=minrl; i<rl; i++) {
            auto start_refine_rl_time = chrono::high_resolution_clock::now();

            bool changed = false;
            unsigned int refoct_int = 0, refoct_ext = 0;

            //the new_pts is a list that holds the coordinates of
            //new points inserted at this iteration. At the end of
            //this bucle, they are inserted in the point vector
            new_pts.clear();

            list<RefinementRegion *>::iterator reg_iter;

            //split the Quadrants as needed
            while(!tmp_Quadrants.empty()) {
                iter=tmp_Quadrants.begin();
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
                    // remove yet processed quad
                    tmp_Quadrants.pop_front();
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
                // remove yet processed quad
                tmp_Quadrants.pop_front();
            }
            
            // don't forget to update list
            std::swap(tmp_Quadrants,new_Quadrants);

            if (!changed) {
                continue;
            }

            //add the new points to the vector
            points.reserve(points.size() + new_pts.size());
            points.insert(points.end(),new_pts.begin(),new_pts.end());

            auto end_refine_rl_time = chrono::high_resolution_clock::now();
            cout << "         * level " << i << " in "
                 << std::chrono::duration_cast<chrono::milliseconds>(end_refine_rl_time-start_refine_rl_time).count();
            cout << " ms"<< endl;
        }

        auto end_refine_quad_time = chrono::high_resolution_clock::now();
        cout << "       * Refine Quad in "
             << std::chrono::duration_cast<chrono::milliseconds>(end_refine_quad_time-start_refine_quad_time).count();
        cout << " ms"<< endl;

        //----------------------------------------------------------
        //produce a one-irregular mesh
        //----------------------------------------------------------

        bool one_irregular = false;
//        new_Quadrants.clear();

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
            while(!tmp_Quadrants.empty()) {
                iter=tmp_Quadrants.begin();
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
                // remove yet processed quad
                tmp_Quadrants.pop_front();
            }

            // don't forget to update list
            std::swap(tmp_Quadrants,new_Quadrants);

            if (one_irregular) {
                break;
            }

            //if no points were added at this iteration, it is no longer
            //necessary to continue the refinement.
            if (new_pts.empty()) {
                break;
            }

            //add the new points to the vector
            points.reserve(points.size() + new_pts.size());
            points.insert(points.end(),new_pts.begin(),new_pts.end());
        }

        auto end_balanced_time = chrono::high_resolution_clock::now();
        cout << "       * Balanced mesh in "
             << std::chrono::duration_cast<chrono::milliseconds>(end_balanced_time-end_refine_quad_time).count();
        cout << " ms"<< endl;

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
                std::cerr << "Error at Mesher::generateQuadtreeMesh";
                std::cerr << " Transition Pattern not found\n";
            }
        }

        //if no points were added at this iteration, it is no longer
        //necessary to continue the refinement.
        if (!new_pts.empty()) {
            //add the new points to the vector
            points.reserve(points.size() + new_pts.size());
            points.insert(points.end(),new_pts.begin(),new_pts.end());
        }

        /*{
            //save pure octree mesh
            FEMesh pure_octree;
            saveOutputMesh(pure_octree,points,tmp_Quadrants);
            string tmp_name = name + "_alloct";
            Services::WriteVTK(tmp_name,pure_octree);
            Services::WriteMixedVolumeMesh(tmp_name,pure_octree);
        }*/

        // put (move!) the Quadrants in a vector
        Quadrants.clear(); //insert will reserve space as well
        Quadrants.insert(Quadrants.end(),make_move_iterator(tmp_Quadrants.begin()),make_move_iterator(tmp_Quadrants.end()));
        tmp_Quadrants.erase(tmp_Quadrants.begin(),tmp_Quadrants.end()); // better to erase as let in a indeterminate state by move

        auto end_time = chrono::high_resolution_clock::now();
        cout << "       * Transition Patterns in "
             << std::chrono::duration_cast<chrono::milliseconds>(end_time-end_balanced_time).count();
        cout << " ms"<< endl;
        cout << "    * splitQuadrants in "
             << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count();
        cout << " ms"<< endl;
    }


    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------

    void Mesher::generateQuadtreeMesh(const unsigned short &rl, Polyline &input,
                                    const list<RefinementRegion *> &all_reg,
                                    const string &name){

        auto start_time = chrono::high_resolution_clock::now();

        //to save m3d files per stage
        Clobscode::Services io;

        //list of temp Quadrants
        list<Quadrant> tmp_Quadrants, new_Quadrants;
        //list of the points added at this refinement iteration:
        list<Point3D> new_pts;
        list<Quadrant>::iterator iter;

        //initialize list with vector, FJA: sure we want a list?
        tmp_Quadrants.assign(Quadrants.begin(),Quadrants.end());

        //create visitors and give them variables
        SplitVisitor sv;
        sv.setPoints(points);
        sv.setEdges(QuadEdges);
        sv.setNewPts(new_pts);

        //----------------------------------------------------------
        //refine each Quadrant until the Refinement Level is reached
        //----------------------------------------------------------

        auto start_refine_quad_time = chrono::high_resolution_clock::now();

        for (unsigned short i=0; i<rl; i++) {
            auto start_refine_rl_time = chrono::high_resolution_clock::now();

            //the new_pts is a list that holds the coordinates of
            //new points inserted at this iteration. At the end of
            //this bucle, they are inserted in the point vector
            new_pts.clear();

            list<RefinementRegion *>::const_iterator reg_iter;

            //split the Quadrants as needed
            while (!tmp_Quadrants.empty()) {
                iter=tmp_Quadrants.begin();

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
                    list<unsigned int> &inter_edges = iter->getIntersectedEdges();

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
                // remove yet processed Quad
                tmp_Quadrants.pop_front();
            } // while

            // don't forget to update list
            std::swap(tmp_Quadrants,new_Quadrants);

            //if no points were added at this iteration, it is no longer
            //necessary to continue the refinement.
            if (new_pts.empty()) {
                cout << "warning at Mesher::generateQuadtreeMesh no new points!!!\n";
                break;
            }
            //add the new points to the vector
            points.reserve(points.size() + new_pts.size());
            points.insert(points.end(),new_pts.begin(),new_pts.end());

            auto end_refine_rl_time = chrono::high_resolution_clock::now();
            cout << "         * level " << i << " in "
                 << std::chrono::duration_cast<chrono::milliseconds>(end_refine_rl_time-start_refine_rl_time).count();
            cout << " ms"<< endl;
        }

        auto end_refine_quad_time = chrono::high_resolution_clock::now();
        cout << "       * Refine Quad in "
             << std::chrono::duration_cast<chrono::milliseconds>(end_refine_quad_time-start_refine_quad_time).count();
        cout << " ms"<< endl;

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
            while (!tmp_Quadrants.empty()) {
                iter=tmp_Quadrants.begin();
                if (!(*iter).accept(&oiv)){
                    //split this Quadrant
                    vector<vector<Point3D> > clipping_coords;
                    sv.setClipping(clipping_coords);
                    vector< vector <unsigned int> > split_elements;
                    sv.setNewEles(split_elements);

                    list<unsigned int> inter_faces = (*iter).getIntersectedEdges();

                    //(*iter).split(points,new_pts,QuadEdges,split_elements,clipping_coords);
                    //split the Quadrant
                    (*iter).accept(&sv);

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
                // remove yet processed Quad
                tmp_Quadrants.pop_front();

            } //while

            // don't forget to update list
            std::swap(tmp_Quadrants,new_Quadrants);

            if (one_irregular) {
                break;
            }

            //if no points were added at this iteration, it is no longer
            //necessary to continue the refinement.
            if (new_pts.empty()) {
                break;
            }
            //add the new points to the vector
            points.reserve(points.size() + new_pts.size());
            points.insert(points.end(),new_pts.begin(),new_pts.end());
        }

        auto end_balanced_time = chrono::high_resolution_clock::now();
        cout << "       * Balanced mesh in "
             << std::chrono::duration_cast<chrono::milliseconds>(end_balanced_time-end_refine_quad_time).count();
        cout << " ms"<< endl;

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
                std::cerr << "Error at Mesher::generateQuadtreeMesh";
                std::cerr << " Transition Pattern not found\n";
            }
        }

        
        
        //Debbuging
        /*{
            //save pure octree mesh
            FEMesh pure_octree;
            saveOutputMesh(pure_octree,points,tmp_Quadrants);
            string tmp_name = name + "_oct";
            Services::WriteVTK(tmp_name,pure_octree);
            //Services::WriteMixedVolumeMesh(tmp_name,pure_octree);
        }*/
        
        
        
        //if no points were added at this iteration, it is no longer
        //necessary to continue the refinement.
        if (!new_pts.empty()) {
            //add the new points to the vector
            points.reserve(points.size()+new_pts.size());
            points.insert(points.end(),new_pts.begin(),new_pts.end());
        }

        // put (move!) the Quadrants into a vector
        //        Quadrants.assign(tmp_Quadrants.begin(),tmp_Quadrants.end());
        Quadrants.clear(); //insert will reserve space as well
        Quadrants.insert(Quadrants.end(),make_move_iterator(tmp_Quadrants.begin()),make_move_iterator(tmp_Quadrants.end()));
        tmp_Quadrants.erase(tmp_Quadrants.begin(),tmp_Quadrants.end()); // better to erase as let in a indeterminate state by move

        auto end_time = chrono::high_resolution_clock::now();
        cout << "       * Transition Patterns in "
             << std::chrono::duration_cast<chrono::milliseconds>(end_time-end_balanced_time).count();
        cout << " ms"<< endl;
        cout << "    * generateQuadtreeMesh in "
             << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count() ;
        cout << " ms"<< endl;

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

    unsigned int Mesher::saveOutputMesh(const shared_ptr<FEMesh> &mesh, bool decoration){
        auto start_time = chrono::high_resolution_clock::now();

        vector<vector<unsigned int> > out_els;
        vector<unsigned short > out_els_ref_level, out_els_surf;
        vector<double > out_els_min_angle;
        
        //even if we don't know the quantity of elements, it will be at least the
        //number of quadrants, so we reserve for the vectors of VTK decoration.
        out_els_min_angle.reserve(Quadrants.size());
        out_els_ref_level.reserve(Quadrants.size());
        out_els_surf.reserve(Quadrants.size());

        //new_idxs will hold the index of used nodes in the outside vector for points.
        //If the a node is not used by any element, its index will be 0 in this vector,
        //therefore the actual index is shiffted in 1. In other words, node 0 is node 1,
        //and node n is node n+1.
        vector<unsigned int> new_idxs (points.size(),0);
        unsigned int out_node_count = 0;
        vector<Point3D> out_pts;
        
        /*Begin debugging:*/
//        for (auto q:Quadrants) {
//            for (auto el:q.getSubElements()) {
//                if (decoration) {
//                    if (q.isDebugging()) {
//                        out_els_ref_level.push_back(1);
//                    }
//                    else {
//                        out_els_ref_level.push_back(3);
//                    }
//                }
//                out_els.push_back(el);
//            }
//        }
//        for (auto p:points) {
//            out_pts.push_back(p.getPoint());
//        }
        
//        mesh.setPoints(out_pts);
//        mesh.setElements(out_els);
//        mesh.setRefLevels(out_els_ref_level);
        
//        return out_els.size();
        /*End debugging:*/
        
        

        //recompute node indexes and update elements with them.
        for (unsigned int i=0; i<Quadrants.size(); i++) {
            const vector<vector<unsigned int> > &sub_els= Quadrants[i].getSubElements();
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
                    //surface regardring quad
                    if (Quadrants[i].isSurface()) {
                        out_els_surf.push_back(1);
                    }
                    else {
                        out_els_surf.push_back(0);
                    }
                    
                    //refinment level herited from quad
                    out_els_ref_level.push_back(Quadrants[i].getRefinementLevel());
                    //compute minAngle
                    unsigned int np=sub_ele_new_idxs.size(); //nb points of the element
                    double minAngle=std::numeric_limits<double>::infinity();
                    for (unsigned int k=0; k<np; ++k) {

                        const Point3D &P0 = out_pts[sub_ele_new_idxs[(k-1+np)%np]];
                        const Point3D &P1 = out_pts[sub_ele_new_idxs[k]];
                        const Point3D &P2 = out_pts[sub_ele_new_idxs[(k+1)%np]];

                        minAngle=std::min(minAngle, P1.angle3Points(P0,P2));
                    }
                    out_els_min_angle.push_back(minAngle);
                }
                out_els.push_back(sub_ele_new_idxs);
            }
        }

        mesh->setPoints(out_pts);
        mesh->setElements(out_els);
        mesh->setRefLevels(out_els_ref_level);
        mesh->setMinAngles(out_els_min_angle);
        mesh->setSurfState(out_els_surf);

        auto end_time = chrono::high_resolution_clock::now();
        cout << "    * SaveOutputMesh in "
             << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count();
        cout << " ms"<< endl;


        return out_els.size();
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------

    unsigned int Mesher::saveOutputMesh(const shared_ptr<FEMesh> &mesh, vector<MeshPoint> &tmp_points,
                                list<Quadrant> &tmp_Quadrants){
        auto start_time = chrono::high_resolution_clock::now();

        vector<Point3D> out_pts;
        list<vector<unsigned int> > tmp_elements;
        vector<vector<unsigned int> > out_els;

        unsigned int n = tmp_points.size();
        out_pts.reserve(n);
        for (unsigned int i=0; i<n; i++) {
            out_pts.push_back(points[i].getPoint());
        }

        list<Quadrant>::const_iterator o_iter;

        for (o_iter=tmp_Quadrants.begin(); o_iter!=tmp_Quadrants.end(); ++o_iter) {

            vector<vector<unsigned int> > sub_els= o_iter->getSubElements();
            for (unsigned int j=0; j<sub_els.size(); j++) {
                tmp_elements.push_back(sub_els[j]);
            }
        }

        out_els.reserve(tmp_elements.size());
        list<vector<unsigned int> >::const_iterator e_iter;

        for (e_iter=tmp_elements.begin(); e_iter!=tmp_elements.end(); ++e_iter) {
            out_els.push_back(*e_iter);
        }

        mesh->setPoints(out_pts);
        mesh->setElements(out_els);

        auto end_time = chrono::high_resolution_clock::now();
        cout << "    * SaveOutputMesh in "
             << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count();
        cout << " ms"<< endl;

        return out_els.size();
    }
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------

    unsigned int Mesher::saveOutputMesh(const shared_ptr<FEMesh> &mesh, const vector<MeshPoint> &tmp_points,
                                        const vector<Quadrant> &tmp_Quadrants){
        auto start_time = chrono::high_resolution_clock::now();

        vector<Point3D> out_pts;
        list<vector<unsigned int> > tmp_elements;
        vector<vector<unsigned int> > out_els;

        unsigned int n = tmp_points.size();
        out_pts.reserve(n);
        for (unsigned int i=0; i<n; i++) {
            out_pts.push_back(points[i].getPoint());
        }

        //list<Quadrant>::const_iterator o_iter;

        for (unsigned int i=0; i<tmp_Quadrants.size(); i++) {
            vector<vector<unsigned int> > sub_els= tmp_Quadrants[i].getSubElements();
            for (unsigned int j=0; j<sub_els.size(); j++) {
                tmp_elements.push_back(sub_els[j]);
            }
        }

        /*for (o_iter=tmp_Quadrants.begin(); o_iter!=tmp_Quadrants.end(); ++o_iter) {

            vector<vector<unsigned int> > sub_els= o_iter->getSubElements();
            for (unsigned int j=0; j<sub_els.size(); j++) {
                tmp_elements.push_back(sub_els[j]);
            }
        }*/

        out_els.reserve(tmp_elements.size());
        list<vector<unsigned int> >::const_iterator e_iter;

        for (e_iter=tmp_elements.begin(); e_iter!=tmp_elements.end(); ++e_iter) {
            out_els.push_back(*e_iter);
        }

        mesh->setPoints(out_pts);
        mesh->setElements(out_els);

        auto end_time = chrono::high_resolution_clock::now();
        cout << "    * SaveOutputMesh in "
        << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count();
        cout << " ms"<< endl;

        return out_els.size();
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    void Mesher::detectFeatureQuadrants(Polyline &input) {
        
        auto start_time = chrono::high_resolution_clock::now();

        unsigned int featCount = 0;
        for (unsigned int i=0; i<Quadrants.size(); i++) {
            if (input.getNbFeatures(Quadrants[i],points)>0) {
                //Quadrants[i].setFeature(); now set by Polyline::getNbFeatures()
                featCount++;
            }
        }
        //cout << "number of quadrants with features: " << featCount << "\n";

        auto end_time = chrono::high_resolution_clock::now();
        cout << "    * detectFeatureQuadrants in "
             << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count();
        cout << " ms"<< endl;

    }
    
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    void Mesher::linkElementsToNodes(){

        auto start_time = chrono::high_resolution_clock::now();

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
        auto end_time = chrono::high_resolution_clock::now();
        cout << "    * linkElementsToNodes in "
             << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count();
        cout << " ms"<< endl;
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------

    void Mesher::detectInsideNodes(Polyline &input){
        auto start_time = chrono::high_resolution_clock::now();

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
        auto end_time = chrono::high_resolution_clock::now();
        cout << "    * detectInsideNodes in "
             << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count();
        cout << " ms"<< endl;
    }
    
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    
    void Mesher::removeOnSurface(Polyline &input){
        auto start_time = chrono::high_resolution_clock::now();

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
        }
        
        if (removed.empty()) {
            return;
        }
        
        //clear removed elements
        removed.clear();
        //now element std::list from Vomule mesh can be cleared, as all remaining
        //elements are still in use and attached to newele std::list.
        Quadrants.clear();
        Quadrants.assign(make_move_iterator(newele.begin()),make_move_iterator(newele.end()));

        auto end_time = chrono::high_resolution_clock::now();
        cout << "    * RemoveOnSurface in "
             << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count();
        cout << " ms"<< endl;
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------

    void Mesher::removeOnSurfaceSafe(Polyline &input){
        auto start_time = chrono::high_resolution_clock::now();

        list<Quadrant> newele,removed;

        //Keep all inside quads and all of th
        for (auto q:Quadrants) {
            if (q.isInside()) {
                newele.push_back(q);
                continue;
            }

            bool onein = false;
            bool feaProj = false;
            Point3D avg;
            for (auto quaNoIdx:q.getPointIndex()) {
                avg+=points[quaNoIdx].getPoint();
                if (points[quaNoIdx].isInside()) {
                    onein = true;
                    break;
                }
                if (points[quaNoIdx].isFeature()) {
                    feaProj = true;
                }
            }
            
            //if it has at least one node inside
            //we don't remove it.
            if (onein) {
                newele.push_back(q);
                continue;
            }
            avg/=q.getPointIndex().size();
            //now we now the Quad has no node inside.
            if (input.pointIsInMesh(avg,q.getIntersectedEdges())) {
                newele.push_back(q);
            }
            else {
                removed.push_back(q);
            }
            
            //the most expensive case: regarding
            //the original intersected edges, does
            //this quad still have a node inside?
            //if yes we keep it. It means that all
            //of its nodes are outside or projected
            //onto the boundary so if we erase it we
            //will loose boundary representation.
            /*bool border = false;
            for (auto qp:q.getPointIndex()) {
                if (input.pointIsInMesh(points[qp].getPoint(),q.getIntersectedEdges())) {
                    newele.push_back(q);
                    border = true;
                    break;
                }
            }
            if (!border) {
                removed.push_back(q);
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
        Quadrants.assign(std::make_move_iterator(newele.begin()),std::make_move_iterator(newele.end()));
        
        auto end_time = chrono::high_resolution_clock::now();
        cout << "    * RemoveOnSurfaceSafe in "
             << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count();
        cout << " ms"<< endl;
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------

    void Mesher::applySurfacePatterns(Polyline &input){
        auto start_time = chrono::high_resolution_clock::now();

        //apply patterns to avoid flat, invalid and
        //poor quality elements.
        SurfaceTemplatesVisitor stv;
        stv.setPoints(points);
        stv.setPolyline(input);


        for (auto &q:Quadrants) {
            
            if (q.getPointIndex().size()!=4) {
                continue;
            }

            if (q.isSurface()) {
                if (!q.accept(&stv)) {
                    cout << "Error in Mesher::applySurfacePatterns: coultd't apply";
                    cout << " a surface pattern\n";
                    cout << q << "\n";
                    continue;
                }
            }

        }
        
        
        auto end_time = chrono::high_resolution_clock::now();
        cout << "    * ApplySurfacePatterns in "
             << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count();
        cout << " ms"<< endl;
    }


    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    //shrink elements intersecting the envelope defined by all
    //input surfaces

    void Mesher::shrinkToBoundary(Polyline &input){
        auto start_time = chrono::high_resolution_clock::now();

        //Slow element removed (but works): from elements intersecting the
        //input domain, detect inner nodes. Project this nodes onto the
        //surface. If after all is done, if an element counts only with "on
        //surface" and "outside" nodes, remove it.
        list<unsigned int> out_nodes;
        list<Quadrant>::iterator oiter;
        
        
        //Manage Quadrants with Features first
        for (const auto &q:Quadrants) {
            if (q.isInside()) {
                continue;
            }
            
            //Save all outside nodes and manage them after
            //dealing with all features quadrants.
            //if a node was projected to a feature
            //it will be skipped during the rest of
            //node projection.
            for (auto pIdx:q.getPointIndex()) {
                if (points[pIdx].isOutside()) {
                    out_nodes.push_back(pIdx);
                }
            }
            
            /*
             All features were already managed.
             
             if (!q.hasFeature()) {
                continue;
            }
            
            list<unsigned int> fs = input.getFeatureProjection(q,points);
            
            if (fs.empty()) {
                cerr << "Error at Mesher::shrinkToBoundary";
                cerr << " Quadrant labeled with feature has none\n";
                continue;
            }
            
            const vector<unsigned int> &epts = q.getPointIndex();
            
            unsigned int fsNum = fs.size(), outNo = 0;
            for (auto pIdx:epts) {
                if (points[pIdx].isOutside() && !points[pIdx].wasProjected()) {
                    outNo++;
                }
            }
            
            list<unsigned int>::const_iterator iter;
            
            for (iter=fs.begin(); iter!=fs.end(); ++iter) {
                
                if (outNo==0) {
                    break;
                }
                
                double best = std::numeric_limits<double>::infinity();
                Point3D projected = input.getPoints()[*iter];
                unsigned int pos = 0;
                bool push = false;
                
                for (auto pIdx:epts) {
                    if (points[pIdx].isOutside() && !points[pIdx].wasProjected()) {
                        
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
                    //Feature projected flag will be used later to
                    //apply surface patterns.
                    points[pos].featureProjected();
                    for (auto pe:points[pos].getElements()) {
                        
                        //this should be studied further.
                        if (Quadrants.at(pe).intersectsSurface()) {
                            Quadrants[pe].setSurface();
                        }
                    }
                    outNo--;
                }
            }*/
        }

        //Manage non Feature Quadrants.
        out_nodes.sort();
        out_nodes.unique();

        for (auto p:out_nodes) {

            if (points.at(p).wasProjected()) {
                continue;
            }
            
            //get the faces of Quadrants sharing this node
            list<unsigned int> p_qInterEdges;

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
        auto end_time = chrono::high_resolution_clock::now();
        cout << "    * ShrinkToBoundary in "
             << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count();
        cout << " ms"<< endl;
    }


    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------

    void Mesher::projectCloseToBoundaryNodes(Polyline &input){
        auto start_time = chrono::high_resolution_clock::now();

        //Slow element removed (but works): from elements intersecting the
        //input domain, detect inner nodes. Project this nodes onto the
        //surface. If after all is done, if an element counts only with "on
        //surface" and "outside" nodes, remove it.
        list<unsigned int> in_nodes;

        //FJA compute everything before to project anything
        for (auto &q:Quadrants) {
            q.computeMaxDistance(points);
        }

        for (auto &q:Quadrants) {

            if (!q.isSurface()) {
                continue;
            }
            
            double md = q.getMaxDistance();
            for (auto pIdx:q.getPointIndex()) {
                points[pIdx].setMaxDistance(md);
            }
        
            //If the quadrant doesn't have a Feature, save the intern nodes
            //to a list to manage them later and continue the process just
            //for feature quadrants first.
            if (!q.hasFeature()) {
                for (auto pIdx:q.getPointIndex()) {
                    if (points[pIdx].isInside()) {
                        in_nodes.push_back(pIdx);
                    }
                }
                continue;
            }
        
            //Manage Quadrants with Features
            list<unsigned int> fs = input.getFeatureProjection(q,points);
            
            unsigned int fsNum = fs.size();
            //we use a list and interator to erase the projected node
            //from the list of features.
            list<unsigned int>::const_iterator itFs;
            for (itFs=fs.begin(); itFs!=fs.end(); ++itFs) {

                if (fsNum==0) {
                    break;
                }

                Point3D featProjected = input.getPoints()[*itFs];

                double best = std::numeric_limits<double>::infinity();
                unsigned int candidate=-1;
                
                bool push = false;
                for (auto pIdx:q.getPointIndex()) {
                    
                    if (points[pIdx].wasProjected()) {
                        continue;
                    }
                    
                    const Point3D &current = points[pIdx].getPoint();
                    double dis = (current - featProjected).Norm();

                    //if(points[pIdx].getMaxDistance()>dis && best>dis){
                    if(best>dis) {
                        best = dis;
                        candidate = pIdx;
                        push = true;
                    }
                }
                
                if (push) {
                    points[candidate].setProjected();
                    points[candidate].setPoint(featProjected);
                    //Feature projected flag will be used later to
                    //apply surface patterns.
                    points[candidate].featureProjected();

                    // First, get index of the feature edges. Will be added next
                    // to intersecting edges of the quads sharing the candidate node
                    // find first edge
//                    vector<PolyEdge>::iterator it=std::find_if(input.getEdges().begin(), input.getEdges().end(),
//                                      [itFs] (const PolyEdge& pe) -> bool {return ((pe.getKey() == *itFs)||(pe.getVal()== *itFs)) ; });
//                    unsigned ie1=std::distance(input.getEdges().begin(),it);
//                    // find second edge
//                    it = std::find_if(++it, input.getEdges().end(),
//                                      [itFs] (const PolyEdge& pe) -> bool {return ((pe.getKey() == *itFs)||(pe.getVal()== *itFs)) ; });
//                    unsigned ie2=std::distance(input.getEdges().begin(),it);

                    for (auto pe:points[candidate].getElements()) {
                        //this should be studied further.
                        if (Quadrants.at(pe).intersectsSurface()) {
                            Quadrants[pe].setSurface();
                        }
//                        // add feature edges to the intersected list
//                        Quadrants[pe].getIntersectedEdges().push_back(ie1);
//                        Quadrants[pe].getIntersectedEdges().push_back(ie2);
//                        // remove duplicates
//                        Quadrants[pe].getIntersectedEdges().sort();
//                        Quadrants[pe].getIntersectedEdges().unique();
                    }

                    fsNum--;
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
                //to avoid topological problems.
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
        auto end_time = chrono::high_resolution_clock::now();
        cout << "    * ProjectCloseToBoundary in "
             << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count();
        cout << " ms"<< endl;

    }
}
