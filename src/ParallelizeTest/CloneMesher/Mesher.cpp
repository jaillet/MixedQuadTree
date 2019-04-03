//Fastest way to have big data ?

#include "../../Point3D.h"
#include "../../Quadrant.h"
#include <list>
#include <stdio>
#include <omp.h>


//Idea :

// As edges are read/write for each quadrant, and we cannot know how openMP
//Slice the for loop, maybe manually slice the loop ?

/**
 * @brief Implementation of algorithm 2 of mesher
 * @details Uses OpenMP.
 * 
 * Output of this function :
 * 
 * tmp_Quadrants (or new_quadrants before the swap) -> Quadrants that are created
 * 
 * 
 * @param nbThread number of thread, must be > 0 and < MAX
 */
void refineMeshParallel(int nbThread) {
	//OpenMP
    if (nbThread > NOMBRE_THREAD || nbThread < 0) {
        std::cout << "Invalid number of threads or not supported by computer" << endl;
        return;
    }

    omp_set_num_threads(NOMBRE_THREAD);


	//list of temp Quadrants
    std::vector<Quadrant> tmp_Quadrants, new_Quadrants;
    //list of the points added at this refinement iteration:
    list<Point3D> new_pts;

    list<Quadrant>::iterator iter;

    //initialising into a list, moving quadrants to save memory
    tmp_Quadrants.assign(make_move_iterator(Quadrants.begin()), make_move_iterator(Quadrants.end()));
    Quadrants.clear();

    //create visitors and give them variables
    SplitVisitor sv;
    sv.setPoints(points);
    sv.setEdges(QuadEdges);
    sv.setNewPts(new_pts);


    //Shared variables

    //points
   	//--- Read in SplitVisitor::visit
    //QuadEdges Warning!
    //--- is filled by SplitVisitor::visit
    //--- and needed and filled in SplitVisitor::splitEdge!

	int rl = 9; //Refinement level
	for (unsigned short i = 0; i < rl; i++) {
        auto start_refine_rl_time = chrono::high_resolution_clock::now();

        //the new_pts is a list that holds the coordinates of
        //new points inserted at this iteration. At the end of
        //this bucle, they are inserted in the point vector
        new_pts.clear();

        list<RefinementRegion *>::const_iterator reg_iter;

        //split the Quadrants as needed
        //begin parallelisation ?
        //new_pts is shared, and new_quadrants
        #pragma parallel for shared(new_points, new_Quadrants) first_private(reg_iter)
        for (int j = 0; j < tmp_Quadrants.size(); ++j) {
   			
            Quadrant* iter = &tmp_Quadrants[j];

            //iter = tmp_Quadrants.begin();


            //Only check, can not modify after treatment
            bool to_refine = false;
            for (reg_iter = all_reg.begin(), reg_iter++; reg_iter != all_reg.end(); ++reg_iter) {

                unsigned short region_rl = (*reg_iter)->getRefinementLevel();
                if (region_rl < i) {
                    continue;
                }

                //If the Quadrant has a greater RL than the region needs, continue
                if (region_rl <= (*iter).getRefinementLevel()) {
                    continue;
                }

                //Get the two extreme nodes of the Quadrant to test intersection with
                //this RefinementRegion. If not, conserve it as it is.
                //unsigned int n_idx1 = (*iter).getPoints()[0];
                //unsigned int n_idx2 = (*iter).getPoints()[2];


                // intersectQaudrant can modify the quadrant with 
                // function Polyline::getNbFeatures in RefinementboundaryRegion
                // maybe no problem for parallelisation, as this information is
                // not used by other thread
                if ((*reg_iter)->intersectsQuadrant(points, *iter)) {
                    to_refine = true;
                }
            }

            //now if refinement is not needed, we add the Quadrant as it was.
            if (!to_refine) {
                new_Quadrants.push_back(*iter); //SHARED VARIABLE
                //End of for loop
            } else {
                list<unsigned int> &inter_edges = iter->getIntersectedEdges();
                unsigned short qrl = (*iter).getRefinementLevel();

                vector<vector<Point3D> > clipping_coords;
                sv.setClipping(clipping_coords);

                vector<vector<unsigned int> > split_elements;
                sv.setNewEles(split_elements);

                //iter->split(points,new_pts,QuadEdges,split_elements,clipping_coords);
                //cout << "Accept" << endl;
                iter->accept(&sv);

                if (inter_edges.empty()) {
                    for (unsigned int j = 0; j < split_elements.size(); j++) {
                        Quadrant o(split_elements[j], qrl + 1);
                        new_Quadrants.push_back(o);
                    }
                } else {
                    for (unsigned int j = 0; j < split_elements.size(); j++) {
                        Quadrant o(split_elements[j], qrl + 1);
                        //the new points are inserted in bash at the end of this
                        //iteration. For this reason, the coordinates must be passed
                        //"manually" at this point (clipping_coords).

                        //select_faces = true
                        IntersectionsVisitor iv(true);
                        //if (o.checkIntersections(input,inter_edges,clipping_coords[j]))
                        iv.setPolyline(input);
                        iv.setEdges(inter_edges);
                        iv.setCoords(clipping_coords[j]);

                        if (o.accept(&iv)) {
                            new_Quadrants.push_back(o);
                        } else {
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

                            if (isItIn(input, inter_edges, clipping_coords[j])) {
                                new_Quadrants.push_back(o);
                            }
                        }
                    }
                }
            }
        } // for all quadrant end parallel region

        // don't forget to update list
        std::swap(tmp_Quadrants, new_Quadrants);

        //if no points were added at this iteration, it is no longer
        //necessary to continue the refinement.
        if (new_pts.empty()) {
            cout << "warning at Mesher::generateQuadtreeMesh no new points!!!\n";
            break;
        }

        //add the new points to the vector
        points.reserve(points.size() + new_pts.size());
        points.insert(points.end(), new_pts.begin(), new_pts.end());

    } //END FOR REFINEMENT LEVEL

    auto end_refine_quad_time = chrono::high_resolution_clock::now();
    cout << "       * Refine Quad in "
         << std::chrono::duration_cast<chrono::milliseconds>(end_refine_quad_time - start_refine_quad_time).count();
    cout << " ms" << endl;

}