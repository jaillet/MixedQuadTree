//Fastest way to have big data ?

#include "../../RefinementRegion.h"
#include "../../Visitors/IntersectionsVisitor.h"
#include "SplitVisitorTest1.h"

#include <list>
#include <chrono>
#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_for.h>
#include <tbb/concurrent_vector.h>
#include <tbb/concurrent_queue.h>
#include <atomic>
#include <mutex>

#include "omp.h"

using namespace Clobscode;


//Idea :

// As edges are read/write for each quadrant, and we cannot know how openMP
//Slice the for loop, maybe manually slice the loop ?

bool isItIn(const Polyline &mesh, const list<unsigned int> &faces, const vector<Point3D> &coords) {
    //this method is meant to be used by Quadrants that don't
    //intersect input domains. If they are inside of at least
    //one input mesh, then they must remain in the output mesh.

    bool first = mesh.pointIsInMesh(coords[0], faces);
    bool second = mesh.pointIsInMesh(coords[1], faces);
    if (first == second) {
        return first;
    }

    //cout << "one inconsistency detected -> hard test\n";
    //return mesh.pointIsInMesh(coords[0],faces);
    return mesh.pointIsInMesh(coords[0]);
}

void refineMeshParallelTest1TBB(int nbThread, list<Quadrant> Quadrants, vector<MeshPoint> points,
                                        set<QuadEdge> QuadEdges,
                                        const list<RefinementRegion *> &all_reg, const unsigned short &rl,
                                        Polyline &input) {

    int NOMBRE_THREAD = tbb::task_scheduler_init::default_num_threads();
    std::cout << NOMBRE_THREAD << std::endl;

    //OpenMP
    if (nbThread > NOMBRE_THREAD || nbThread < 0) {
        std::cout << "Invalid number of threads or not supported by computer" << std::endl;
        return;
    }

    tbb::task_scheduler_init test(nbThread);


    //list of temp Quadrants
    vector<Quadrant> tmp_Quadrants;
    tbb::concurrent_vector<Quadrant> new_Quadrants;

    //list of the points added at this refinement iteration:
    // TODO !! INSERT IN SPLITVISITOR !!
    list<Point3D> new_pts;

    //initialising into a list, moving quadrants to save memory
    tmp_Quadrants.assign(make_move_iterator(Quadrants.begin()), make_move_iterator(Quadrants.end()));
    Quadrants.clear();

    tbb::atomic<int> counterRefine = 0;

    //Shared variables

    //points
    //--- Read in SplitVisitor::visit
    //QuadEdges Warning!
    //--- is filled by SplitVisitor::visit
    //--- and needed and filled in SplitVisitor::splitEdge!


    tbb::atomic<int> nb_points;
    std::mutex mtx_new_pts;
    std::mutex mtx_new_edges;

    for (unsigned short i = 0; i < rl; i++) {
        auto start_refine_rl_time = chrono::high_resolution_clock::now();

        //the new_pts is a list that holds the coordinates of
        //new points inserted at this iteration. At the end of
        //this bucle, they are inserted in the point vector
        new_pts.clear();

        nb_points = points.size();

        //split the Quadrants as needed
        //begin parallelisation ?
        //new_pts is shared, and new_quadrants
        tbb::parallel_for(tbb::blocked_range<std::size_t>(0, tmp_Quadrants.size()),
                          [&](const tbb::blocked_range<std::size_t> &range) {

              //create visitors and give them variables
              SplitVisitorTest1 sv;
              sv.setPoints(points); // READ
              sv.setEdges(
                      QuadEdges); // TODO INSERT / REMOVE / READ A faire en priorite car le plus contraignant
              sv.setNewPts(new_pts); // TODO INSERT / READ
              sv.setCounterPointst(&nb_points);
              sv.setMutexForPoints(&mtx_new_pts);
              sv.setMutexForEdges(&mtx_new_edges);

              list<RefinementRegion *>::const_iterator reg_iter;

              for (auto j = range.begin(); j != range.end(); ++j) {

                  Quadrant &iter = tmp_Quadrants[j];

                  //Only check, can not modify after treatment
                  bool to_refine = false;

                  for (reg_iter = all_reg.begin(), reg_iter++;
                       reg_iter != all_reg.end(); ++reg_iter) {


                      unsigned short region_rl = (*reg_iter)->getRefinementLevel();
                      if (region_rl < i) {
                          continue;
                      }

                      //If the Quadrant has a greater RL than the region needs, continue
                      if (region_rl <= iter.getRefinementLevel()) {
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
                      if ((*reg_iter)->intersectsQuadrant(points, iter)) {
                          to_refine = true;
                          //counterRefine.fetch_and_increment();
                          break;
                      }
                  }


                  //now if refinement is not needed, we add the Quadrant as it was.
                  if (!to_refine) {
                      new_Quadrants.push_back(iter); //SHARED VARIABLE
                      //End of for loop
                  } else {

                      list<unsigned int> &inter_edges = iter.getIntersectedEdges();
                      unsigned short qrl = iter.getRefinementLevel();

                      vector<vector<Point3D> > clipping_coords;
                      sv.setClipping(clipping_coords);

                      vector<vector<unsigned int> > split_elements;
                      sv.setNewEles(split_elements);

                      iter.accept(&sv);

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
              } //END FOR QUADRANTS
          
          }); //END TBB TASK

        // don't forget to update list
        tmp_Quadrants.assign(make_move_iterator(new_Quadrants.begin()), make_move_iterator(new_Quadrants.end()));
        new_Quadrants.clear();
        //PAUL: Why not swap??

        //if no points were added at this iteration, it is no longer
        //necessary to continue the refinement.
        if (new_pts.empty()) {
            cout << "warning at Mesher::generateQuadtreeMesh no new points!!!\n";
            break;
        }

        //add the new points to the vector
        points.reserve(points.size() + new_pts.size());
        points.insert(points.end(), new_pts.begin(), new_pts.end());

        auto end_refine_rl_time = chrono::high_resolution_clock::now();
        cout << "         * level " << i << " in "
             << std::chrono::duration_cast<chrono::milliseconds>(end_refine_rl_time - start_refine_rl_time).count();
        cout << " ms" << endl;

    } //END FOR REFINEMENT LEVEL

    std::cout << counterRefine << std::endl;
}

void refineMeshParallelTest1(int nbThread, list<Quadrant> Quadrants, vector<MeshPoint> points,
                                     set<QuadEdge> QuadEdges,
                                     const list<RefinementRegion *> &all_reg, const unsigned short &rl,
                                     Polyline &input) {

    int NOMBRE_THREAD = omp_get_max_threads();

    //OpenMP
    if (nbThread > NOMBRE_THREAD || nbThread < 0) {
        std::cout << "Invalid number of threads or not supported by computer" << std::endl;
        return;
    }

    omp_set_num_threads(NOMBRE_THREAD);


    //list of temp Quadrants
    vector<Quadrant> tmp_Quadrants;
    tbb::concurrent_vector<Quadrant> new_Quadrants;

    //list of the points added at this refinement iteration:
    // TODO !! INSERT IN SPLITVISITOR !!
    list<Point3D> new_pts;

    //initialising into a list, moving quadrants to save memory
    tmp_Quadrants.assign(make_move_iterator(Quadrants.begin()), make_move_iterator(Quadrants.end()));
    Quadrants.clear();

    tbb::atomic<int> counterRefine = 0;

    //Shared variables

    //points
    //--- Read in SplitVisitor::visit
    //QuadEdges Warning!
    //--- is filled by SplitVisitor::visit
    //--- and needed and filled in SplitVisitor::splitEdge!


    tbb::atomic<int> nb_points;
    std::mutex mtx_new_pts;
    std::mutex mtx_new_edges;

    //create visitors and give them variables
    SplitVisitorTest1 sv;
    sv.setPoints(points); // READ
    sv.setEdges(QuadEdges); // TODO INSERT / REMOVE / READ A faire en priorite car le plus contraignant
    sv.setNewPts(new_pts); // TODO INSERT / READ
    sv.setCounterPointst(&nb_points);
    sv.setMutexForPoints(&mtx_new_pts);
    sv.setMutexForEdges(&mtx_new_edges);

    for (unsigned short i = 0; i < rl; i++) {
        auto start_refine_rl_time = chrono::high_resolution_clock::now();

        //the new_pts is a list that holds the coordinates of
        //new points inserted at this iteration. At the end of
        //this bucle, they are inserted in the point vector
        new_pts.clear();

        nb_points = points.size();

        list<RefinementRegion *>::const_iterator reg_iter;

        //split the Quadrants as needed
        //begin parallelisation ?
        //new_pts is shared, and new_quadrants

#pragma parallel for shared(new_points, new_Quadrants, nb_points, mtx_new_pts, mtx_new_edges) first_private(reg_iter, sv)
        for (int j = 0; j < tmp_Quadrants.size(); ++j) {

            Quadrant &iter = tmp_Quadrants[j];


            //Only check, can not modify after treatment
            bool to_refine = false;

            for (reg_iter = all_reg.begin(), reg_iter++; reg_iter != all_reg.end(); ++reg_iter) {


                unsigned short region_rl = (*reg_iter)->getRefinementLevel();
                if (region_rl < i) {
                    continue;
                }

                //If the Quadrant has a greater RL than the region needs, continue
                if (region_rl <= iter.getRefinementLevel()) {
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
                if ((*reg_iter)->intersectsQuadrant(points, iter)) {
                    to_refine = true;
                    //counterRefine.fetch_and_increment();
                    break;
                }
            }


            //now if refinement is not needed, we add the Quadrant as it was.
            if (!to_refine) {
                new_Quadrants.push_back(iter); //SHARED VARIABLE
                //End of for loop
            } else {

                list<unsigned int> &inter_edges = iter.getIntersectedEdges();
                unsigned short qrl = iter.getRefinementLevel();

                vector<vector<Point3D> > clipping_coords;
                sv.setClipping(clipping_coords);

                vector<vector<unsigned int> > split_elements;
                sv.setNewEles(split_elements);

                iter.accept(&sv);

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
        }

        // don't forget to update list
        tmp_Quadrants.assign(make_move_iterator(new_Quadrants.begin()), make_move_iterator(new_Quadrants.end()));
        new_Quadrants.clear();

        //if no points were added at this iteration, it is no longer
        //necessary to continue the refinement.
        if (new_pts.empty()) {
            cout << "warning at Mesher::generateQuadtreeMesh no new points!!!\n";
            break;
        }

        //add the new points to the vector
        points.reserve(points.size() + new_pts.size());
        points.insert(points.end(), new_pts.begin(), new_pts.end());

        auto end_refine_rl_time = chrono::high_resolution_clock::now();
        cout << "         * level " << i << " in "
             << std::chrono::duration_cast<chrono::milliseconds>(end_refine_rl_time - start_refine_rl_time).count();
        cout << " ms" << endl;

    } //END FOR REFINEMENT LEVEL

    std::cout << counterRefine <<
              std::endl;
}