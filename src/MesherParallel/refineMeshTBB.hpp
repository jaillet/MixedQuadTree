#include "../Mesher.h"

#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/task_group.h>
#include <tbb/concurrent_vector.h>
#include <tbb/concurrent_unordered_set.h>
#include "ParallelTest1TBB/SplitVisitorTest1TBB.h"
#include "ParallelReductionTBB/ParallelReduceTBB.hpp"
#include "ParallelCustomReductionV2TBB/ParallelReduceTBBV2.hpp"
#include <atomic>

namespace std {
    template<>
    struct hash<Clobscode::QuadEdge> {
        size_t operator()(const Clobscode::QuadEdge &k) const {
            // Compute individual hash values for two data members and combine them using XOR and bit shifting
            return ((hash<int>()(k[0]) ^ (hash<int>()(k[1]) << 1)) >> 1);
        }
    };
}


namespace Clobscode {


    void Mesher::refineMeshReductionTBB(int nbThread, list<Quadrant> &tmp_Quadrants, vector<MeshPoint> &points,
                                        set<QuadEdge> &QuadEdges,
                                        const list<RefinementRegion *> &all_reg, const unsigned short &rl,
                                        Polyline &input) {

        int NOMBRE_THREAD = tbb::task_scheduler_init::default_num_threads();
        std::cout << NOMBRE_THREAD << std::endl;

        if (nbThread > NOMBRE_THREAD || nbThread < 0) {
            std::cout << "Invalid number of threads or not supported by computer" << std::endl;
            return;
        }

        tbb::task_scheduler_init test(nbThread);

        // TEST REDUCTION
        list<Point3D> new_pts;
        vector<MeshPoint> tmp_points(points.begin(), points.end());
        set<QuadEdge> tmp_edges(QuadEdges.begin(), QuadEdges.end());
        vector<Quadrant> tmp_quadrants(tmp_Quadrants.begin(), tmp_Quadrants.end());


        for (unsigned short i = 0; i < rl; i++) {
            auto start_refine_rl_time = chrono::high_resolution_clock::now();

            new_pts.clear();

            int split = tmp_quadrants.size() / nbThread + 1;
            split = std::max(split, 10000);

            RefineMeshReduction rmr(i, tmp_quadrants, tmp_edges, input, tmp_points, all_reg, true);
            parallel_reduce(tbb::blocked_range<size_t>(0, tmp_quadrants.size(), split), rmr);


            std::swap(tmp_quadrants, rmr.getNewQuadrants());

            //if no points were added at this iteration, it is no longer
            //necessary to continue the refinement.
            if (rmr.getNewPts().empty()) {
                cout << "warning at Mesher::generateQuadtreeMesh no new points!!!\n";
                break;
            }

            //add the new points to the vector
            tmp_points.reserve(tmp_points.size() + rmr.getNewPts().size());
            tmp_points.insert(tmp_points.end(), rmr.getNewPts().begin(), rmr.getNewPts().end());

            //add the new edges to the vector
            for (auto edge : rmr.getNewEdges()) {
                auto found = tmp_edges.find(edge);
                if (found != tmp_edges.end()) {
                    tmp_edges.erase(found);
                    tmp_edges.insert(edge);
                } else {
                    tmp_edges.insert(edge);
                }
            }

            auto end_refine_rl_time = chrono::high_resolution_clock::now();
            long total = std::chrono::duration_cast<chrono::milliseconds>(end_refine_rl_time - start_refine_rl_time).count();
            cout << "         * level " << i << " in "
                 << total;
            cout << " ms" << endl;

            std::cout << "           ---- Points : " << tmp_points.size() << std::endl;
            std::cout << "           ---- QuadEdges : " << tmp_edges.size() << std::endl;
            std::cout << "           ---- Quadrants : " << tmp_quadrants.size() << std::endl;


        }

        std::cout << "----------------------------------------------------------" << std::endl;
        std::cout << "----------------------------------------------------------" << std::endl;
        std::cout << "---------------------END OF REDUCTION TBB--------------" << std::endl;
        std::cout << "----------------------------------------------------------" << std::endl;
        std::cout << "----------------------------------------------------------" << std::endl;


    }

    void Mesher::refineCustomMeshReductionTBBV2(int nbThread, list<Quadrant> &tmp_Quadrants, vector<MeshPoint> &points,
                                                set<QuadEdge> &QuadEdges,
                                                const list<RefinementRegion *> &all_reg, const unsigned short &rl,
                                                Polyline &input) {

        int NOMBRE_THREAD = tbb::task_scheduler_init::default_num_threads();
        std::cout << NOMBRE_THREAD << std::endl;

        if (nbThread > NOMBRE_THREAD || nbThread < 0) {
            std::cout << "Invalid number of threads or not supported by computer" << std::endl;
            return;
        }

        tbb::task_scheduler_init test(nbThread);
        tbb::task_group tg;


        // TEST REDUCTION
        list<Point3D> new_pts;
        vector<MeshPoint> tmp_points(points.begin(), points.end());
        set<QuadEdge> tmp_edges(QuadEdges.begin(), QuadEdges.end());
        vector<Quadrant> tmp_quadrants(tmp_Quadrants.begin(), tmp_Quadrants.end());


        for (unsigned short i = 0; i < rl; i++) {
            auto start_refine_rl_time = chrono::high_resolution_clock::now();

            int split = tmp_quadrants.size() / (nbThread) + 1;
            split = std::max(split, 5000);

            vector<RefineMeshReductionV2 *> threads;

            threads.push_back(new RefineMeshReductionV2(i, tmp_quadrants, tmp_edges, input, tmp_points, all_reg, true));

            int remainingQuads = tmp_quadrants.size();
            int prevStart = 0;

            for (int j = 0; j < nbThread && remainingQuads > 0; j++) {
                if (j != 0) {
                    threads.push_back(new RefineMeshReductionV2(i, tmp_quadrants, tmp_edges, input, tmp_points, all_reg, false));
                }

                if (remainingQuads < split) split = remainingQuads;

                remainingQuads -= split;


                tg.run([&threads, &tmp_quadrants, j, split, prevStart]{ // run in task group
                    threads[j]->operator()(tbb::blocked_range<size_t>(prevStart, (prevStart + split)));
                });

                prevStart += split;
            }

            tg.wait();

            threads[0]->doMasterJoin();

            for (int j = 1; j < threads.size(); j++) {
                threads[0]->join(*threads[j]);
            }

            RefineMeshReductionV2 & rmr = *threads[0];

            std::swap(tmp_quadrants, rmr.getNewQuadrants());

            //if no points were added at this iteration, it is no longer
            //necessary to continue the refinement.

            if (rmr.getNewPts().empty()) {
                cout << "warning at Mesher::generateQuadtreeMesh no new points!!!\n";
                break;
            }

            auto end_refine_rl_time = chrono::high_resolution_clock::now();
            long total = std::chrono::duration_cast<chrono::milliseconds>(end_refine_rl_time - start_refine_rl_time).count();
            cout << "         * level " << i << " in "
                 << total;
            cout << " ms" << endl;

            std::cout << "           ---- Points : " << tmp_points.size() << std::endl;
            std::cout << "           ---- QuadEdges : " << tmp_edges.size() << std::endl;
            std::cout << "           ---- Quadrants : " << tmp_quadrants.size() << std::endl;

        }

        std::cout << "----------------------------------------------------------" << std::endl;
        std::cout << "----------------------------------------------------------" << std::endl;
        std::cout << "---------------------END OF REDUCTION V2 TBB--------------" << std::endl;
        std::cout << "----------------------------------------------------------" << std::endl;
        std::cout << "----------------------------------------------------------" << std::endl;

    }

    void Mesher::refineMeshCustomReductionTBB(int nbThread, list<Quadrant> &tmp_Quadrants, vector<MeshPoint> &points,
                                              set<QuadEdge> &QuadEdges,
                                              const list<RefinementRegion *> &all_reg, const unsigned short &rl,
                                              Polyline &input) {

        int NOMBRE_THREAD = tbb::task_scheduler_init::default_num_threads();
        std::cout << NOMBRE_THREAD << std::endl;

        if (nbThread > NOMBRE_THREAD || nbThread < 0) {
            std::cout << "Invalid number of threads or not supported by computer" << std::endl;
            return;
        }

        tbb::task_scheduler_init test(nbThread);
        tbb::task_group tg;

        std::cout << "Start refine mesh custom reduction V1 IntelTBB with " << nbThread << " threads." << std::endl;

        // TEST REDUCTION
        list<Point3D> new_pts;
        vector<MeshPoint> tmp_points(points.begin(), points.end());
        set<QuadEdge> tmp_edges(QuadEdges.begin(), QuadEdges.end());
        vector<Quadrant> tmp_quadrants(tmp_Quadrants.begin(), tmp_Quadrants.end());


        for (unsigned short i = 0; i < rl; i++) {
            auto start_refine_rl_time = chrono::high_resolution_clock::now();

            new_pts.clear();

            int split = tmp_quadrants.size() / (nbThread) + 1;
            split = std::max(split, 5000);
            //std::cout << split << std::endl;

            vector<RefineMeshReduction *> threads;


            //threads.emplace_back(i, tmp_quadrants, tmp_edges, input, tmp_points, all_reg);
            //Create the master thread (can fill final structure in join)
            threads.push_back(new RefineMeshReduction(i, tmp_quadrants, tmp_edges, input, tmp_points, all_reg, true));
            //RefineMeshReduction rmr1 (i, tmp_quadrants, tmp_edges, input, tmp_points, all_reg);
            //RefineMeshReduction rmr2 (i, tmp_quadrants, tmp_edges, input, tmp_points, all_reg);

            int remainingQuads = tmp_quadrants.size();
            int prevStart = 0;

            for (int j = 0; j < nbThread && remainingQuads > 0; j++) {
                if (j != 0) {
                    //threads.emplace_back(i, tmp_quadrants, tmp_edges, input, tmp_points, all_reg);
                    //Create thread not master
                    threads.push_back(
                            new RefineMeshReduction(i, tmp_quadrants, tmp_edges, input, tmp_points, all_reg, false));
                }

                if (remainingQuads < split) split = remainingQuads;

                remainingQuads -= split;


                tg.run([&threads, &tmp_quadrants, j, split, prevStart] { // run in task group
                    //std::cout << "Start from " << prevStart << " to " << prevStart + split << " / " << tmp_quadrants.size() << std::endl;
                    threads[j]->operator()(tbb::blocked_range<size_t>(prevStart, (prevStart + split)));
                });

                prevStart += split;
            }

            tg.wait();

            for (int j = 1; j < threads.size(); j++) {
                threads[0]->join(*threads[j]);
            }

            RefineMeshReduction &rmr = *threads[0];

            //parallel_reduce(tbb::blocked_range<size_t>(0, tmp_quadrants.size(), split), rmr);


            std::swap(tmp_quadrants, rmr.getNewQuadrants());

            //if no points were added at this iteration, it is no longer
            //necessary to continue the refinement.
            if (rmr.getNewPts().empty()) {
                cout << "warning at Mesher::generateQuadtreeMesh no new points!!!\n";
                break;
            }

            //add the new points to the vector
            tmp_points.reserve(tmp_points.size() + rmr.getNewPts().size());
            tmp_points.insert(tmp_points.end(), rmr.getNewPts().begin(), rmr.getNewPts().end());

            //add the new edges to the vector
            for (auto edge : rmr.getNewEdges()) {
                auto found = tmp_edges.find(edge);
                if (found != tmp_edges.end()) {
                    tmp_edges.erase(found);
                    tmp_edges.insert(edge);
                } else {
                    tmp_edges.insert(edge);
                }

            }
            //tmp_edges.insert(rmr.getNewEdges().begin(), rmr.getNewEdges().end());



            auto end_refine_rl_time = chrono::high_resolution_clock::now();
            long total = std::chrono::duration_cast<chrono::milliseconds>(
                    end_refine_rl_time - start_refine_rl_time).count();
            cout << "         * level " << i << " in "
                 << total;
            cout << " ms" << endl;

            //long outside = std::chrono::duration_cast<chrono::milliseconds>(end_outside_block_time - start_outside_block_time).count();
            //cout << "TBB for outside / inside " << outside << " ms (" << (outside * 100.0 / total) << "%) ";
            //cout << " split visitor " << time_split_visitor << " ms (" << (time_split_visitor * 100.0 / outside) << "% of time) ";
            //cout << endl;

            std::cout << "           ---- Points : " << tmp_points.size() << std::endl;
            std::cout << "           ---- QuadEdges : " << tmp_edges.size() << std::endl;
            std::cout << "           ---- Quadrants : " << tmp_quadrants.size() << std::endl;

        }

        std::cout << "----------------------------------------------------------" << std::endl;
        std::cout << "----------------------------------------------------------" << std::endl;
        std::cout << "---------------------END OF CUSTOM REDUCTION TBB-----------" << std::endl;
        std::cout << "----------------------------------------------------------" << std::endl;
        std::cout << "----------------------------------------------------------" << std::endl;
    }

    void Mesher::refineMeshParallelTest1TBB(int nbThread, list<Quadrant> Quadrants, vector<MeshPoint> points,
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
        //moving quadrants to save memory
        tmp_Quadrants.assign(make_move_iterator(Quadrants.begin()), make_move_iterator(Quadrants.end()));
        Quadrants.clear();

        tbb::concurrent_vector<Quadrant> new_Quadrants;

        //list of the points added at this refinement iteration:
        tbb::concurrent_vector<Point3D> new_pts;

        tbb::concurrent_unordered_set<QuadEdge, std::hash<QuadEdge>> quadEdges;
        quadEdges.insert(QuadEdges.begin(), QuadEdges.end());


        //Shared variables

        //points
        //--- Read in SplitVisitor::visit
        //QuadEdges Warning!
        //--- is filled by SplitVisitor::visit
        //--- and needed and filled in SplitVisitor::splitEdge!


        unsigned int nb_points = 0;

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

                                  //auto start_block_time = chrono::high_resolution_clock::now();


                                  //create visitors and give them variables
                                  SplitVisitorTest1TBB sv;
                                  sv.setPoints(points); // READ
                                  sv.setEdges(
                                          quadEdges); // TODO INSERT / REMOVE / READ A faire en priorite car le plus contraignant
                                  sv.setNewPts(new_pts); // TODO INSERT / READ
                                  sv.setCounterPoints(nb_points);

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

                                          // intersectQaudrant can modify the quadrant with
                                          // function Polyline::getNbFeatures in RefinementboundaryRegion
                                          // maybe no problem for parallelisation, as this information is
                                          // not used by other thread
                                          if ((*reg_iter)->intersectsQuadrant(points, iter)) {
                                              to_refine = true;
                                              break;
                                          }
                                      }


                                      //now if refinement is not needed, we add the Quadrant as it was.
                                      if (!to_refine) {
                                          new_Quadrants.push_back(iter); //SHARED VARIABLE
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

                                                  IntersectionsVisitor iv(true);
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
            long total = std::chrono::duration_cast<chrono::milliseconds>(
                    end_refine_rl_time - start_refine_rl_time).count();
            cout << "         * level " << i << " in "
                 << total;
            cout << " ms" << endl;

            std::cout << "           ---- Points : " << points.size() << std::endl;
            std::cout << "           ---- QuadEdges : " << quadEdges.size() << std::endl;
            std::cout << "           ---- Quadrants : " << tmp_Quadrants.size() << std::endl;


        } //END FOR REFINEMENT LEVEL

        std::cout << "----------------------------------------------------------" << std::endl;
        std::cout << "----------------------------------------------------------" << std::endl;
        std::cout << "---------------------END OF TEST 1 INTELTBB---------------" << std::endl;
        std::cout << "----------------------------------------------------------" << std::endl;
        std::cout << "----------------------------------------------------------" << std::endl;

    }


}