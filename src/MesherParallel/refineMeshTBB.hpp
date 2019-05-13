#include "../Mesher.h"

#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/task_group.h>
#include <tbb/concurrent_vector.h>
#include <tbb/concurrent_unordered_set.h>
#include "ParallelTest1TBB/SplitVisitorTest1TBB.h"
#include "ParallelReductionTBB/ParallelReduceTBB.hpp"
#include "ParallelCustomReductionTBBV2/ParallelReduceTBBV2.hpp"
#include "ParallelCustomReductionTBBV3/ParallelReduceTBBV3.hpp"
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

    string Mesher::refineMeshParallelTest1TBB(int nbThread, list<Quadrant> Quadrants, vector<MeshPoint> points,
                                            set<QuadEdge> QuadEdges,
                                            const list<RefinementRegion *> &all_reg, const unsigned short &rl,
                                            Polyline &input) {

        int NOMBRE_THREAD = tbb::task_scheduler_init::default_num_threads();

        //OpenMP
        if (nbThread > NOMBRE_THREAD || nbThread < 0) {
            std::cout << "Invalid number of threads or not supported by computer" << std::endl;
            return;
        }

        tbb::task_scheduler_init test(nbThread);

        std::cout << "Start refine mesh parallel test 1 IntelTBB with " << nbThread << " threads." << std::endl;

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
            cout << "\t* level " << i << " in "
                 << total;
            cout << " ms" << endl;

            std::cout << "\t---- Points : " << points.size() << std::endl;
            std::cout << "\t---- QuadEdges : " << quadEdges.size() << std::endl;
            std::cout << "\t---- Quadrants : " << tmp_Quadrants.size() << std::endl;


        } //END FOR REFINEMENT LEVEL

        std::cout << "----------------------------------------------------------" << std::endl;
        std::cout << "----------------------------------------------------------" << std::endl;
        std::cout << "---------------------END OF TEST 1 INTELTBB---------------" << std::endl;
        std::cout << "----------------------------------------------------------" << std::endl;
        std::cout << "----------------------------------------------------------" << std::endl;

    }


    string Mesher::refineMeshReductionTBB(int nbThread, list<Quadrant> &tmp_Quadrants, vector<MeshPoint> &points,
                                        set<QuadEdge> &QuadEdges,
                                        const list<RefinementRegion *> &all_reg, const unsigned short &rl,
                                        Polyline &input) {

        int NOMBRE_THREAD = tbb::task_scheduler_init::default_num_threads();

        if (nbThread > NOMBRE_THREAD || nbThread < 0) {
            std::cout << "Invalid number of threads or not supported by computer" << std::endl;
            return;
        }

        tbb::task_scheduler_init test(nbThread);

        std::cout << "Start refine mesh reduction IntelTBB with " << nbThread << " threads." << std::endl;

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
            long total = std::chrono::duration_cast<chrono::milliseconds>(
                    end_refine_rl_time - start_refine_rl_time).count();
            cout << "\t* level " << i << " in "
                 << total;
            cout << " ms" << endl;

            std::cout << "\t---- Points : " << tmp_points.size() << std::endl;
            std::cout << "\t---- QuadEdges : " << tmp_edges.size() << std::endl;
            std::cout << "\t---- Quadrants : " << tmp_quadrants.size() << std::endl;


        }

        std::cout << "----------------------------------------------------------" << std::endl;
        std::cout << "----------------------------------------------------------" << std::endl;
        std::cout << "---------------------END OF REDUCTION TBB--------------" << std::endl;
        std::cout << "----------------------------------------------------------" << std::endl;
        std::cout << "----------------------------------------------------------" << std::endl;


    }


    string Mesher::refineMeshCustomReductionTBB(int nbThread, list<Quadrant> &tmp_Quadrants, vector<MeshPoint> &points,
                                              set<QuadEdge> &QuadEdges,
                                              const list<RefinementRegion *> &all_reg, const unsigned short &rl,
                                              Polyline &input) {

        int NOMBRE_THREAD = tbb::task_scheduler_init::default_num_threads();

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

            int remainingQuads = tmp_quadrants.size();
            int prevStart = 0;

            for (int j = 0; j < nbThread && remainingQuads > 0; j++) {

                //Create thread not master
                threads.push_back(
                        new RefineMeshReduction(i, tmp_quadrants, tmp_edges, input, tmp_points, all_reg, j == 0));

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
                delete threads[j];
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
                auto found = tmp_edges.insert(edge);
                if (!found.second) {
                    tmp_edges.erase(found.first);
                    tmp_edges.insert(edge);
                }
            }
            //tmp_edges.insert(rmr.getNewEdges().begin(), rmr.getNewEdges().end());


            delete threads[0];

            auto end_refine_rl_time = chrono::high_resolution_clock::now();
            long total = std::chrono::duration_cast<chrono::milliseconds>(
                    end_refine_rl_time - start_refine_rl_time).count();
            cout << "\t* level " << i << " in "
                 << total;
            cout << " ms" << endl;

            //long outside = std::chrono::duration_cast<chrono::milliseconds>(end_outside_block_time - start_outside_block_time).count();
            //cout << "TBB for outside / inside " << outside << " ms (" << (outside * 100.0 / total) << "%) ";
            //cout << " split visitor " << time_split_visitor << " ms (" << (time_split_visitor * 100.0 / outside) << "% of time) ";
            //cout << endl;

            std::cout << "\t---- Points : " << tmp_points.size() << std::endl;
            std::cout << "\t---- QuadEdges : " << tmp_edges.size() << std::endl;
            std::cout << "\t---- Quadrants : " << tmp_quadrants.size() << std::endl;

        }

        std::cout << "----------------------------------------------------------" << std::endl;
        std::cout << "----------------------------------------------------------" << std::endl;
        std::cout << "---------------------END OF CUSTOM REDUCTION TBB-----------" << std::endl;
        std::cout << "----------------------------------------------------------" << std::endl;
        std::cout << "----------------------------------------------------------" << std::endl;
    }

    string Mesher::refineCustomMeshReductionTBBV2(int nbThread, list<Quadrant> &tmp_Quadrants, vector<MeshPoint> &points,
                                                set<QuadEdge> &QuadEdges,
                                                const list<RefinementRegion *> &all_reg, const unsigned short &rl,
                                                Polyline &input) {

        int NOMBRE_THREAD = tbb::task_scheduler_init::default_num_threads();

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

            //threads.push_back(new RefineMeshReductionV2(i, tmp_quadrants, tmp_edges, input, tmp_points, all_reg, true));

            int remainingQuads = tmp_quadrants.size();
            int prevStart = 0;

            for (int j = 0; j < nbThread && remainingQuads > 0; j++) {
                threads.push_back(
                        new RefineMeshReductionV2(i, tmp_quadrants, tmp_edges, input, tmp_points, all_reg, j == 0));

                if (remainingQuads < split) split = remainingQuads;

                remainingQuads -= split;


                tg.run([&threads, &tmp_quadrants, j, split, prevStart] { // run in task group
                    threads[j]->operator()(tbb::blocked_range<size_t>(prevStart, (prevStart + split)));
                });

                prevStart += split;
            }

            tg.wait();

            threads[0]->doMasterJoin();

            for (int j = 1; j < threads.size(); j++) {
                threads[0]->join(*threads[j]);
                delete threads[j];
            }

            RefineMeshReductionV2 &rmr = *threads[0];

            std::swap(tmp_quadrants, rmr.getNewQuadrants());

            //if no points were added at this iteration, it is no longer
            //necessary to continue the refinement.

            if (rmr.getNewPts().empty()) {
                cout << "warning at Mesher::generateQuadtreeMesh no new points!!!\n";
                break;
            }

            delete threads[0];

            auto end_refine_rl_time = chrono::high_resolution_clock::now();
            long total = std::chrono::duration_cast<chrono::milliseconds>(
                    end_refine_rl_time - start_refine_rl_time).count();
            cout << "\t* level " << i << " in "
                 << total;
            cout << " ms" << endl;

            std::cout << "\t---- Points : " << tmp_points.size() << std::endl;
            std::cout << "\t---- QuadEdges : " << tmp_edges.size() << std::endl;
            std::cout << "\t---- Quadrants : " << tmp_quadrants.size() << std::endl;

        }

        std::cout << "----------------------------------------------------------" << std::endl;
        std::cout << "----------------------------------------------------------" << std::endl;
        std::cout << "---------------------END OF REDUCTION V2 TBB--------------" << std::endl;
        std::cout << "----------------------------------------------------------" << std::endl;
        std::cout << "----------------------------------------------------------" << std::endl;

    }


    string Mesher::refineCustomMeshReductionTBBV3(int nbThread, list<Quadrant> &tmp_Quadrants, vector<MeshPoint> &points,
                                                set<QuadEdge> &QuadEdges,
                                                const list<RefinementRegion *> &all_reg, const unsigned short &rl,
                                                Polyline &input) {

        string result = "";

        int NOMBRE_THREAD = tbb::task_scheduler_init::default_num_threads();

        if (nbThread > NOMBRE_THREAD || nbThread < 0) {
            result = "Invalid number of threads or not supported by computer\n";
            return result;
        }

        // TEST REDUCTION
        vector<MeshPoint> tmp_points(points.begin(), points.end());
        tbb::concurrent_unordered_set<QuadEdge, std::hash<QuadEdge>> tmp_edges(QuadEdges.begin(), QuadEdges.end());
        tbb::concurrent_vector<Quadrant> tmp_quadrants(tmp_Quadrants.begin(), tmp_Quadrants.end());

        tbb::task_scheduler_init test(nbThread);
        tbb::task_group tg;

        for (unsigned short i = 0; i < rl; i++) {

            // !! START SPLIT
            auto start_refine_rl_time = chrono::high_resolution_clock::now();
            auto start_split = chrono::high_resolution_clock::now();

            long split = tmp_quadrants.size() / (nbThread) + 1;
            split = std::max(split, 5000l);
            //std::cout << split << "/" << tmp_quadrants.size() << std::endl;

            vector<RefineMeshReductionV3 *> threads;

            //threads.push_back(new RefineMeshReductionV3(i, tmp_quadrants, tmp_edges, input, tmp_points, all_reg, true));

            long remainingQuads = tmp_quadrants.size();
            int prevStart = 0;

            for (int j = 0; j < nbThread && remainingQuads > 0; j++) {
                threads.push_back(
                        new RefineMeshReductionV3(i, tmp_quadrants, tmp_edges, input, tmp_points, all_reg, j == 0));

                if (remainingQuads < split) split = remainingQuads;
                remainingQuads -= split;

                tg.run([&threads, j, split, prevStart] { // run in task group
                    threads[j]->operator()(tbb::blocked_range<size_t>(prevStart, (prevStart + split)));
                });

                prevStart += split;
            }

            tg.wait();

            //auto end_split = chrono::high_resolution_clock::now();

            //long total1 = std::chrono::duration_cast<chrono::milliseconds>(end_split - start_split).count();
            //cout << " time split " << total1 << endl;


            // !! START JOIN
            //auto start_join = chrono::high_resolution_clock::now();

            // init map vector to map local point index to global
            std::vector<std::tr1::unordered_map<size_t, unsigned int>> threadToGlobal(threads.size() - 1);

            unsigned long old_points_size = tmp_points.size();

            // compute new points index sequentially (for now ?)
            for (int t = 0; t < threads.size(); t++) {
                if (t == 0) {
                    //add the new points to the vector
                    tmp_points.reserve(tmp_points.size() + threads[0]->getNewPts().size());
                    tmp_points.insert(tmp_points.end(), threads[0]->getNewPts().begin(), threads[0]->getNewPts().end());
                } else {
                    // compute new index and add if necessary
                    int counter = 0;
                    for (const Point3D &point : threads[t]->getNewPts()) {
                        std::tr1::unordered_map<size_t, unsigned int> &threadMap = threadToGlobal[t - 1];
                        std::tr1::unordered_map<size_t, unsigned int> &masterMap = threads[0]->getNewMaps();

                        // compute hash
                        size_t hashPoint = point.operator()(point);

                        auto found = masterMap.insert(std::pair<size_t, unsigned int>(hashPoint, tmp_points.size()));

                        if (!found.second) {
                            // point already exists
                            threadMap[counter++ + old_points_size] = (*found.first).second; // get index
                        } else {
                            // point inserted in map
                            threadMap[counter++ + old_points_size] = tmp_points.size();
                            tmp_points.emplace_back(point);
                        }
                    }
                }
            }

            //if no points were added at this iteration, it is no longer
            //necessary to continue the refinement.

            if (old_points_size == tmp_points.size()) {
                result += "warning at Mesher::generateQuadtreeMesh no new points!!!\n";
                break;
            }

            //auto end_join_points = chrono::high_resolution_clock::now();

            //long total4 = std::chrono::duration_cast<chrono::milliseconds>(end_join_points - start_join).count();
            //cout << " time points " << total4 << endl;

            //tbb::task_scheduler_init test2(nbThread);

            // now compute new edges and quads in parallel
            tmp_quadrants.clear();


            for (unsigned int t = 1; t < threads.size(); t++) {
                tg.run([&, t] { // run in task group
                    //std::cout << "Edge start" << std::endl;
                    //auto start_edge = chrono::high_resolution_clock::now();

                    RefineMeshReductionV3 *rmr = threads[t];
                    std::tr1::unordered_map<size_t, unsigned int> &threadMap = threadToGlobal[t - 1];

                    for (const QuadEdge &local_edge : rmr->getNewEdges()) {
                        // build new edge with right index
                        vector<unsigned int> index(3, 0);

                        for (unsigned int i = 0; i < 3; i++) {
                            if (local_edge[i] < old_points_size) {
                                // index refer point not created during this refinement level
                                index[i] = local_edge[i];
                            } else {
                                // point created locally, need to update the point with correct index
                                index[i] = threadMap[local_edge[i]];
                            }
                        }

                        QuadEdge edge(index[0], index[1], index[2]);

                        auto found = tmp_edges.insert(edge); // try insert


                        // if edge already exists
                        if (!found.second) {
                            if (edge[2] != 0 && edge[2] != (*found.first)[2]) {
                                // since all points have been replaced, if it's different then midpoint has been created
                                (found.first)->updateMidPoint(edge[2]);
                            }
                        }


                    }
                    //auto end_edge = chrono::high_resolution_clock::now();

                    //long total2 = std::chrono::duration_cast<chrono::milliseconds>(end_edge - start_edge).count();
                    //cout << " time edges " << total2  << endl;

                    //std::cout << "Edge end" << std::endl;

                    //std::cout << "Quad start" << std::endl;

                    //auto start_quad = chrono::high_resolution_clock::now();
                    for (const Quadrant &local_quad : rmr->getNewQuadrants()) {
                        // build new quad with right index
                        vector<unsigned int> new_pointindex(4, 0);

                        for (unsigned int j = 0; j < 4; j++) {
                            if (local_quad.getPointIndex(j) < old_points_size) {
                                // index refer point not created during this refinement level
                                new_pointindex[j] = local_quad.getPointIndex(j);
                            } else {
                                // point created, need to update the point with correct index
                                new_pointindex[j] = threadMap[local_quad.getPointIndex(j)];
                            }
                        }

                        tmp_quadrants.emplace_back(new_pointindex, local_quad);
                    }

                    //auto end_quad = chrono::high_resolution_clock::now();

                    //long total3 = std::chrono::duration_cast<chrono::milliseconds>(end_quad - start_quad).count();
                    //cout << " time quad " << total3 << endl;

                    //std::cout << "Quad end" << std::endl;
                });
            }

            //add the new edges of "master" to the vector
            for (const QuadEdge &edge : threads[0]->getNewEdges()) {
                auto found = tmp_edges.insert(edge); // try insert

                // if edge already exists
                if (!found.second) {
                    if (edge[2] != 0 && edge[2] != (*found.first)[2]) {
                        // since all points have been replaced, if it's different then midpoint has been created
                        (found.first)->updateMidPoint(edge[2]);
                    }
                }
            }

            // add new quads of master
            for (const Quadrant &quad : threads[0]->getNewQuadrants()) {
                tmp_quadrants.push_back(quad);
            }

            // Wait for completion of the task group
            tg.wait();

            for (auto &thread : threads) {
                delete thread;
            }

            //auto end_join = chrono::high_resolution_clock::now();

            //long total2 = std::chrono::duration_cast<chrono::milliseconds>(end_join - end_split).count();
            //cout << " time join " << total2 << endl;

            auto end_refine_rl_time = chrono::high_resolution_clock::now();
            long total = std::chrono::duration_cast<chrono::milliseconds>(
                    end_refine_rl_time - start_refine_rl_time).count();
            
            result += "Level " + i + " in " + total + " ms\n";
            result += "Points : " + std::to_string(tmp_points.size()) + "\n";
            result += "QuadEdges : " + std::to_string(tmp_edges.size()) + "\n";
            result += "Quadrants : " + std::to_string(tmp_quadrants.size()) + "\n";
        }

        return result;

        // std::cout << "----------------------------------------------------------" << std::endl;
        // std::cout << "----------------------------------------------------------" << std::endl;
        // std::cout << "---------------------END OF REDUCTION V3 TBB--------------" << std::endl;
        // std::cout << "----------------------------------------------------------" << std::endl;
        // std::cout << "----------------------------------------------------------" << std::endl;

    }

    string Mesher::refineMeshCustomReductionTBBV4(int nbThread, list<Quadrant> &tmp_Quadrants, vector<MeshPoint> &points,
                                                set<QuadEdge> &QuadEdges,
                                                const list<RefinementRegion *> &all_reg, const unsigned short &rl,
                                                Polyline &input) {
        string result = "";

        int NOMBRE_THREAD = tbb::task_scheduler_init::default_num_threads();

        if (nbThread > NOMBRE_THREAD || nbThread < 0) {
            result = "Invalid number of threads or not supported by computer\n";
            return result;
        }

        //std::cout << "Start refine mesh custom reduction V4 IntelTBB with " << nbThread << " threads." << std::endl;

        // TEST REDUCTION
        vector<MeshPoint> tmp_points(points.begin(), points.end());
        set<QuadEdge> tmp_edges(QuadEdges.begin(), QuadEdges.end());
        vector<Quadrant> tmp_quadrants(tmp_Quadrants.begin(), tmp_Quadrants.end());

        tbb::task_scheduler_init test(nbThread);
        tbb::task_group tg;

        for (unsigned short i = 0; i < rl; i++) {

            // !! START SPLIT
            auto start_refine_rl_time = chrono::high_resolution_clock::now();
            auto start_split = chrono::high_resolution_clock::now();

            long split = tmp_quadrants.size() / (nbThread) + 1;
            split = std::max(split, 5000l);
            //std::cout << split << "/" << tmp_quadrants.size() << std::endl;

            vector<RefineMeshReductionV2 *> threads;

            threads.push_back(new RefineMeshReductionV2(i, tmp_quadrants, tmp_edges, input, tmp_points, all_reg, true));

            long remainingQuads = tmp_quadrants.size();
            int prevStart = 0;

            for (int j = 0; j < nbThread && remainingQuads > 0; j++) {
                if (j != 0)
                    threads.push_back(
                            new RefineMeshReductionV2(i, tmp_quadrants, tmp_edges, input, tmp_points, all_reg, false));

                if (remainingQuads < split) split = remainingQuads;
                remainingQuads -= split;

                tg.run([&threads, j, split, prevStart] { // run in task group
                    threads[j]->operator()(tbb::blocked_range<size_t>(prevStart, (prevStart + split)));
                });

                prevStart += split;
            }

            tg.wait();

            //auto end_split = chrono::high_resolution_clock::now();

            //long total1 = std::chrono::duration_cast<chrono::milliseconds>(end_split - start_split).count();
            //cout << " time split " << total1 << endl;


            // !! START JOIN
            //auto start_join = chrono::high_resolution_clock::now();

            tmp_quadrants.clear();

            // init map vector to map local point index to global
            std::vector<std::tr1::unordered_map<unsigned int, unsigned int>> threadToGlobal(threads.size());

            unsigned long old_points_size = tmp_points.size();

            // compute new points index sequentially
            for (int t = 0; t < threads.size(); t++) {

                auto start_edge = chrono::high_resolution_clock::now();

                // compute new index and add if necessary
                std::tr1::unordered_map<unsigned int, unsigned int> &threadMap = threadToGlobal[t];

                for (const QuadEdge &local_edge : threads[t]->getNewEdges()) {

                    if (local_edge[0] < old_points_size && local_edge[1] < old_points_size) {
                        // midpoint created on existing thread
                        auto found = tmp_edges.find(local_edge);

                        if ((*found)[2] != 0) {
                            // midpoint already created by another thread
                            threadMap.insert(pair<unsigned long, unsigned long>(local_edge[2], (*found)[2]));
                        } else {
                            // midpoint created only by this thread for now
                            // add midpoint in point and update edge midpoint
                            //(*found).updateMidPoint(tmp_points.size());
                            QuadEdge quadEdge(local_edge);
                            quadEdge.updateMidPoint(tmp_points.size());
                            auto pos = tmp_edges.erase(found);
                            tmp_edges.insert(pos, quadEdge);

                            threadMap.insert(pair<unsigned long, unsigned long>(local_edge[2], tmp_points.size()));
                            tmp_points.emplace_back(threads[t]->getNewPts()[local_edge[2] - old_points_size]);
                        }
                    } else {
                        // it not about midpoint created on "old" edge but totally new edge
                        // check if index need to be created or change

                        // build new edge with correct index
                        vector<unsigned int> index(3, 0);

                        for (unsigned int j = 0; j < 3; j++) {
                            if (local_edge[j] < old_points_size) {
                                // index refer point not created during this refinement level
                                index[j] = local_edge[j];
                            } else {
                                // point created locally, need to update the point with correct index
                                if (threadMap.find(local_edge[j]) != threadMap.end()) {
                                    index[j] = threadMap[local_edge[j]];
                                } else {
                                    index[j] = tmp_points.size();
                                    threadMap.insert(
                                            pair<unsigned int, unsigned int>(local_edge[j], tmp_points.size()));
                                    tmp_points.emplace_back(threads[t]->getNewPts()[local_edge[j] - old_points_size]);
                                }
                            }
                        }

                        QuadEdge edge(index[0], index[1], index[2]);

                        tmp_edges.insert(edge); // try insert
                    }
                }

                //auto end_edge = chrono::high_resolution_clock::now();

                //long total2 = std::chrono::duration_cast<chrono::milliseconds>(end_edge - start_edge).count();
                //cout << " time edges & points " << total2  << endl;


                //auto start_quad = chrono::high_resolution_clock::now();
                for (Quadrant &local_quad : threads[t]->getNewQuadrants()) {
                    // build new quad with right index
                    vector<unsigned int> new_pointindex(4, 0);

                    for (unsigned int j = 0; j < 4; j++) {
                        if (local_quad.getPointIndex(j) < old_points_size) {
                            // index refer point not created during this refinement level
                            new_pointindex[j] = local_quad.getPointIndex(j);
                        } else {
                            // point created, need to update the point with correct index
                            new_pointindex[j] = threadMap[local_quad.getPointIndex(j)];
                        }
                    }

                    tmp_quadrants.emplace_back(new_pointindex, local_quad);
                }

                //auto end_quad = chrono::high_resolution_clock::now();

                //long total3 = std::chrono::duration_cast<chrono::milliseconds>(end_quad - start_quad).count();
                //cout << " time quad " << total3 << endl;

            }

            //auto end_join = chrono::high_resolution_clock::now();

            //long total2 = std::chrono::duration_cast<chrono::milliseconds>(end_join - end_split).count();
            //cout << " time join " << total2 << endl;

            auto end_refine_rl_time = chrono::high_resolution_clock::now();
            long total = std::chrono::duration_cast<chrono::milliseconds>(
                    end_refine_rl_time - start_refine_rl_time).count();

            result += "Level " + i + " in " + total + " ms\n";
            result += "Points : " + std::to_string(tmp_points.size()) + "\n";
            result += "QuadEdges : " + std::to_string(tmp_edges.size()) + "\n";
            result += "Quadrants : " + std::to_string(tmp_quadrants.size()) + "\n";

        }

        // std::cout << "----------------------------------------------------------" << std::endl;
        // std::cout << "----------------------------------------------------------" << std::endl;
        // std::cout << "---------------------END OF REDUCTION V4 TBB--------------" << std::endl;
        // std::cout << "----------------------------------------------------------" << std::endl;
        // std::cout << "----------------------------------------------------------" << std::endl;

    }


}