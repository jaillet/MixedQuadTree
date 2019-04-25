#include "CustomSplitVisitor.h"
#include <tbb/blocked_range.h>
#include <tr1/unordered_map>
#include <tbb/task_group.h>
#include <tbb/task_scheduler_init.h>

#include "../../Visitors/IntersectionsVisitor.h"
#include "../../RefinementRegion.h"
#include "../../Polyline.h"
#include "../../Quadrant.h"

/*
using Clobscode::CustomSplitVisitor;
using Clobscode::QuadEdge;
using Clobscode::RefinementRegion;
using Clobscode::Quadrant;
using Clobscode::MeshPoint;
using Clobscode::Polyline;
using Clobscode::Point3D;*/
using std::vector;
using std::list;
using std::set;

namespace Clobscode {

    class RefineMeshReductionV2 {

        //Private variables
        unsigned int m_rl;

        //Need to be read in reduction :
        vector<Quadrant> m_new_Quadrants;
        vector<Point3D> m_new_pts; //Local new pts (not merged)
        std::tr1::unordered_map<size_t, unsigned int> m_map_new_pts;
        set<QuadEdge> m_new_edges;

        //Read only
        const Polyline &input;
        vector<MeshPoint> &points;
        const list<RefinementRegion *> &all_reg;
        vector<Quadrant> &tmp_Quadrants;
        set<QuadEdge> &edges;

        CustomSplitVisitor csv;

        int numberOfJoint;
        bool master;
        bool joinDone;

        // usefull constante
        unsigned long nb_points;

        void setSplitVisitor() {
            csv.setPoints(points);
            csv.setEdges(edges);
            csv.setNewPts(m_new_pts);
            csv.setNewEdges(m_new_edges);
            csv.setMapPts(m_map_new_pts);
        }


    public:

        RefineMeshReductionV2(unsigned int refinementLevel, vector<Quadrant> &tmp_Quadrants, set<QuadEdge> &quadEdges,
                            Polyline &input, vector<MeshPoint> &points, const list<RefinementRegion *> &all_reg, const bool master) :
                m_rl(refinementLevel), input(input), points(points), all_reg(all_reg), tmp_Quadrants(tmp_Quadrants),
                edges(quadEdges), master(master){
            setSplitVisitor();
            numberOfJoint = 0;
            nb_points = points.size();
            joinDone = false;
        }

        /**
         * @brief Splitting constructor. Must be able to run concurrently with operator() and method join.
         * @details split is a dummy argument of type split, distinguishes the splitting constructor from a copy constructor.
         */
        RefineMeshReductionV2(RefineMeshReductionV2 &x, tbb::split) :
                m_rl(x.m_rl), input(x.input), points(x.points), all_reg(x.all_reg), tmp_Quadrants(x.tmp_Quadrants),
                edges(x.edges) {
            setSplitVisitor();
            numberOfJoint = 0;
            nb_points = points.size();
            joinDone = false;
        }

        /**
         * @brief insert points and edges in global structures
         */
        void doMasterJoin() {
            if (!joinDone) {

                //add the new points to the vector
                points.reserve(points.size() + m_new_pts.size());
                points.insert(points.end(), m_new_pts.begin(), m_new_pts.end());

                //add the new edges to the vector
                for (auto edge : m_new_edges) {
                    auto found = edges.find(edge);
                    if (found != edges.end()) {
                        found = edges.erase(found);
                        edges.insert(found, edge);
                    } else {
                        edges.insert(edge);
                    }

                }

                joinDone = true;
            }
        }

        /**
         * @brief Reduction.
         * @details Join results. The result in rmr should be merged into the result of this.
         */
        void join(const RefineMeshReductionV2 &rmr) {

            //std::cout << "Start join" << (master ? " master" : "") << std::endl;

            numberOfJoint += rmr.numberOfJoint + 1;

            // best than map for insert and access
            std::tr1::unordered_map<unsigned int, unsigned int> taskToGlobal;

            // allow master to insert directly in final structures
            if (master) {
                doMasterJoin();

                auto start_points = chrono::high_resolution_clock::now();

                int i = 0;
                for (const Point3D &point : rmr.m_new_pts) {

                    //auto found = m_new_pts.find(point);
                    size_t hashPoint = point.operator()(point);

                    auto found = m_map_new_pts.find(hashPoint);
                    if (found != m_map_new_pts.end()) {
                        taskToGlobal[i++ + nb_points] = m_map_new_pts[hashPoint];
                    } else {
                        m_map_new_pts[hashPoint] = points.size();
                        points.push_back(point);
                        taskToGlobal[i++ + nb_points] = points.size() - 1;
                    }
                }

                auto end_points = chrono::high_resolution_clock::now();

                long total1 = std::chrono::duration_cast<chrono::milliseconds>(end_points - start_points).count();
                cout << " time points " << total1 << endl;

                //tbb::task_scheduler_init def_init; // Use the default number of threads.
                tbb::task_group tg;

                tg.run([&] { // run in task group
                    //std::cout << "Edge start" << std::endl;
                    auto start_quad = chrono::high_resolution_clock::now();

                    long timeInit = 0;
                    long timeInsert = 0;
                    long timeInsert1 = 0;
                    long timeInsertErase = 0;
                    long timeFound = 0;

                    for (const QuadEdge &local_edge : rmr.m_new_edges) {
                        // build new edge with right index
                        auto start_quad = chrono::high_resolution_clock::now();

                        vector<unsigned long> index(3, 0);

                        for (unsigned int i = 0; i < 3; i++) {
                            if (local_edge[i] < nb_points) {
                                // index refer point not created during this refinement level
                                index[i] = local_edge[i];
                            } else {
                                // point created, need to update the point with correct index
                                index[i] = taskToGlobal[local_edge[i]];
                            }
                        }

                        QuadEdge edge(index[0], index[1], index[2]);

                        auto end_quad = chrono::high_resolution_clock::now();

                        timeInit += std::chrono::duration_cast<chrono::nanoseconds>(end_quad - start_quad).count();

                        start_quad = chrono::high_resolution_clock::now();

                        //auto found = edges.find(edge);
                        auto found = edges.insert(edge);

                        timeFound += std::chrono::duration_cast<chrono::nanoseconds>( chrono::high_resolution_clock::now() - start_quad).count();

                        if (!found.second) {
                            if (edge[2] != 0 && edge[2] != (*found.first)[2]) {
                                auto start = chrono::high_resolution_clock::now();
                                // since all points have been replaced, if it's different then midpoint has been created
                                // is it possible ?
                                auto hint = edges.erase(found.first);
                                edges.insert(hint, edge);
                                timeInsertErase += std::chrono::duration_cast<chrono::nanoseconds>( chrono::high_resolution_clock::now() - start).count();
                            }
                        }

                        end_quad = chrono::high_resolution_clock::now();

                        timeInsert += std::chrono::duration_cast<chrono::nanoseconds>(end_quad - start_quad).count();
                    }
                    auto end_quad = chrono::high_resolution_clock::now();

                    long total2 = std::chrono::duration_cast<chrono::nanoseconds>(end_quad - start_quad).count();
                    cout << " time edges " << total2 << " ns => init " << (timeInit * 100 / total2) << " %, insert " << (timeInsert * 100 / total2) << " % (" << (100 * timeFound / timeInsert) << "% try insert, " << (100 * timeInsertErase / timeInsert) << "% erase insert)" << endl;

                    //std::cout << "Edge end" << std::endl;
                });

                // Run another job concurrently with the loop above.
                // It can use up to the default number of threads.
                //tg.run([&] { // run in task group
                    //std::cout << "Quad start" << std::endl;

                auto start_quad = chrono::high_resolution_clock::now();
                    for (const Quadrant &local_quad : rmr.m_new_Quadrants) {
                        // build new quad with right index

                        vector<unsigned int> new_pointindex(4, 0);
                        for (unsigned int i = 0; i < 4; i++) {
                            if (local_quad.getPointIndex(i) < nb_points) {
                                // index refer point not created during this refinement level
                                new_pointindex[i] = local_quad.getPointIndex(i);
                            } else {
                                // point created, need to update the point with correct index
                                new_pointindex[i] = taskToGlobal[local_quad.getPointIndex(i)];
                            }
                        }

                        Quadrant quad(new_pointindex, m_rl);
                        m_new_Quadrants.push_back(quad);

                    }

                auto end_quad = chrono::high_resolution_clock::now();

                long total2 = std::chrono::duration_cast<chrono::milliseconds>(end_quad - start_quad).count();
                cout << " time quad " << total2 << endl;

                    //std::cout << "Quad end" << std::endl;
                //});

                // Wait for completion of the task group
                tg.wait();

                //std::cout << "End join" << (master ? " master" : "") << std::endl;

            } else {


                int i = 0;
                for (const Point3D &point : rmr.m_new_pts) {

                    //auto found = m_new_pts.find(point);
                    size_t hashPoint = point.operator()(point);

                    auto found = m_map_new_pts.find(hashPoint);
                    if (found != m_map_new_pts.end()) {
                        taskToGlobal[i++ + nb_points] = m_map_new_pts[hashPoint];
                    } else {
                        m_map_new_pts[hashPoint] = nb_points + m_new_pts.size();
                        m_new_pts.push_back(point);
                        taskToGlobal[i++ + nb_points] = nb_points + m_new_pts.size() - 1;
                    }
                }

                tbb::task_scheduler_init def_init; // Use the default number of threads.
                tbb::task_group tg;

                tg.run([&] { // run in task group
                    //std::cout << "Edge start" << std::endl;
                    for (const QuadEdge &local_edge : rmr.m_new_edges) {
                        // build new edge with right index
                        vector<unsigned long> index(3, 0);

                        for (unsigned int i = 0; i < 3; i++) {
                            if (local_edge[i] < nb_points) {
                                // index refer point not created during this refinement level
                                index[i] = local_edge[i];
                            } else {
                                // point created, need to update the point with correct index
                                index[i] = taskToGlobal[local_edge[i]];
                            }
                        }

                        QuadEdge edge(index[0], index[1], index[2]);

                        auto found = m_new_edges.find(edge);

                        if (found == m_new_edges.end()) {
                            m_new_edges.insert(edge);
                        } else {
                            if (edge[2] != 0 && edge[2] != (*found)[2]) {
                                // since all points have been replaced, if it's different then midpoint has been created
                                // is it possible ?
                                m_new_edges.erase(found);
                                m_new_edges.insert(edge);
                            }
                        }
                    }

                    //std::cout << "Edge end" << std::endl;
                });

                // Run another job concurrently with the loop above.
                // It can use up to the default number of threads.
                tg.run([&] { // run in task group
                    //std::cout << "Quad start" << std::endl;
                    for (const Quadrant &local_quad : rmr.m_new_Quadrants) {
                        // build new quad with right index

                        vector<unsigned int> new_pointindex(4, 0);
                        for (unsigned int i = 0; i < 4; i++) {
                            if (local_quad.getPointIndex(i) < nb_points) {
                                // index refer point not created during this refinement level
                                new_pointindex[i] = local_quad.getPointIndex(i);
                            } else {
                                // point created, need to update the point with correct index
                                new_pointindex[i] = taskToGlobal[local_quad.getPointIndex(i)];
                            }
                        }

                        Quadrant quad(new_pointindex, m_rl);
                        m_new_Quadrants.push_back(quad);

                    }
                    //std::cout << "Quad end" << std::endl;
                });

                // Wait for completion of the task group
                tg.wait();

                //std::cout << "End join" << (master ? " master" : "") << std::endl;
            }
        }

        /**
         * @brief Accumulate result for subrange.
         */
        void operator()(const tbb::blocked_range<size_t> &range) {

            list<RefinementRegion *>::const_iterator reg_iter;

            for (auto i = range.begin(); i != range.end(); ++i) {

                Quadrant &iter = tmp_Quadrants[i];

                //Only check, can not modify after treatment
                bool to_refine = false;

                for (reg_iter = all_reg.begin(), reg_iter++; reg_iter != all_reg.end(); ++reg_iter) {


                    unsigned short region_rl = (*reg_iter)->getRefinementLevel();
                    if (region_rl < m_rl) {
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
                    m_new_Quadrants.push_back(iter);
                    //End of for loop
                } else {
                    //(paul) Idea : add a task here (only if to refined, check if faster..)

                    list<unsigned int> &inter_edges = iter.getIntersectedEdges();
                    unsigned short qrl = iter.getRefinementLevel();

                    vector<vector<Point3D> > clipping_coords;
                    csv.setClipping(clipping_coords);

                    vector<vector<unsigned int> > split_elements;
                    csv.setNewEles(split_elements);

                    //auto start_sv_time = chrono::high_resolution_clock::now();
                    iter.accept(&csv);
                    //auto end_sv_time = chrono::high_resolution_clock::now();
                    //time_split_visitor += std::chrono::duration_cast<chrono::milliseconds>(end_sv_time - start_sv_time).count();

                    if (inter_edges.empty()) {
                        for (unsigned int j = 0; j < split_elements.size(); j++) {
                            Quadrant o(split_elements[j], qrl + 1);
                            m_new_Quadrants.push_back(o);
                        }
                    } else {
                        for (unsigned int j = 0; j < split_elements.size(); j++) {
                            Quadrant o(split_elements[j], qrl + 1);

                            //the new points are inserted in bash at the end of this
                            //iteration. For this reason, the coordinates must be passed
                            //"manually" at this point (clipping_coords).

                            //(paul) TODO add as attribute ?
                            IntersectionsVisitor iv(true);
                            iv.setPolyline(input);
                            iv.setEdges(inter_edges);
                            iv.setCoords(clipping_coords[j]);

                            if (o.accept(&iv)) {
                                m_new_Quadrants.push_back(o);
                            } else {
                                //The element doesn't intersect any input face.
                                //It must be checked if it's inside or outside.
                                //Only in the first case add it to m_new_Quadrants.
                                //Test this with parent Quadrant faces only.

                                //Comment the following lines of this 'else' if
                                //only intersecting Quadrants are meant to be
                                //displayed.

                                //note: inter_edges is quite enough to check if
                                //element is inside input, no Quadrant needed,
                                //so i moved the method to mesher  --setriva

                                if (isItIn(input, inter_edges, clipping_coords[j])) {
                                    m_new_Quadrants.push_back(o);
                                }
                            }
                        }
                    }
                }
            } //END FOR QUADRANTS
        }

        inline vector<Quadrant> &getNewQuadrants() { return m_new_Quadrants; }

        inline vector<Point3D> &getNewPts() { return m_new_pts; }

        inline set<QuadEdge> &getNewEdges() { return m_new_edges; }

        bool isItIn(const Polyline &mesh, const list<unsigned int> &faces, const vector<Point3D> &coords) const {
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

    };

}

//How to use it :

// https://www.threadingbuildingblocks.org/docs/help/tbb_userguide/parallel_reduce.html

// RefineMeshReduction rmr(a);
// parallel_reduce( blocked_range<size_t>(0,n), rmr );
// return rmr.my_sum;