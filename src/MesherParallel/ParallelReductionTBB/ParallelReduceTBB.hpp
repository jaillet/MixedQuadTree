#include "../../ParallelizeTest/Join/CustomSplitVisitorV3.h"
#include <tbb/blocked_range.h>
#include <tr1/unordered_map>

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

    class RefineMeshReduction {

        //Private variables
        unsigned int m_rl;

        //Need to be read in reduction :
        vector<Quadrant> m_new_Quadrants;
        vector<Point3D> m_new_pts; //Local new pts (not merged)
        std::tr1::unordered_map<size_t, unsigned int> m_map_new_pts;
        set<QuadEdge> m_new_edges;

        //Read only
        const Polyline &input;
        const vector<MeshPoint> &points;
        const list<RefinementRegion *> &all_reg;
        vector<Quadrant> &tmp_Quadrants;
        set<QuadEdge> &edges;

        CustomSplitVisitor csv;

        int numberOfJoint;
        bool master;

        void setSplitVisitor() {
            csv.setPoints(points);
            csv.setEdges(edges);
            csv.setNewPts(m_new_pts);
            csv.setNewEdges(m_new_edges);
            csv.setMapPts(m_map_new_pts);
        }


    public:

        RefineMeshReduction(unsigned int refinementLevel, vector<Quadrant> &tmp_Quadrants, set<QuadEdge> &quadEdges,
                            Polyline &input, vector<MeshPoint> &points, const list<RefinementRegion *> &all_reg, const bool master) :
                m_rl(refinementLevel), input(input), points(points), all_reg(all_reg), tmp_Quadrants(tmp_Quadrants),
                edges(quadEdges), master(master){
            setSplitVisitor();
            numberOfJoint = 0;
        }

        /**
         * @brief Splitting constructor. Must be able to run concurrently with operator() and method join.
         * @details split is a dummy argument of type split, distinguishes the splitting constructor from a copy constructor.
         */
        RefineMeshReduction(RefineMeshReduction &x, tbb::split) :
                m_rl(x.m_rl), input(x.input), points(x.points), all_reg(x.all_reg), tmp_Quadrants(x.tmp_Quadrants),
                edges(x.edges) {
            setSplitVisitor();
            numberOfJoint = 0;
            master = false;
        }

        /**
         * @brief Reduction.
         * @details Join results. The result in rmr should be merged into the result of this.
         */
        void join(const RefineMeshReduction &rmr) {

            numberOfJoint += rmr.numberOfJoint + 1;

            // best than map for insert and access
            std::tr1::unordered_map<unsigned int, unsigned int> taskToGlobal;

            int i = 0;
            for (const Point3D &point : rmr.m_new_pts) {

                //auto found = m_new_pts.find(point);
                size_t hashPoint = point.operator()(point);

                auto found = m_map_new_pts.find(hashPoint);
                if (found != m_map_new_pts.end()) {
                    taskToGlobal[i++ + points.size()] = m_map_new_pts[hashPoint];
                } else {
                    m_map_new_pts[hashPoint] = points.size() + m_new_pts.size();
                    m_new_pts.push_back(point);
                    taskToGlobal[i++ + points.size()] = points.size() + m_new_pts.size() - 1;
                }
            }

            tbb::task_scheduler_init def_init; // Use the default number of threads.
            tbb::task_group tg;

            tg.run([&] { // run in task group
                for (const QuadEdge &local_edge : rmr.m_new_edges) {
                    // build new edge with right index
                    vector<unsigned long> index(3, 0);

                    for (unsigned int i = 0; i < 3; i++) {
                        if (local_edge[i] < points.size()) {
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
            });

            // Run another job concurrently with the loop above.
            // It can use up to the default number of threads.
            tg.run([&] { // run in task group
                for (const Quadrant &local_quad : rmr.m_new_Quadrants) {
                    // build new quad with right index

                    vector<unsigned int> new_pointindex(4, 0);
                    for (unsigned int i = 0; i < 4; i++) {
                        if (local_quad.getPointIndex(i) < points.size()) {
                            // index refer point not created during this refinement level
                            new_pointindex[i] = local_quad.getPointIndex(i);
                        } else {
                            // point created, need to update the point with correct index
                            new_pointindex[i] = taskToGlobal[local_quad.getPointIndex(i)];
                        }
                    }

                    Quadrant quad(new_pointindex, m_rl);
                    quad.intersected_edges = local_quad.intersected_edges;
                    quad.intersected_features = local_quad.intersected_features;
                    m_new_Quadrants.push_back(quad);

                }
            });

            // Wait for completion of the task group
            tg.wait();
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
                        break;
                    }
                }


                //now if refinement is not needed, we add the Quadrant as it was.
                if (!to_refine) {
                    m_new_Quadrants.push_back(iter);
                } else {
                    //(paul) Idea : add a task here (only if to refined, check if faster..)

                    list<unsigned int> &inter_edges = iter.getIntersectedEdges();
                    unsigned short qrl = iter.getRefinementLevel();

                    vector<vector<Point3D> > clipping_coords;
                    csv.setClipping(clipping_coords);

                    vector<vector<unsigned int> > split_elements;
                    csv.setNewEles(split_elements);

                    iter.accept(&csv);

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