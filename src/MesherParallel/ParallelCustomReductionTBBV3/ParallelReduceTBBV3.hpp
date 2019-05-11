#include "CustomSplitVisitorV3.h"
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

    class RefineMeshReductionV3 {

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
        tbb::concurrent_vector<Quadrant> &tmp_Quadrants;
        tbb::concurrent_unordered_set<QuadEdge, std::hash<QuadEdge>> &edges;

        CustomSplitVisitorV3 csv;
        bool master;

        void setSplitVisitor() {
            csv.setPoints(points);
            csv.setEdges(edges);
            csv.setNewPts(m_new_pts);
            csv.setNewEdges(m_new_edges);
            csv.setMapPts(m_map_new_pts);
        }


    public:

        RefineMeshReductionV3(unsigned int refinementLevel, tbb::concurrent_vector<Quadrant> &tmp_Quadrants, tbb::concurrent_unordered_set<QuadEdge, std::hash<QuadEdge>> &quadEdges,
                              Polyline &input, vector<MeshPoint> &points, const list<RefinementRegion *> &all_reg,
                              const bool master) :
                m_rl(refinementLevel), input(input), points(points), all_reg(all_reg), tmp_Quadrants(tmp_Quadrants),
                edges(quadEdges), master(master) {
            setSplitVisitor();
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
                    //End of for loop
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

        inline std::tr1::unordered_map<size_t, unsigned int> &getNewMaps() { return m_map_new_pts; }

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