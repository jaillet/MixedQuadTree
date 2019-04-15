#include "SplitVisitorNoCounterWithConcurrentSet.h"
#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_for.h>
#include <tbb/concurrent_vector.h>
#include <tbb/concurrent_unordered_set.h>
#include <tbb/parallel_reduce.h>


#include "../../Visitors/IntersectionsVisitor.h"
#include "../../Polyline.h"
#include "../../Quadrant.h"

using Clobscode::SplitVisitorNoCounterWithConcurrentSet;
using Clobscode::QuadEdge;
using Clobscode::Quadrant;
using Clobscode::MeshPoint;
using Clobscode::Polyline;
using Clobscode::Point3D;
using std::vector;
using std::list;
using std::set;


class RefineMeshReduction {
    
	//Private variables
    unsigned int m_rl;

    //list of temp Quadrants, only read
    vector<Quadrant> m_tmp_Quadrants;

    //Need to be read in reduction :
    vector<Quadrant> m_new_Quadrants;
    vector<Point3D> m_new_pts; //Local new pts (not merged)

    //Read only
    const Polyline *m_input;
    const vector<MeshPoint> m_points;

    //Only one shared, edges
    //Take ptr to share it without fear of non existing constructor
    tbb::concurrent_unordered_set<QuadEdge, std::hash<QuadEdge>>* m_quadEdges;

    SplitVisitorNoCounterWithConcurrentSet m_sv;

    void setSplitVisitor() {
		m_sv.setPoints(m_points);
		m_sv.setEdges(*m_quadEdges);
		m_sv.setNewPts(m_new_pts);
    }


public:

	//Output variables, filled in reduction
    vector<Point3D> new_pts; //Need to be filled from all rmr.m_new_pts

    //End output variables

    /**
     * @brief Splitting constructor. Must be able to run concurrently with operator() and method join.
     * @details split is a dummy argument of type split, distinguishes the splitting constructor from a copy constructor. 
     */
    RefineMeshReduction( RefineMeshReduction& x, split ) : 
    m_tmp_Quadrants(x.m_tmp_Quadrants) {

    	m_quadEdges = x.m_quadEdges;
    	setSplitVisitor();
    }
    
    /**
     * @brief Reduction.
     * @details Join results. The result in rmr should be merged into the result of this.
     */
    void join( const RefineMeshReduction& rmr ) {


    	//TODO :

    	//fill new_pts from all rmr.m_new_pts


    }
    
    RefineMeshReduction(unsigned int refinementLevel, const vector<Quadrant>& tmp_Quadrants, tbb::concurrent_unordered_set<QuadEdge, std::hash<QuadEdge>> *quadEdges, 
    					const Polyline *input, const vector<MeshPoint>& points) :
    m_rl(refinementLevel), m_tmp_Quadrants(tmp_Quadrants), m_input(input), m_points(points) {
    	m_quadEdges = quadEdges;
    	setSplitVisitor();
    }

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

    /**
     * @brief Accumulate result for subrange.
     */
    void operator()( const blocked_range<size_t>& range ) {

        nb_points = points.size();


        list<RefinementRegion *>::const_iterator reg_iter;

        for (auto j = range.begin(); j != range.end(); ++j) {

            Quadrant &iter = m_tmp_Quadrants[j];

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
                m_new_Quadrants.push_back(iter); //SHARED VARIABLE
                //End of for loop
            } else {
                //(paul) Idea : add a task here (only if to refined, check if faster..)

                list<unsigned int> &inter_edges = iter.getIntersectedEdges();
                unsigned short qrl = iter.getRefinementLevel();

                vector<vector<Point3D> > clipping_coords;
                sv.setClipping(clipping_coords);

                vector<vector<unsigned int> > split_elements;
                sv.setNewEles(split_elements);

                //auto start_sv_time = chrono::high_resolution_clock::now();
                iter.accept(&sv);
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

                            if (isItIn(*input, inter_edges, clipping_coords[j])) {
                                m_new_Quadrants.push_back(o);
                            }
                        }
                    }
                }
            }
        } //END FOR QUADRANTS    
    }
};


//How to use it :

// https://www.threadingbuildingblocks.org/docs/help/tbb_userguide/parallel_reduce.html

// RefineMeshReduction rmr(a);
// parallel_reduce( blocked_range<size_t>(0,n), rmr );
// return rmr.my_sum;