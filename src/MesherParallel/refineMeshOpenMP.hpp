#include "../Mesher.h"
#include <omp.h>

#include "SplitVisitorOpenMP.h"

using namespace std;

namespace Clobscode {

	void Mesher::refineMeshParallelOpenMP(int nbThread, list<Quadrant> Quadrants, vector<MeshPoint> points,
                                        set<QuadEdge> QuadEdges,
                                        const list<RefinementRegion *> &all_reg, const unsigned short &rl,
                                        Polyline &input) {
            int MAX_THREAD = omp_get_max_threads();

            if (nbThread > MAX_THREAD || nbThread < 0) {
                std::cout << "Invalid number of threads or not supported by computer" << std::endl;
                return;
            }
            omp_set_num_threads(nbThread);


            std::cout << "Start refine mesh parallel Open MP with " << nbThread << " threads." << std::endl;

            //Convert list into vector of temp Quadrants
            vector<Quadrant> tmp_Quadrants;
            tmp_Quadrants.assign(make_move_iterator(Quadrants.begin()), make_move_iterator(Quadrants.end()));


            //Shared variables
            vector<Quadrant> new_Quadrants;
            vector<Point3D> new_pts;
            //Also QuadEdges

            tbb::atomic<int> counterRefine = 0;

            unsigned int nb_points;

            for (unsigned short i = 0; i < rl; i++) {
                auto start_refine_rl_time = chrono::high_resolution_clock::now();

                //atomic<long> time_split_visitor;
                //atomic<long> time_inside_block;
                //time_inside_block = 0;
                //time_split_visitor = 0;

                //the new_pts is a list that holds the coordinates of
                //new points inserted at this iteration. At the end of
                //this bucle, they are inserted in the point vector
                new_pts.clear();

                nb_points = points.size();

                //auto start_outside_block_time = chrono::high_resolution_clock::now();

                //split the Quadrants as needed
				  //auto start_block_time = chrono::high_resolution_clock::now();


                

			  #pragma omp parallel shared(nb_points, new_Quadrants, new_pts, QuadEdges, points)
			  	{
				//create visitors and give them variables
				SplitVisitorOpenMP sv(nb_points);
				sv.setPoints(points); // READ
				sv.setEdges(QuadEdges); // INSERT / REMOVE / READ
				sv.setNewPts(new_pts); // INSERT / READ

				//schedule(dynamic)
				#pragma omp for 
			  	for (unsigned int j = 0; j < tmp_Quadrants.size(); ++j) {
			  		Quadrant &iter = tmp_Quadrants[j];

					//Only check, can not modify after treatment
					bool to_refine = false;
					list<RefinementRegion *>::const_iterator reg_iter;

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
			      		#pragma omp critical(new_quad)
			        	new_Quadrants.push_back(iter);
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
			                  #pragma omp critical(new_quad)
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
			                  		#pragma omp critical(new_quad)
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
			                      		#pragma omp critical(new_quad)
			                        	new_Quadrants.push_back(o);
			                      }
			                  }
			              }
			          }
			      }
			  	} //END FOR QUADRANTS


				  //auto end_block_time = chrono::high_resolution_clock::now();
				  //time_inside_block += std::chrono::duration_cast<chrono::milliseconds>(end_block_time - start_block_time).count();
				  //cout << time << " / " << std::chrono::duration_cast<chrono::milliseconds>(end_block_time - start_block_time).count();
				  //cout << " ms for " << range.size() << " quadrants" << endl;

			  	}//END PARALLEL

                //auto end_outside_block_time = chrono::high_resolution_clock::now();

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
                long total = std::chrono::duration_cast<chrono::milliseconds>(end_refine_rl_time - start_refine_rl_time).count();
                
             	cout << "         * level " << i << " in "
                     << total;
                cout << " ms" << endl;

                std::cout << "           ---- Points : " << points.size() << std::endl;
	            std::cout << "           ---- QuadEdge : " <<  QuadEdges.size() << std::endl;
	            std::cout << "           ---- Quadrants : " << tmp_Quadrants.size() << std::endl; 

                // inside / split / outside stats

                //long outside = std::chrono::duration_cast<chrono::milliseconds>(end_outside_block_time - start_outside_block_time).count();
                //cout << "TBB for outside " << outside << " ms (" << (outside * 100.0 / total) << "%) ";
                //cout << "TBB for inside " << time_inside_block << " ms (all threads cumulated) ";
                //cout << " split visitor " << time_split_visitor << " ms (" << (time_split_visitor * 100.0 / time_inside_block) << "% of time inside) ";
                //cout << endl;




	            
            } //END FOR REFINEMENT LEVEL

            // output result to mesher ! comment if not needed
            //Quadrants.clear();
            //Quadrants.assign(tmp_Quadrants.begin(), tmp_Quadrants.end());

            //QuadEdges.clear();
            //QuadEdges.insert(quadEdges.begin(), quadEdges.end());

            std::cout << "----------------------------------------------------------" << std::endl;
            std::cout << "----------------------------------------------------------" << std::endl;
            std::cout << "---------------------END OF OPENMP------------------------" << std::endl;
            std::cout << "----------------------------------------------------------" << std::endl;
            std::cout << "----------------------------------------------------------" << std::endl;
        }

        /*
        void Mesher::refineMeshReductionOpenMP(int nbThread, list<Quadrant> Quadrants, vector<MeshPoint> points,
                                        set<QuadEdge> QuadEdges,
                                        const list<RefinementRegion *> &all_reg, const unsigned short &rl,
                                        Polyline &input) {
            int MAX_THREAD = omp_get_max_threads();

            if (nbThread > MAX_THREAD || nbThread < 0) {
                std::cout << "Invalid number of threads or not supported by computer" << std::endl;
                return;
            }
            omp_set_num_threads(nbThread);


            std::cout << "Start refine mesh reduction Open MP with " << nbThread << " threads." << std::endl;

            //Convert list into vector of temp Quadrants
            vector<Quadrant> tmp_Quadrants;
            tmp_Quadrants.assign(make_move_iterator(Quadrants.begin()), make_move_iterator(Quadrants.end()));


            //Shared variables
            vector<Quadrant> new_Quadrants;
            vector<Point3D> new_pts;
            //Also QuadEdges

            tbb::atomic<int> counterRefine = 0;

            unsigned int nb_points;

            for (unsigned short i = 0; i < rl; i++) {
                auto start_refine_rl_time = chrono::high_resolution_clock::now();

                //atomic<long> time_split_visitor;
                //atomic<long> time_inside_block;
                //time_inside_block = 0;
                //time_split_visitor = 0;

                //the new_pts is a list that holds the coordinates of
                //new points inserted at this iteration. At the end of
                //this bucle, they are inserted in the point vector
                new_pts.clear();

                nb_points = points.size();

                //auto start_outside_block_time = chrono::high_resolution_clock::now();

                //split the Quadrants as needed
				  //auto start_block_time = chrono::high_resolution_clock::now();


                //create visitors and give them variables
				SplitVisitorOpenMP sv(nb_points);
				sv.setPoints(points); // READ
				sv.setEdges(QuadEdges); // INSERT / REMOVE / READ
				sv.setNewPts(new_pts); // INSERT / READ


			  #pragma omp parallel for private(sv) shared(nb_points, new_Quadrants, new_pts, QuadEdges)
			  	for (unsigned int j = 0; j < tmp_Quadrants.size(); ++j) {
			  		Quadrant &iter = tmp_Quadrants[j];

					//Only check, can not modify after treatment
					bool to_refine = false;
					list<RefinementRegion *>::const_iterator reg_iter;

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
			      		#pragma omp critical(new_quad)
			        	new_Quadrants.push_back(iter);
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
			                  #pragma omp critical(new_quad)
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
			                  		#pragma omp critical(new_quad)
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
			                      		#pragma omp critical(new_quad)
			                        	new_Quadrants.push_back(o);
			                      }
			                  }
			              }
			          }
			      }
			  } //END FOR QUADRANTS


				  //auto end_block_time = chrono::high_resolution_clock::now();
				  //time_inside_block += std::chrono::duration_cast<chrono::milliseconds>(end_block_time - start_block_time).count();
				  //cout << time << " / " << std::chrono::duration_cast<chrono::milliseconds>(end_block_time - start_block_time).count();
				  //cout << " ms for " << range.size() << " quadrants" << endl;


                //auto end_outside_block_time = chrono::high_resolution_clock::now();

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
                long total = std::chrono::duration_cast<chrono::milliseconds>(end_refine_rl_time - start_refine_rl_time).count();
                
             	cout << "         * level " << i << " in "
                     << total;
                cout << " ms" << endl;

                std::cout << "           ---- Points : " << points.size() << std::endl;
	            std::cout << "           ---- QuadEdge : " <<  QuadEdges.size() << std::endl;
	            std::cout << "           ---- Quadrants : " << tmp_Quadrants.size() << std::endl; 

                // inside / split / outside stats

                //long outside = std::chrono::duration_cast<chrono::milliseconds>(end_outside_block_time - start_outside_block_time).count();
                //cout << "TBB for outside " << outside << " ms (" << (outside * 100.0 / total) << "%) ";
                //cout << "TBB for inside " << time_inside_block << " ms (all threads cumulated) ";
                //cout << " split visitor " << time_split_visitor << " ms (" << (time_split_visitor * 100.0 / time_inside_block) << "% of time inside) ";
                //cout << endl;

	            
            } //END FOR REFINEMENT LEVEL

            // output result to mesher ! comment if not needed
            //Quadrants.clear();
            //Quadrants.assign(tmp_Quadrants.begin(), tmp_Quadrants.end());

            //QuadEdges.clear();
            //QuadEdges.insert(quadEdges.begin(), quadEdges.end());

            std::cout << "----------------------------------------------------------" << std::endl;
            std::cout << "----------------------------------------------------------" << std::endl;
            std::cout << "---------------------END OF OPENMP------------------------" << std::endl;
            std::cout << "----------------------------------------------------------" << std::endl;
            std::cout << "----------------------------------------------------------" << std::endl;
        }

        */


}