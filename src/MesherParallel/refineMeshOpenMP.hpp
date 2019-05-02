#include "../Mesher.h"
#include <omp.h>

#include <iostream>
#include "SplitVisitorOpenMP.h"
#include "SplitVisitorReductionOpenMP.h"

#include <unordered_map>

using namespace std;

namespace Clobscode {

	void make_reduceV1(	vector<Point3D> *threads_new_pts, set<QuadEdge> *threads_new_edges, vector<Quadrant> *threads_new_quadrants, 
					unsigned int number_of_threads, const unsigned int total_nb_points_before_reduce, unsigned int rl);


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

        
        void Mesher::refineMeshReductionOpenMP(const int nbThread, list<Quadrant> Quadrants, vector<MeshPoint> points,
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


            //Each thread has it's own data set
            vector<Quadrant> thread_new_Quadrants[nbThread];
            vector<Point3D> thread_new_pts[nbThread];
            set<QuadEdge> thread_new_edges[nbThread];

            //Result of reduce :
            vector<Quadrant> reduced_new_Quadrants;


            unsigned int nb_points;

            for (unsigned short i = 0; i < 1; i++) {
                auto start_refine_rl_time = chrono::high_resolution_clock::now();


                nb_points = points.size();

				//Begin of parallel region
				#pragma omp parallel shared(nb_points, thread_new_Quadrants, thread_new_pts, thread_new_edges)
				{

					//std::cout << "THREAD NUMBER : " <<  << std::endl;

					//Create reference to the dataset of actual thread
					vector<Quadrant> &new_Quadrants = thread_new_Quadrants[omp_get_thread_num()];
					vector<Point3D> &new_pts = thread_new_pts[omp_get_thread_num()];
					set<QuadEdge> &new_edges = thread_new_edges[omp_get_thread_num()];

					new_Quadrants.clear();
					new_pts.clear();
					new_edges.clear();

					//create visitors and give them variables
					//Each thread link its own dataset
					SplitVisitorReductionOpenMP sv;
					sv.setPoints(points); // READ only
					//sv.setEdges(new_edges); // INSERT / REMOVE / READ
					//Split the method to have a dataset for reading only (current_edges)
					//And one only for new edges -> empty at the beginning, filled by sv
					sv.setCurrentEdges(QuadEdges);
					sv.setNewEdges(new_edges);

					sv.setNewPts(new_pts); // INSERT / READ
					sv.setThreadNum(omp_get_thread_num());

					//split the Quadrants as needed
			  		#pragma omp parallel for schedule(dynamic)
			  		for (unsigned int j = 0; j < tmp_Quadrants.size(); ++j) {
			  			Quadrant &iter = tmp_Quadrants[j];

			  		//Check if to refine
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
								break;
							}
						}

					//END check if to refine


						//now if refinement is not needed, we add the Quadrant as it was.
				      	if (!to_refine) {
				      		
				        	new_Quadrants.push_back(iter);
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

									} // END ELSE ACCEPT
								} // END FOR SPLIT ELEMENTS
							} //END ELSE INTER_EDGES EMPTY
						} //END ELSE TO REFINE
					} //END FOR QUADRANTS
				

					// Begin reduce


					//Wait all thread
					#pragma omp barrier

					//Only one thread does the reduce
					if (omp_get_thread_num() == 0) {
						//Master thread, use #pragma omp master ?

						// Fill index 0 in arrays thread_new_pts, new_edges, new_Quadrants
						// with a size=nbThread the result of the merge of all thread's new elements
						// points.size() is needed to know at which index threads have started
						make_reduceV1(thread_new_pts, thread_new_edges, thread_new_Quadrants,
									nbThread, points.size(), i);

						// Output of reduce : 
						//thread_new_pts[0]
						//thread_new_edges[0]
						//thread_new_Quadrants[0]


					}

					// End reduce


				} // END PARALLEL REGION



				//Update edges
				for (set<QuadEdge>::iterator i = thread_new_edges[0].begin(); i != thread_new_edges[0].end(); ++i)
				{
					QuadEdges.insert(*i);
				}


				// don't forget to update list
                tmp_Quadrants.assign(make_move_iterator(thread_new_Quadrants[0].begin()), make_move_iterator(thread_new_Quadrants[0].end()));
                thread_new_Quadrants[0].clear();

                //if no points were added at this iteration, it is no longer
                //necessary to continue the refinement.
                if (thread_new_pts[0].empty()) {
                    cout << "warning at Mesher::generateQuadtreeMesh no new points!!!\n";
                    break;
                }

                //add the new points to the vector
                points.reserve(points.size() + thread_new_pts[0].size());
                points.insert(points.end(), thread_new_pts[0].begin(), thread_new_pts[0].end());

                auto end_refine_rl_time = chrono::high_resolution_clock::now();
                long total = std::chrono::duration_cast<chrono::milliseconds>(end_refine_rl_time - start_refine_rl_time).count();
                
             	cout << "         * level " << i << " in " << total << " ms" << endl;

                cout << "           ---- Points : " << points.size() << endl;
	            cout << "           ---- QuadEdge : " <<  QuadEdges.size() << endl;
	            cout << "           ---- Quadrants : " << tmp_Quadrants.size() << endl; 

			} // END FOR REFINEMENT LEVEL


            // output result to mesher ! comment if not needed
            //Quadrants.clear();
            //Quadrants.assign(tmp_Quadrants.begin(), tmp_Quadrants.end());

            //QuadEdges.clear();
            //QuadEdges.insert(quadEdges.begin(), quadEdges.end());

            std::cout << "----------------------------------------------------------" << std::endl;
            std::cout << "----------------------------------------------------------" << std::endl;
            std::cout << "---------------------END OF OPENMP REDUCTION------------------------" << std::endl;
            std::cout << "----------------------------------------------------------" << std::endl;
            std::cout << "----------------------------------------------------------" << std::endl;
        }


	void make_reduceV1(	vector<Point3D> *threads_new_pts, set<QuadEdge> *threads_new_edges, vector<Quadrant> *threads_new_quadrants, 
						unsigned int number_of_threads, const unsigned int total_nb_points_before_reduce, unsigned int rl) {



		std::cout << "BEGIN REDUCE at level rl=" << rl << std::endl;

		//The result :
		vector<Point3D> &new_pts = threads_new_pts[0];
		set<QuadEdge> &new_edges = threads_new_edges[0];
		vector<Quadrant> &new_quadrants = threads_new_quadrants[0];

		//The current number of points
		unsigned int total_nb_points = total_nb_points_before_reduce;
		// The map to link real index of new_pts
		unordered_map<unsigned int, unsigned int> map_new_pts;
		// The task to global
		unordered_map<unsigned int, unsigned int> taskToGlobal;
		for (int thread_number = 1; thread_number < number_of_threads; ++thread_number)
		{
			

			std::cout << "Jointure thread number i=" << thread_number << std::endl;

			std::cout << "Points start" << std::endl;

			vector<Point3D> &thread_new_pts = threads_new_pts[thread_number];
			for(const Point3D& point : thread_new_pts) {
				size_t hashPoint = point.operator()(point);
		        auto found = map_new_pts.find(hashPoint);

		        if (found == map_new_pts.end()) {
		            //We did not found the point
		            //We add it in the global new_pts
		            new_pts.push_back(point);
		            //We add
		            map_new_pts[hashPoint] = total_nb_points;
		            
		        }
		        //Else we found it
		        taskToGlobal[total_nb_points] = map_new_pts[hashPoint];
		        total_nb_points++;

			}

			std::cout << "Points end" << std::endl;

			std::cout << "Edge start" << std::endl;

			set<QuadEdge> &thread_new_edges = threads_new_edges[thread_number];
			
			for(const QuadEdge& local_edge : thread_new_edges)
			{

				std::cout << "Edge : " << local_edge << std::endl;

				// build new edge with right index
		        vector<unsigned long> index(3, 0);

		        bool all_created = true;
		        for (unsigned int j = 0; j < 3; j++) {
		            if (local_edge[j] < total_nb_points_before_reduce || local_edge[j] == 0) {
		            	//Don't forget midpoint = 0
		                // index refer point not created during this refinement level
		                index[j] = local_edge[j];
		            } else {
		                // point created, need to update the point with correct index
		                index[j] = taskToGlobal[local_edge[j]];
		                all_created = false;
		            }
		        }

		        QuadEdge edge(index[0], index[1], index[2]);


		        //If all created, no need to search, we can insert directly
		        if (all_created) {
		        	new_edges.insert(edge);
		        }
		        else {
			        auto found = new_edges.find(edge);
			        if (found == new_edges.end()) {
			            new_edges.insert(edge);
			        } else {
			            if (edge[2] != 0 && edge[2] != (*found)[2]) {
			                // since all points have been replaced, if it's different then midpoint has been created
			                // is it possible ?
			                new_edges.erase(found);
			                new_edges.insert(edge);
			            }
			        }
			    }
			}

			std::cout << "Edge end" << std::endl;

			std::cout << "Quad start" << std::endl;
			vector<Quadrant> &thread_new_quadrants = threads_new_quadrants[thread_number];
			
			for(const Quadrant &local_quad : thread_new_quadrants) {
				// build new quad with right index
		        vector<unsigned int> new_pointindex(4, 0);
		        for (unsigned int j = 0; j < 4; j++) {
		            if (local_quad.getPointIndex(j) < total_nb_points_before_reduce) {
		                // index refer point not created during this refinement level
		                new_pointindex[j] = local_quad.getPointIndex(j);
		            } else {
		                // point created, need to update the point with correct index
		                new_pointindex[j] = taskToGlobal[local_quad.getPointIndex(j)];
		            }
		        }

		        Quadrant quad(new_pointindex, rl);
		        new_quadrants.push_back(quad);
			
			}
	    	std::cout << "Quad end" << std::endl;


		} // END FOR ALL THREAD

		std::cout << "END REDUCE at level rl=" << rl << std::endl;
		std::cout << "---------------------------------------" << std::endl;


	}

}