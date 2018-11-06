/*
 <Mix-mesher: region type. This program generates a mixed-elements 2D mesh>

 Copyright (C) <2013,2018>  <Claudio Lobos> All rights reserved.

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/lgpl.txt>
 */
/**
* @file Main.cpp
* @author Claudio Lobos, Fabrice Jaillet
* @version 0.1
* @brief
**/

#include "Mesher.h"
#include "Polyline.h"
#include "FEMesh.h"
#include "Services.h"
#include "RefinementCubeRegion.h"
#include "RefinementSurfaceRegion.h"
#include "RefinementInputSurfaceRegion.h"
#include "RefinementBoundaryRegion.h"
#include "RefinementAllRegion.h"
#include "Point3D.h"
#include <string>
#include <cctype>
#include <time.h>
#include <chrono>

using std::atoi;
using std::cout;
using std::cerr;
using std::vector;
using std::string;
using std::endl;
using Clobscode::RefinementRegion;
using Clobscode::RefinementCubeRegion;
using Clobscode::RefinementSurfaceRegion;
using Clobscode::RefinementBoundaryRegion;
using Clobscode::RefinementAllRegion;
using Clobscode::RefinementInputSurfaceRegion;
using Clobscode::Point3D;
using Clobscode::Services;

namespace chrono = std::chrono;

//-------------------------------------------------------------------
//-------------------------------------------------------------------

void endMsg(){
    cout << std::flush;
    cerr << "use: ./mesher [-p] input.poly [-u] output\n";
    cerr << "              [-s] ref_level [-a] ref_level [-b] file.reg\n";
    cerr << "              [-r] input_surface rl [-g] [-v]\n";
    cerr << "where:\n";
    cerr << "  one of the parameters must be an input surface mesh in\n";
    cerr << "  mdl or off format. If output name is not provided it\n";
    cerr << "  will be saved in input_name.m3d. Options:\n";
    cerr << "    -s Refine Quadrants intersecting the input surface.\n";
    cerr << "       Parameter ref_level is the refinement level\n";
    cerr << "    -a Refine all elements in the input domain.\n";
    cerr << "       Parameter ref_level is the refinement level\n";
    cerr << "    -b Refine block regions provided in file file.reg\n";
    cerr << "    -r Refine surface region. Will refine all the elements\n";
    cerr << "       in the provided input_surface at level rl\n";
    cerr << "    -q if supported (only VTK by now), write quality attributes to output file.\n";
    cerr << "    -g save output mesh in GetFem format (gmf) \n";
    cerr << "    -v save output mesh + input in VTK ASCII format (vtk)\n";
    cerr << "    -m save output mesh in M3D ASCII format (m3d)\n";
    cerr << "    -x save output mesh in GMSH ASCII format (gmsh)\n";
    cerr << "    -i save output mesh in MVM ASCII format (mvm)\n";
    cerr << "    -o save output mesh in OFF ASCII format (off)\n";
}

//-------------------------------------------------------------------
//-------------------------------------------------------------------

int main(int argc,char** argv){
	
    if (argc<4) {
        endMsg();
        return EXIT_FAILURE;
    }
    
	//const int n_meshes = 1;
    string in_name, out_name;
	bool out_name_given = false, in_name_given = false;
//	bool edge_projection = false;
	
	unsigned short ref_level = 0, rl = 0, cminrl=0, omaxrl=0;
    //cminrl: current min refinement level (used when starting from an Quadtree mesh)
    //omaxrl: old max refinement level (used when starting from an Quadtree mesh)
    list<unsigned int> roctli;
    //this list contains the index of the Quadrants previously generated that need
    //one extra level of refinement.
    
	list<RefinementRegion *> all_regions;
    RefinementRegion *rr;

    vector<double> bounds;
    Point3D pmin,pmax;
    
    vector<Clobscode::Polyline> inputs;
//    inputs.reserve(4);
    //Clobscode::Services io;
    
    bool getfem=false, gmshformat=false, vtkformat=false, Quadrant_start=false, m3dfor=false;
    bool mvmfor=false, offfor=false;
    bool decoration=false; //if supported, write quality attributes to output file
    
    //for reading an Quadrant mesh as starting point.
    vector<MeshPoint> oct_points;
    vector<Quadrant> oct_Quadrants;
    set<QuadEdge> oct_edges;
    vector<unsigned int> oct_ele_link;
    GeometricTransform gt;
    
    for (int i=1; i<argc; i++) {
        
		if (argv[i][0]!='-') {
			cout << "Error: expected option -X and got " << argv[i] << "\n";
			endMsg();
            return EXIT_FAILURE;
		}
        
        bool inout = false;
        switch (argv[i][1]) {
            case 'g':
            case 'v':
            case 'm':
            case 'x':
            case 'i':
            case 'o':
            case 'q':
                inout = true;
                break;
            default:
                break;
        }
        
        if (argc==i+1 && !inout) {
            std::cerr << "Error: expected argument for option " << argv[i] << "\n" << std::flush;
			endMsg();
            exit (EXIT_FAILURE);
		}
        switch (argv[i][1]) {
            case 't': //test polyline
                for (uint i=0;i<inputs.size();++i) {
                    std::cerr << inputs[i].crossingNumber(Point3D(.5,0.5,0.0));
                    std::cerr << inputs[i].crossingNumber(Point3D(1.0,0.0,0.0));
                    std::cerr << inputs[i].crossingNumber(Point3D(0.42,0.49,0.0));
                    std::cerr << inputs[i].crossingNumber(Point3D(1.0000000001,0.0,0.0));
                    std::cerr << inputs[i].crossingNumber(Point3D(.99999999999,0.0,0.0)) << std::endl;
                    std::cerr << inputs[i].windingNumber(Point3D(0.5,0.5,1.0));
                    std::cerr << inputs[i].windingNumber(Point3D(0.2,1.2,1.0));
                    std::cerr << inputs[i].windingNumber(Point3D(1.2,0.2,1.0));
                    std::cerr << inputs[i].windingNumber(Point3D(0.4,0.3,1.0));
                    std::cerr << inputs[i].windingNumber(Point3D(.9999999999999999,1.0,0.0));
                    std::cerr << inputs[i].windingNumber(Point3D(1.0000000001,0.000000001,0.0));
                    std::cerr << inputs[i].windingNumber(Point3D(0.9999999999,0.999999999,0.0)) << std::endl;
                    std::cerr << inputs[i].getProjection(Point3D(0.5,0.50000001,0.0))<< std::endl;
                    std::cerr << inputs[i].getProjection(Point3D(0.2,1.2,0.0))<< std::endl;
                    std::cerr << inputs[i].getProjection(Point3D(1.2,0.2,0.0))<< std::endl;
                    std::cerr << inputs[i].getProjection(Point3D(1.2,0.2,1.0))<< std::endl;
                    std::cerr << inputs[i].getProjection(Point3D(.9999999999999999,1.0,0.0))<< std::endl;
                    std::cerr << inputs[i].getProjection(Point3D(1.0000000001,0.000000001,0.0))<< std::endl;
                    std::cerr << inputs[i].getProjection(Point3D(0.9999999999,0.999999999,0.0)) << std::endl;
                    std::cerr << inputs[i].getCentroid() << std::endl;
                }
                break;
            case 'd':
                in_name = argv[i+1];

                if (!Services::ReadMdlMesh(in_name,inputs)) {
                    std::cerr << "couldn't read file " << argv[i+1] << std::endl;
                    return EXIT_FAILURE;
                }
                in_name_given = true;
                i++;
                break;
            case 'p': //read a .poly file (Triangle)
                in_name = argv[i+1];

                if (!Services::ReadPolyFile(in_name,inputs)) {
                    std::cerr << "couldn't read file " << argv[i+1] << std::endl;
                    return EXIT_FAILURE;
                }
                in_name_given = true;
                i++;
                break;
            case 'u':
                out_name = argv[i+1];
                out_name_given = true;
                i++;
                break;
            case 'g':
                getfem = true;
                break;
            case 'x':
                gmshformat = true;
                break;
            case 'v':
                vtkformat = true;
                break;
            case 'q':
                decoration = true;
                break;
            case 'm':
                m3dfor = true;
                break;
            case 'i':
                mvmfor = true;
                break;
            case 'o':
                offfor = true;
                break;
            case 'a':
                rl = atoi(argv[i+1]);
                if (ref_level<rl) {
                    ref_level = rl;
                }
                //+-10 is an arbitrary number to ensure the Bbox contains
                //the entire input mesh
                rr = new RefinementAllRegion(rl);
                
                //see if force rotation enable
                if (argv[i][2]=='r') {
                    rr->forceInputRotation();
                }
                
                all_regions.push_back(rr);
                i++;
                break;
            case 's':
                rl = atoi(argv[i+1]);
                if (ref_level<rl) {
                    ref_level = rl;
                }
                rr = new RefinementInputSurfaceRegion(rl);
                
                //see if force rotation enable
                if (argv[i][2]=='r') {
                    rr->forceInputRotation();
                }
                
                all_regions.push_back(rr);
                i++;
                break;
            case 'b':
                unsigned short max_reg;
                max_reg= Services::readRefinementRegions(argv[i+1],all_regions);
                if (ref_level<max_reg) {
                    ref_level = max_reg;
                }
                i++;
                break;
            case 'r':
                rl = atoi(argv[i+2]);
                Services::readSurfaceRefinementRegion(argv[i+1],all_regions,rl);
                if (ref_level<rl) {
                    ref_level = rl;
                }
                i+=2;
                break;
            case 'c':
                Quadrant_start = true;
                Services::ReadQuadMesh(argv[i+1], oct_points, oct_Quadrants,
                                         oct_edges,oct_ele_link,gt,cminrl,omaxrl);
                if (ref_level<omaxrl) {
                    ref_level = omaxrl;
                }
                i++;
                break;
            case 'l':
                if (Quadrant_start) {
                    Services::ReadQuadrantList(argv[i+1],roctli,oct_ele_link);
                    list<unsigned int>::iterator oeiter;
                    for (oeiter=roctli.begin(); oeiter!=roctli.end(); oeiter++) {
                        rl = oct_Quadrants[*oeiter].getRefinementLevel();
                        if (ref_level<=rl) {
                            ref_level = rl+1;
                        }
                    }
                }
                else {
                    cerr << "Warning: option -l needs a previously provided Quadrant";
                    cerr << " mesh (option -o) skipping\n";
                }
                i++;
                break;
            default:
                cerr << "Warning: unknown option " << argv[i] << " skipping\n";
                break;
        }
    }
    
    //cout << "max target rl: " << ref_level << endl;
    
    if (!in_name_given) {
        cerr << "No input domain surface mesh provided. Aborting\n";
        list<RefinementRegion *>::iterator rriter;
        for (rriter = all_regions.begin(); rriter!=all_regions.end(); rriter++) {
            delete *rriter;
        }
        exit (EXIT_FAILURE);
    }
	
	//give default output name if non is provided
	if (!out_name_given) {
		unsigned int last_point = in_name.find_last_of(".");
		out_name = in_name.substr(0,last_point);
	}

    // Uncomment this if you want to export input as .poly file
    //Services::WritePolyFile(out_name,inputs);
	
    cout << "  Starting generation/refinement\n";
    auto start_time = chrono::high_resolution_clock::now();
    
    //Generate the mesh following the above constraints.
	Clobscode::Mesher mesher;
    std::shared_ptr<Clobscode::FEMesh> output;
    
    // and next proceed with mesh generation or refinement
    if (!Quadrant_start) {

        //Boundary: add in first position a Region dedicated to handle correctly the boundary
        //see if force rotation enable
        RefinementBoundaryRegion *rb =new RefinementBoundaryRegion(inputs.at(0),0);//std::numeric_limits<unsigned short>::max());
        //    if (argv[i][2]=='r') {
        //        rb->forceInputRotation();
        //    }
        all_regions.push_front(rb);

        output = mesher.generateMesh(inputs.at(0),ref_level,out_name,all_regions,decoration);
    }
    else {
        mesher.setInitialState(oct_points,oct_Quadrants,oct_edges);
        if (omaxrl<ref_level) {
            omaxrl = ref_level;
        }
        output = mesher.refineMesh(inputs.at(0),ref_level,out_name,roctli,
                                   all_regions,gt,cminrl,omaxrl,decoration);
    }
    
    auto gen_time = chrono::high_resolution_clock::now();
    cout << "  Generation done in " << std::chrono::duration_cast<chrono::milliseconds>(gen_time-start_time).count();
    cout << " ms"<< endl;
    
    if (getfem) {
        Services::WriteMeshGetfem(out_name,output);
    }

    if (gmshformat) {
        Services::WriteGMSH(out_name,output);
    }

    if (vtkformat) {
        Services::WriteVTK(out_name,output);
        Services::WritePolyVTK(out_name,inputs);
    }

    if (m3dfor) {
        Services::WriteOutputMesh(out_name,output);
    }
    if (mvmfor) {
        Services::WriteMixedVolumeMesh(out_name,output);
    }
    if (offfor) {
        Services::WriteOFF(out_name,output);
    }

    auto end_time = chrono::high_resolution_clock::now();
    cout << "  Write done in " << std::chrono::duration_cast<chrono::milliseconds>(end_time-gen_time).count();
    cout << " ms"<< endl;
    cout << "  All done in " << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count();
    cout << " ms"<< endl;
	
    list<RefinementRegion *>::iterator rriter;
    for (rriter = all_regions.begin(); rriter!=all_regions.end(); rriter++) {
        delete *rriter;
    }

	return EXIT_SUCCESS;
}

