/*
 <Mix-mesher: region type. This program generates a mixed-elements 2D mesh>

 Copyright (C) <2013,2020>  <Claudio Lobos> All rights reserved.

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
* @file Services.cpp
* @author Claudio Lobos, Fabrice Jaillet
* @version 0.1
* @brief
**/

#include "Services.h"
#include <string.h>
#include <iostream>
#include <fstream>
#include <algorithm>    // std::find_if
#include <chrono>

using Clobscode::Point3D;
using Clobscode::Polyline;
using Clobscode::MeshPoint;
using Clobscode::Quadrant;
using Clobscode::QuadEdge;
using std::vector;
using std::string;

namespace Clobscode
{

//-------------------------------------------------------------------
//-------------------------------------------------------------------

unsigned short Services::readRefinementRegions(string name,
                                               list<RefinementRegion *> &regions){

    char word [256];
    int cant;
    double x,y,z;
    unsigned short max_refinement = 0;

    vector<Point3D> cube_pts;

    FILE *file = fopen(name.c_str(),"r");

    if (file==NULL) {
        std::cout << "File " << name << " doesn't exist\n";
        return max_refinement;
    }

    //read number of nodes
    while(true){
        if(std::fscanf(file,"%s",word)==EOF){
            fclose(file);
            return max_refinement;
        }
        if(!strcmp(word,"n_regions\0"))
            break;
    }
    std::fscanf(file,"%i",&cant);

    //read each node
    cube_pts.reserve(cant*2);

    for( int i=0;i<cant;i++){

        //Read region's first point
        std::fscanf(file,"%s",word);
        x = atof(word);
        std::fscanf(file,"%s",word);
        y = atof(word);
        std::fscanf(file,"%s",word);
        z = atof(word);
        Point3D p1 (x,y,z);

        //Read region's second point
        std::fscanf(file,"%s",word);
        x = atof(word);
        std::fscanf(file,"%s",word);
        y = atof(word);
        std::fscanf(file,"%s",word);
        z = atof(word);
        Point3D p2 (x,y,z);

        //Read Refinement Level
        unsigned short rrl = 0;
        std::fscanf(file,"%s",word);
        rrl = atoi(word);

        if (rrl>max_refinement) {
            max_refinement = rrl;
        }

        RefinementRegion *rr = new RefinementCubeRegion(p1,p2,rrl);

        regions.push_back(rr);
    }

    fclose(file);

    return max_refinement;
}

//-------------------------------------------------------------------
//-------------------------------------------------------------------
bool Services::readSurfaceRefinementRegion(string name,
                                           list<RefinementRegion *> &regions,
                                           const unsigned short &rrl){

    vector<Polyline> tmp;
    if (!ReadPolyFile(name,tmp)) {
        return false;
    }
    RefinementRegion *rr = new RefinementSurfaceRegion(tmp[0],rrl);
    regions.push_back(rr);

    return true;
}

//-------------------------------------------------------------------
//-------------------------------------------------------------------
bool Services::readDrawingRefinementRegion(string name,
                                           list<RefinementRegion *> &regions,
                                           const unsigned short &rrl){

    vector<Polyline> tmp;
    if (!ReadPolyFile(name,tmp)) {
        return false;
    }
    RefinementRegion *rr = new RefinementDrawingRegion(tmp[0],rrl);
    regions.push_back(rr);

    return true;
}

//-------------------------------------------------------------------
//-------------------------------------------------------------------
bool Services::ReadOffMesh(string name,
                           vector<Clobscode::Polyline> &clobs_inputs){

    char word [256];
    int np, nf; //nb pts, nfaces=nb polylines
    double x,y,z;
    bool skel = false;
    vector<vector<unsigned int> > alledges;
    vector<Point3D> pts;

    FILE *file = fopen(name.c_str(),"r");

    if (file==NULL) {
        std::cout << "File " << name << " doesn't exist\n";
        fclose(file);
        return false;
    }

    while(true){
        if( fscanf(file,"%s",word) == EOF) {
            fclose(file);
            return false;
        }
        if (!strcmp(word,"")) {
            continue;
        }
        if(!strcmp(word,"OFF\0"))
            break;
        if(!strcmp(word,"2OFF\0"))
            break;
        if(!strcmp(word,"SKEL\0")){
            skel = true;
            break;
        }
        fclose(file);
        return false;
    }

    //read number of points
    fscanf(file,"%i",&np);
    //read number of faces
    fscanf(file,"%i",&nf);

    if (!skel) {
        //read number of edges [info not needed].
        fscanf(file,"%s",word);
    }

    //read each node
    pts.reserve(np);

    for( int i=0;i<np;i++){
        std::fscanf(file,"%s",word);
        x=atof(word);
        std::fscanf(file,"%s",word);
        y=atof(word);
        std::fscanf(file,"%s",word);
        z=atof(word);
        Point3D p (x,y,z);
        pts.push_back(p);
        fgets(word,256,file);
    }

    //number of face points
    for( int i=0;i<nf;i++){
        unsigned int nfp;
        std::fscanf(file,"%u",&nfp);
        alledges.resize(nfp, std::vector<unsigned int>(2));

        int fpts;
        for(unsigned int j=0;j<nfp;j++){
            std::fscanf(file,"%i",&fpts);
            alledges[(j-1)%nfp][1]=fpts;
            alledges[j][0]=fpts;
        }
        //read any other data in the line
        fgets(word,256,file);

        Polyline pm (pts, alledges);
        clobs_inputs.push_back(pm);
        alledges.resize(0);
    }
    fclose(file);


    return true;
}

//-------------------------------------------------------------------
//-------------------------------------------------------------------
bool Services::ReadMdlMesh(std::string name,
                           vector<Clobscode::Polyline> &clobs_inputs){

    char word [256];
    int cant;
    bool isEdges=false; //true if mdl edge format
    bool isPolyline=false; //true if mdl polyline format
    double x,y,z;
    vector<vector<unsigned int> > alledges;
    vector<Point3D> pts;

    FILE *file = fopen(name.c_str(),"r");

    if (file==NULL) {
        std::cout << "File " << name << " doesn't exist\n";
        return false;
    }

    //read number of nodes
    while(true){
        if(std::fscanf(file,"%s",word)==EOF){
            fclose(file);
            return false;
        }
        if(!strcmp(word,"ARRAY1<POINT3D>]\0"))
            break;
    }
    std::fscanf(file,"%i",&cant);

    if(cant<=0)
        return false;
    //read each node
    pts.reserve(cant);

    for( int i=0;i<cant;i++){
        std::fscanf(file,"%s",word);
        x=atof(word);
        std::fscanf(file,"%s",word);
        y=atof(word);
        std::fscanf(file,"%s",word);
        z=atof(word);
        Point3D p (x,y,z);
        pts.push_back(p);
    }

    //read number of "Edges"
    cant = 0;
    while(1){
        if(std::fscanf(file,"%s",word) == EOF){
            std::cout << "didn't find polylines\n";
            fclose(file);
            return false;
        }

        if(strstr(word,"Edges"))
            isEdges=true;
        else if(strstr(word,"Polyline"))
            isPolyline=true;
        else if(!strcmp(word,"ARRAY1<STRING>]\0")){
            //std::fscanf(file,"%s",word);
            std::fscanf(file,"%i",&cant);
            break;
        }
    }

    if (isEdges==true) { // list of edges, assume CCW
        alledges.reserve(cant);
        //read each edge (assuming they have 2 endpoints)
        int dust;
        for( int i=0;i<cant;i++){
            std::vector<unsigned int> edgpts(2,0);
            for(unsigned int j=0;j<2;j++){
                std::fscanf(file,"%u",&edgpts[j]);
            }
            //read some unnecessary for the moment integers
            for(unsigned int j=0;j<1;j++)
                std::fscanf(file,"%i",&dust);

            alledges.push_back(edgpts);
        }
    } else if (isPolyline==true) { // one polyline, nodes 0 1 2 3 4, assume CCW
        exit(5); // TODO
    } else { // no edges nor polyline, just circulate nodes, assuming CCW
        exit(6); // TODO
    }
    fclose(file);

    Polyline pm (pts, alledges);
    clobs_inputs.push_back(pm);

    return true;
}


//-------------------------------------------------------------------
//-------------------------------------------------------------------
//        First line: <# of vertices> <dimension (must be 2)> <# of attributes> <# of boundary markers (0 or 1)>
//        Following lines: <vertex #> <x> <y> [attributes] [boundary marker]
//        One line: <# of segments> <# of boundary markers (0 or 1)>
//        Following lines: <segment #> <endpoint> <endpoint> [boundary marker]
//        One line: <# of holes>
//        Following lines: <hole #> <x> <y>
//        Optional line: <# of regional attributes and/or area constraints>
//        Optional following lines: <region #> <x> <y> <attribute> <maximum area>

bool Services::ReadPolyFile(std::string name,
                            vector<Clobscode::Polyline> &clobs_inputs){

    char word [256];
    int cant;
    int shift=-1; //trick for determining shift from first point indice (0 or 1)
    double x,y,z;
    vector<vector<unsigned int> > alledges;
    vector<Point3D> pts;

    FILE *file = fopen(name.c_str(),"r");

    if (file==NULL) {
        std::cout << "File " << name << " doesn't exist\n";
        return false;
    }

    //read first line
    while(true) {
        if(std::fscanf(file,"%s",word)==EOF){
            fclose(file);
            return false;
        }
        if(!strncmp(word,"#",1)) {
            //skip end of line
            fscanf(file,"%2000[^\n]", word);
        }
        else
            break;
    }
    cant=atoi(word);
    std::fscanf(file,"%s",word);
    int dim=atoi(word);
    std::fscanf(file,"%s",word);
    int nattribs=atoi(word);
    std::fscanf(file,"%s",word);
    int nmarkers=atoi(word);

    if(cant<=0)
        return false;
    //read each node
    pts.reserve(cant);
    int i=0;
    do {
        if(std::fscanf(file,"%s",word)==EOF){
            fclose(file);
            return false;
        }
        if(!strncmp(word,"#",1)) {
            //skip end of line
            fscanf(file,"%2000[^\n]", word);
        }
        else {
            if (shift==-1) // if first point
                shift=atoi(word);
            // else, discard indice of all following nodes contained in word
            std::fscanf(file,"%s",word);
            x=atof(word);
            std::fscanf(file,"%s",word);
            y=atof(word);
            if (dim>=3) { // should not be, but one never knows...
                std::fscanf(file,"%s",word);
                z=atof(word);
            } else {
                z=0.;
            }
            for (int j=0; j<nattribs;++j) {
                std::fscanf(file,"%s",word);
                //std::cerr << "Skipping this node attribute " << word << "\n";
            }
            for (int j=0; j<nmarkers;++j) {
                std::fscanf(file,"%s",word);
                //std::cerr << "Skipping this node marker " << word << "\n";
            }

            pts.push_back( Point3D(x,y,z));

            ++i;
        }
    } while(i!=cant);

    //read number of "Edges"
    cant = 0;
    //read first line
    while(true) {
        if(std::fscanf(file,"%s",word)==EOF){
            fclose(file);
            return false;
        }
        if(!strncmp(word,"#",1)) {
            //skip end of line
            fscanf(file,"%2000[^\n]", word);
        }
        else
            break;
    }
    cant=atoi(word);
    std::fscanf(file,"%s",word);
    nmarkers=atoi(word);

    alledges.reserve(cant);
    if(cant<=0)
        return false;
    //read each edge
    i=0;
    do {
        if(std::fscanf(file,"%s",word)==EOF){
            fclose(file);
            return false;
        }
        if(!strncmp(word,"#",1)) {
            //skip end of line
            fscanf(file,"%2000[^\n]", word);
        }
        else {
            // here, discard indice of the edge contained in word

            //read each edge (assuming they have 2 endpoints)
            std::vector<unsigned int> edgpts(2,0);
            for(unsigned int j=0;j<2;j++){
                std::fscanf(file,"%u",&edgpts[j]); // 1-j if inverted edges
                edgpts[j]-=shift;                  // 1-j if inverted edges
            }
            //read some unnecessary for the moment integers
            for(int j=0;j<nmarkers;j++) {
                std::fscanf(file,"%s",word);
                //std::cerr << "Skipping this edge marker " << word << "\n";
            }
            alledges.push_back(edgpts);
            ++i;
        }
    } while (i!=cant);

    //read number of "Holes"
    cant = 0;
    //read first line
    while(true) {
        if(std::fscanf(file,"%s",word)==EOF){
            fclose(file);
            return true;
        }
        if(!strncmp(word,"#",1)) {
            //skip end of line
            fscanf(file,"%2000[^\n]", word);
        }
        else
            break;
    }
    cant=atoi(word);
    if (cant>0) {
        std::cerr << "ReadPolyFile(): hole feature not yet supported... "
                  << cant << " hole(s) skipped.\n";
        std::cerr << "  and possibly some regions too.\n";
    }

    fclose(file);

    Polyline pm (pts, alledges);
    clobs_inputs.push_back(pm);

    return true;
}


//-------------------------------------------------------------------
//-------------------------------------------------------------------
bool Services::ReadQuadrantList(std::string name, list<unsigned int> &olist,
                                vector<unsigned int> &ele_quadt_ref) {
    FILE *file = fopen(name.c_str(),"r");
    unsigned int idx = 0;
    char word [256];
    while (std::fscanf(file,"%s",word) != EOF) {
        bool num = true;
        for (unsigned int i=0; word[i]!= '\0'; i++) {
            if (!isdigit(word[i])) {
                num = false;
                break;
            }
        }
        if(num) {
            idx = atoi(word);
            if (ele_quadt_ref.size()<=idx) {
                cerr << "Invalid element index while reading list of elements to refine\n";
                cerr << "Quadtree mesh must be readed before the list is provided\n";
                cerr << "Aborting!!!\n";
                std::abort();
            }

            olist.push_back(ele_quadt_ref[idx]);
        }
    }
    olist.sort();
    olist.unique();
    fclose(file);
    return true;
}



//-------------------------------------------------------------------
//-------------------------------------------------------------------
bool Services::ReadQuadMesh(std::string name, vector<MeshPoint> &points,
                            vector<Quadrant> &Quadrants,
                            map<QuadEdge, EdgeInfo> &edges,
                            vector<unsigned int> &ele_oct_ref,
                            GeometricTransform &gt,
                            unsigned short &minrl,
                            unsigned short &maxrl) {
    auto start_time = chrono::high_resolution_clock::now();

    char word [256];
    double x,y,z;
    unsigned int e1,e2,e3,e4,e5;
    unsigned int elem;
    unsigned int np=0, ne=0, no=0, nl=0;

    FILE *file = fopen(name.c_str(),"r");

    if (file==NULL) {
        std::cout << "File " << name << " doesn't exist\n";
        return false;
    }

    //read header
    std::fscanf(file,"%u",&np);
    std::fscanf(file,"%u",&ne);
    std::fscanf(file,"%u",&no);

    //read each node
    points.reserve(np);

    for(unsigned int i=0;i<np;i++){
        std::fscanf(file,"%s",word);
        x=atof(word);
        std::fscanf(file,"%s",word);
        y=atof(word);
        std::fscanf(file,"%s",word);
        z=atof(word);
        points.push_back(Point3D (x,y,z));
    }

    //read edges
    for(unsigned int i=0;i<ne;i++){
        std::fscanf(file,"%u",&e1);
        std::fscanf(file,"%u",&e2);
        std::fscanf(file,"%u",&e3);
        //Neighbor Quads must be re-computed in method: Mesher::setInitialState
        //we leave them as 0 index for the moment.
        unsigned int non = std::numeric_limits<unsigned int>::max();
        edges.emplace(QuadEdge (e1,e2),EdgeInfo (e3,non,non));
    }

    //read the Quadrants, its refinement level and
    //the input faces intersected by it.
    Quadrants.reserve(no);
    unsigned int nop = 0, nof = 0, ni=0;
    int orl = 0;

    minrl = 100;
    maxrl = 0;

    //unsigned int noregular=0;
    for (unsigned int i=0; i<no; i++) {
        vector<unsigned int> opts;
        list<unsigned int> ofcs;
        std::fscanf(file,"%s",word);
        nop = atoi(word);
        opts.reserve(nop);

        if (nop!=4) {
            cerr << "warning at Services::ReadQuadMesh\n";
            cerr << "         Quadrant hasn't 4 nodes: " << nop << "\n";
            cout << "Quadrant index " << i << "\n";
            continue;
        }

        for (unsigned int j=0; j<nop; j++) {
            std::fscanf(file,"%u",&ni);
            opts.push_back(ni);
        }
        std::fscanf(file,"%i",&orl);
        std::fscanf(file,"%u",&nof);
        for (unsigned int j=0; j<nof; j++) {
            std::fscanf(file,"%u",&ni);
            ofcs.push_back(ni);
        }
        Quadrant Quadrant (opts,orl,i);
        Quadrant.setIntersectedEdges(ofcs);
        Quadrants.push_back(Quadrant);

        if (orl<minrl) {
            minrl = short(orl);
        }
        if (orl>maxrl) {
            maxrl = short(orl);
        }
    }

    std::fscanf(file,"%s %s",word,word);
    std::fscanf(file,"%s",word);
    x=atof(word);
    std::fscanf(file,"%s",word);
    y=atof(word);
    std::fscanf(file,"%s",word);
    z=atof(word);
    Point3D c(x,y,z);
    std::fscanf(file,"%s",word);
    x=atof(word);
    std::fscanf(file,"%s",word);
    y=atof(word);
    std::fscanf(file,"%s",word);
    z=atof(word);
    gt.setCentroid(c);
    gt.setXAxis(x);
    gt.setYAxis(y);
    gt.setZAxis(z);
    
    //read the element Quadrant link
    ele_oct_ref.reserve(no);
    //            unsigned int checksum = 0;
    for (unsigned int i=0; i<no; i++) {
        std::fscanf(file,"%u",&elem);
        for (unsigned int j=0; j<elem; j++) {
            ele_oct_ref.push_back(i);
        }
    }

    fclose(file);

    auto end_time = chrono::high_resolution_clock::now();
    cout << "  Read  done in "
         << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count() ;
    cout << " ms"<< endl;

    return true;
}

//-------------------------------------------------------------------
//-------------------------------------------------------------------
bool Services::WriteQuadtreeMesh(std::string name, const vector<MeshPoint> &points,
                                 const vector<Quadrant> &Quadrants,
                                 const map<QuadEdge, EdgeInfo> &edges,
                                 const GeometricTransform &gt){

    auto start_time = chrono::high_resolution_clock::now();

    //            QuadEdge qe;
    //<QuadEdge>::const_iterator my_edge;

    string vol_name = name+".oct";

    //write the surface mesh
    FILE *f = fopen(vol_name.c_str(),"wt");
    unsigned int np = points.size();
    unsigned int no = Quadrants.size();
    unsigned int ne = edges.size();

    fprintf(f,"%u %u %u\n\n", np, ne, no);

    //write points
    for(unsigned int i=0;i<np;i++){
        Point3D p = points[i].getPoint();
        fprintf(f,"%+1.8E %+1.8E %+1.8E\n",p[0],p[1],p[2]);
    }
    fprintf(f,"\n");

    //write edges
    for (const auto qe: edges) {
        //index of both nodes defining the edge, as well as its mid-node
        fprintf(f,"%u %u %u\n",qe.first[0],qe.first[1],(qe.second)[0]);
        //the two quad neighbor index must be recomputed at reading procees,
        //therefore we wont write them in the file.
    }
    
    fprintf(f,"\n");

    //pair sub-elements with Quadrant index.
    //this info is printed per Quadrant and the elements are
    //printed in order in the mesh file so we can compute
    //for each element to which Quadrant it belongs.
    /*for (unsigned int i=0; i<Quadrants.size(); i++) {
        unsigned int nse = Quadrants[i].getSubElements().size();
        fprintf(f,"%u ",nse);
    }
    fprintf(f,"\n\n");*/

    for (unsigned int i=0; i<Quadrants.size(); i++) {
        const vector<unsigned int> &opts = Quadrants[i].getPointIndex();
        unsigned int nopts = opts.size();
        if (nopts<4) {
            cerr << "warning at Services::WriteQuadMesh\n";
            cerr << "        Quadrant has less than 4 nodes\n";
            fprintf(f,"%u ",nopts);
        }
        else {
            nopts=4;
            fprintf(f,"4 ");
        }
        
        for (unsigned int j=0; j<nopts; j++) {
            fprintf(f,"%u ",opts[j]);
        }
        
        fprintf(f,"%u\n",Quadrants[i].getRefinementLevel());
        list<unsigned int> ofcs = Quadrants[i].getIntersectedEdges();
        list<unsigned int>::iterator fiter;
        nopts = ofcs.size();
        fprintf(f,"%u",nopts);
        for (fiter=ofcs.begin(); fiter!=ofcs.end(); fiter++) {
            fprintf(f," %u",*fiter);
        }
        fprintf(f,"\n");
    }

    fprintf(f,"\nGeometric Transform\n");
    Point3D c = gt.getCentroid();
    fprintf(f,"%f %f %f\n",c[0],c[1],c[2]);
    fprintf(f,"%f %f %f\n\n",gt.getXAxis(),gt.getYAxis(),gt.getZAxis());

    fclose(f);

    auto end_time = chrono::high_resolution_clock::now();
    cout << "    * WriteQuadtreeMesh in "
         << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count() ;
    cout << " ms"<< endl;

    return true;
}
    
    
    //-------------------------------------------------------------------
    //-------------------------------------------------------------------
    bool Services::addOctElemntInfo(std::string name, vector<Quadrant> &quadrants,
                                    map<unsigned, bool> &removedquads,
                                    const list<unsigned int> &quadmeshidx){
        
        string vol_name = name+".oct";
        
        //write the volume mesh
        FILE *f = fopen(vol_name.c_str(),"a");
        unsigned int no = quadrants.size();
        
        
        
        
        //Now we must update the .oct file with information linking
        //elements to octants. To this purpose we update the removedoct map
        //setting to false any octant remaining in the octant vector.
        //We also need a map from octant index -> octant position in vector
        map<unsigned int, unsigned int> quadpos;
        for (unsigned int i=0;i<quadrants.size();i++) {
            unsigned int oi = quadrants[i].getIndex();
            removedquads[oi] = false;
            quadpos[oi] = i;
        }
        
        
        
        //pair sub-elements with octant index.
        //this info is printed per octant and the elements are
        //printed in order in the mesh file so we can compute
        //for each element to which octant it belongs.
        //If 0 is printed means that the octant has 0 sub-elements.
        //This octant was later removed in the mesh generation process
        //due to proximity to the boundary. However, to refine a mesh
        //from an existing one, it is still necessary.
        for (auto o: quadmeshidx) {
            if (removedquads[o]) {
                fprintf(f,"%u ",0);
            }
            else {
                unsigned int nse = quadrants[quadpos[o]].getSubElements().size();
                fprintf(f,"%u ",nse);
            }
        }
        
        fprintf(f,"\n");
        
        fclose(f);
        
        return true;
    }


//-------------------------------------------------------------------
// http://gmsh.info/doc/texinfo/gmsh.html#MSH-ASCII-file-format
//-------------------------------------------------------------------
bool Services::WriteGMSH(std::string name, const shared_ptr<FEMesh> &output){

    vector<Point3D> points = output->getPoints();
    vector<vector<unsigned int> > elements = output->getElements();

    if (elements.empty()) {
        std::cerr << "no output elements\n";
        return false;
    }

    string vol_name = name+".gmsh";

    //write the volume mesh
    FILE *f = fopen(vol_name.c_str(),"wt");

    fprintf(f,"$MeshFormat\n2.2 0 8\n$EndMeshFormat\n");

    fprintf(f,"$Nodes\n%i",(int)points.size());

    //write points
    for(unsigned int i=0;i<points.size();i++){
        fprintf(f,"\n%i",i+1);
        fprintf(f," %+1.8E",points[i][0]);
        fprintf(f," %+1.8E",points[i][1]);
        fprintf(f," %+1.8E",points[i][2]);
    }
    fprintf(f,"\n$EndNodes\n");

    fprintf(f,"$Elements\n%i",(int)elements.size());

    //get all the elements in a std::vector
    for (unsigned int i=0; i<elements.size(); i++) {
        fprintf(f,"\n%i",i+1);
        std::vector<unsigned int> epts = elements[i];
        unsigned int np = epts.size();
        if (np == 3) {
            fprintf(f," 2 2 99 2"); //TRIANGLE
        }
        else if (np == 4){ //QUAD
            fprintf(f," 3 2 99 2");
        }

        for (unsigned int j= 0; j<np; j++) {
            fprintf(f," %i", epts[j]+1);
        }
    }
    fprintf(f,"\n$EndElements\n");

    fclose(f);

    return true;
}

//-------------------------------------------------------------------
//-------------------------------------------------------------------
bool Services::WriteVTK(std::string name, const shared_ptr<FEMesh> &output){

    vector<Point3D> points = output->getPoints();
    vector<vector<unsigned int> > elements = output->getElements();

    if (elements.empty()) {
        std::cerr << "no output elements\n";
        return false;
    }

    string vol_name = name+".vtk";

    //write the volume mesh
    FILE *f = fopen(vol_name.c_str(),"wt");

    fprintf(f,"# vtk DataFile Version 2.0\nvtk file: 2D Unstructured Grid %s\nASCII",name.c_str());
    fprintf(f,"\n\nDATASET UNSTRUCTURED_GRID\nPOINTS %i float",(int)points.size());

    //write points
    for(unsigned int i=0;i<points.size();i++){
        if (i%2==0) {
            fprintf(f,"\n");
        }
        fprintf(f," %+1.8E",points[i][0]);
        fprintf(f," %+1.8E",points[i][1]);
        fprintf(f," %+1.8E",points[i][2]);
    }

    //count conectivity index.
    unsigned int conectivity = 0;
    for (unsigned int i=0; i<elements.size(); i++) {
        conectivity+=elements[i].size()+1;
    }

    fprintf(f,"\n\nCELLS %i %i\n",(int)elements.size(),conectivity);

    //get all the elements in a std::vector
    for (unsigned int i=0; i<elements.size(); i++) {
        std::vector<unsigned int> epts = elements[i];
        unsigned int np = epts.size();
        fprintf(f,"%i", np);

        for (unsigned int j= 0; j<np; j++) {
            fprintf(f," %i", epts[j]);
        }

        fprintf(f,"\n");
    }

    fprintf(f,"\nCELL_TYPES %i",(int)elements.size());
    for (unsigned int i=0; i<elements.size(); i++) {
        if (i%30==0) {fprintf(f,"\n"); }
        unsigned int np = elements[i].size();
        if (np == 3) {
            fprintf(f," 5"); //VTK_TRIANGLE
        }
        else if (np == 4){ //VTK_QUAD
            fprintf(f," 9");
        }
    }

    fprintf(f,"\nCELL_DATA %i\n",(int)elements.size());
    fprintf(f,"SCALARS elemType int 1\n");
    fprintf(f,"LOOKUP_TABLE my_color");
    for (unsigned int i=0; i<elements.size(); i++) {
        if (i%30==0) {fprintf(f,"\n");}
        fprintf(f," %lu",elements[i].size()-3);
    }
    fprintf(f,"\nLOOKUP_TABLE my_color 2\n");
    fprintf(f,"1.0 0.0 1.0 1.0\n");
    fprintf(f,"0.0 1.0 1.0 1.0\n");

    
    //write surface state if computed before (-q option, decoration==true)
    const vector <unsigned short> &surfele= output->getSurfState();
    if (surfele.size()>0) {
        fprintf(f,"\nSCALARS surfState int 1");
        fprintf(f,"\nLOOKUP_TABLE default");
        for (unsigned int i=0; i<surfele.size(); i++) {
            if (i%30==0) {fprintf(f,"\n");}
            fprintf(f," %u",surfele[i]);
        }
    }
    
    //write refinement levels if computed before (-q option, decoration==true)
    const vector <unsigned short> &reflevels= output->getRefLevels();
    if (reflevels.size()>0) {
        fprintf(f,"\nSCALARS refLevel int 1");
        fprintf(f,"\nLOOKUP_TABLE default");
        for (unsigned int i=0; i<reflevels.size(); i++) {
            if (i%30==0) {fprintf(f,"\n");}
            fprintf(f," %u",reflevels[i]);
        }
    }
    //write minAngles if computed before (-q option, decoration==true)
    const vector <double> &minAngles= output->getMinAngles();
    double minAngleTri=std::numeric_limits<int>::max();
    double minAngleQuad=std::numeric_limits<int>::max();
    double maxAngleTri=std::numeric_limits<int>::min();
    double maxAngleQuad=std::numeric_limits<int>::min();
    if (minAngles.size()>0) {
        fprintf(f,"\n\nSCALARS minAngle int 1");
        fprintf(f,"\nLOOKUP_TABLE min");
        for (unsigned int i=0; i<minAngles.size(); i++) {
            int angle;
            if (i%30==0) {fprintf(f,"\n");}
            if (elements[i].size()==3) { //triangle
                fprintf(f," %d", (int) (60.-minAngles[i]));
                minAngleTri=min(minAngleTri,minAngles[i]);
            }
            else { //quad
                fprintf(f," %d", (int) (90.-minAngles[i]));
                minAngleQuad=min(minAngleQuad,minAngles[i]);
            }
        }
        cout << "    min Tri and Quad angles "
        << minAngleTri << " " << minAngleQuad << endl;
    }

    //write maxAngles if computed before (-q option, decoration==true)
    const vector <double> &maxAngles= output->getMaxAngles();
    if (maxAngles.size()>0) {
        fprintf(f,"\n\nSCALARS maxAngle int 1");
        fprintf(f,"\nLOOKUP_TABLE max");
        for (unsigned int i=0; i<maxAngles.size(); i++) {
            int angle;
            if (i%30==0) {fprintf(f,"\n");}
            if (elements[i].size()==3) { //triangle
                fprintf(f," %d", (int) (maxAngles[i]-60.));
                maxAngleTri=max(maxAngleTri,maxAngles[i]);
            }
            else { //quad
                fprintf(f," %d", (int) (maxAngles[i]-90.));
                maxAngleQuad=max(maxAngleQuad,maxAngles[i]);
            }
        }
        cout << "    max Tri and Quad angles "
        << maxAngleTri << " " << maxAngleQuad << endl;
    }

    //write debugging state if computed before (-q option, decoration==true)
    const vector <unsigned short> &debugging= output->getDebugging();
    if (debugging.size()>0) {
        fprintf(f,"\n\nSCALARS debugging int 1");
        fprintf(f,"\nLOOKUP_TABLE default");
        for (unsigned int i=0; i<surfele.size(); i++) {
            if (i%30==0) {fprintf(f,"\n");}
            fprintf(f," %u",debugging[i]);
        }
    }

    fclose(f);

    return true;
}

//-------------------------------------------------------------------
//-------------------------------------------------------------------
bool Services::WriteOFF(std::string name, const shared_ptr<FEMesh> &output){

    const vector<Point3D> &points = output->getPoints();
    const vector<vector<unsigned int> > &elements = output->getElements();
    const vector<unsigned int> &colored = output->getColoredCells();

    if (elements.empty()) {
        std::cerr << "no output elements\n";
        return false;
    }

    string vol_name = name+".off";

    //write the volume mesh
    FILE *f = fopen(vol_name.c_str(),"wt");

    fprintf(f,"OFF %i %i 0\n\n",(int)points.size(),(int)elements.size());

    //write points
    for(unsigned int i=0;i<points.size();i++){
        fprintf(f," %f",points[i][0]);
        fprintf(f," %f",points[i][1]);
        fprintf(f," %f\n",points[i][2]);
    }

    fprintf(f,"\n");

    if (colored.empty()) {
        for (unsigned int i=0; i<elements.size(); i++) {
            const std::vector<unsigned int> &epts = elements[i];
            unsigned int np = epts.size();
            fprintf(f,"%i", np);

            for (unsigned int j= 0; j<np; j++) {
                fprintf(f," %i", epts[j]);
            }

            if (np==3) {
                fprintf(f," 0 1 0");
            }

            fprintf(f,"\n");
        }
        return true;
    }

    //get all the elements in a std::vector
    for (unsigned int i=0; i<elements.size(); i++) {
        const std::vector<unsigned int> &epts = elements[i];
        unsigned int np = epts.size();
        fprintf(f,"%i", np);

        for (unsigned int j= 0; j<np; j++) {
            fprintf(f," %i", epts[j]);
        }

        if (colored[i]!=0) {
            fprintf(f," 1 0 0");
        }
        else {
            if (np==3) {
                fprintf(f," 0 1 0");
            }
        }

        fprintf(f,"\n");
    }

    fclose(f);

    return true;
}
//-------------------------------------------------------------------
//-------------------------------------------------------------------
bool Services::WritePolyVTK(std::string name, vector<Polyline> inputs) {

    string vol_name = name+".poly.vtk";

    //write the volume mesh
    FILE *f = fopen(vol_name.c_str(),"wt");

    // count the total number of Points, if many Polylines in Input
    vector <unsigned int> noPts(inputs.size());
    unsigned int noPtsTot=0, shiftNoPts=0;
    for (unsigned int i=0;i<inputs.size();i++) {
        noPts[i]=inputs[i].getPoints().size();
        noPtsTot+=noPts[i];
    }
    if (noPtsTot==0) {
        std::cout << "How strange, no inputs points to write...\n";
        return false;
    }

    fprintf(f,"# vtk DataFile Version 2.0\nPolydata %s\nASCII",name.c_str());
    fprintf(f,"\n\nDATASET POLYDATA\nPOINTS %u float", noPtsTot);

    //for each input polyline
    for (const Polyline& ply:inputs) {

        //write points
        const vector<Point3D> &points=ply.getPoints();
        for(unsigned int i=0;i<points.size();i++){
            if (i%2==0) {
                fprintf(f,"\n");
            }
            fprintf(f," %+1.8E",points[i][0]);
            fprintf(f," %+1.8E",points[i][1]);
            fprintf(f," %+1.8E",points[i][2]);
        }
    }

    //tricky part to recover connected subpolylines (polygons) in a polyline
    vector <vector< unsigned int> > polygons;
    uint iPly=0, noEdgesTot=0;
    for (const Polyline &ply:inputs) {
        // get connected polylines from input, assumming they are in order
        // make a copy of edges
        list<PolyEdge> edges(ply.getEdges().begin(),ply.getEdges().end());
        noEdgesTot+=edges.size(); //data used at final to write VTK info
        while (edges.size()>0) {
            polygons.push_back(vector<unsigned int>());
            auto e=edges.front(); edges.pop_front(); // get first edge and remove it from list
            uint start=e.getKey();
            do {
                // keep the first point
                polygons[iPly].push_back(shiftNoPts+e.getKey());
                // and find to next edge that has the second pointas first point
                uint val=e.getVal();
                auto it=std::find_if(edges.begin(), edges.end(),
                                     [&val] (const PolyEdge& pe) -> bool {return pe.getKey() == val; });
                // go to that edge and remove it from list
                e=*it;
                edges.erase(it);
            } while (start!=e.getVal());
            polygons[iPly].push_back(shiftNoPts+e.getKey()); // insert last point
            ++iPly; //that will start a new polygon
        }
        shiftNoPts+=ply.getPoints().size(); //don't forget to shift local numbering to global
    }


    fprintf(f,"\n\nPOLYGONS %lu %lu\n",polygons.size(),polygons.size()+noEdgesTot);
    for (auto ply:polygons) {

        //write index of the polygon points, as processed right above
        fprintf(f,"%lu", ply.size());
        for (unsigned int i=0; i<ply.size(); i++) {
            fprintf(f," %u", ply[i]);
        }
        fprintf(f,"\n");
    }

    fclose(f);

    return true;
}
//-------------------------------------------------------------------
//-------------------------------------------------------------------
bool Services::WritePolyFile(std::string name, vector<Polyline> inputs) {

    string vol_name = name+".poly.poly";

    //write the volume mesh
    FILE *f = fopen(vol_name.c_str(),"wt");

    // count the total number of Points, if many Polylines in Input
    vector <unsigned int> noPts(inputs.size());
    unsigned int noPtsTot=0, shiftNoPts=0;
    for (unsigned int i=0;i<inputs.size();i++) {
        noPts[i]=inputs[i].getPoints().size();
        noPtsTot+=noPts[i];
    }
    if (noPtsTot==0) {
        std::cout << "How strange, no inputs points to write...\n";
        return false;
    }

    fprintf(f,"# Triangle\n# %s\n",name.c_str());
    fprintf(f,"%u 2 0 0\n", noPtsTot);

    //for each input polyline
    for (const Polyline& ply:inputs) {

        //write points
        const vector<Point3D> &points=ply.getPoints();
        for(unsigned int i=0;i<points.size();i++){
            fprintf(f," %u",i);
            fprintf(f," %+1.8E",points[i][0]);
            fprintf(f," %+1.8E\n",points[i][1]);
            //fprintf(f," %+1.8E",points[i][2]);
        }
    }

    uint noEdgesTot=0;
    for (const Polyline &ply:inputs) {
        noEdgesTot+=ply.getEdges().size(); //data used at final to write VTK info
    }

    fprintf(f,"%u 0\n",noEdgesTot);

    //tricky part to recover connected subpolylines (polygons) in a polyline
    vector <vector< unsigned int> > polygons;
    for (const Polyline &ply:inputs) {
        const vector<PolyEdge> &qedges = ply.getEdges();
        for (unsigned int i=0; i<qedges.size(); ++i) {
            fprintf(f," %u",i);
            fprintf(f," %u",shiftNoPts+qedges[i].getKey());
            fprintf(f," %u\n",shiftNoPts+qedges[i].getVal());
        }
        shiftNoPts+=ply.getPoints().size(); //don't forget to shift local numbering to global
    }

    fprintf(f,"0\n"); //no holes

    fclose(f);

    return true;
}

//-------------------------------------------------------------------
//-------------------------------------------------------------------
bool Services::WriteMixedVolumeMesh(std::string name, const shared_ptr<FEMesh> &output){

    const vector<Point3D> &points = output->getPoints();
    const vector<vector<unsigned int> > &elements = output->getElements();

    if (elements.empty()) {
        std::cerr << "no output elements\n";
        return false;
    }

    string vol_name = name+".mvm";

    //write the volume mesh
    FILE *f = fopen(vol_name.c_str(),"wt");

    //            unsigned int n = points.size();

    fprintf(f,"MIXED\n%i %i\n\n",(int)points.size(),(int)elements.size());

    //write points
    for(unsigned int i=0;i<points.size();i++){
        fprintf(f,"%+1.8E",points[i][0]);
        fprintf(f," %+1.8E",points[i][1]);
        fprintf(f," %+1.8E\n",points[i][2]);
    }

    //get all the elements in a std::vector
    for (unsigned int i=0; i<elements.size(); i++) {
        std::vector<unsigned int> epts = elements[i];
        unsigned int np = epts.size();
        if (np == 3) {
            fprintf(f,"T");
        }
        else if (np == 4){
            fprintf(f,"Q");
        }

        for (unsigned int j= 0; j<np; j++) {
            fprintf(f," %i", epts[j]);
        }

        fprintf(f,"\n");
    }
    fclose(f);

    return true;
}

//-------------------------------------------------------------------
//-------------------------------------------------------------------
bool Services::WriteOutputMesh(std::string name, const shared_ptr<FEMesh> &output){

    vector<Point3D> points = output->getPoints();
    vector<vector<unsigned int> > elements = output->getElements();

    if (elements.empty()) {
        std::cerr << "no output elements\n";
        return false;
    }

    string vol_name = name+".m2d";

    //write the volume mesh
    FILE *f = fopen(vol_name.c_str(),"wt");

    unsigned int n = points.size();

    fprintf(f,"%s\n","[Nodes, ARRAY1<STRING>]");
    fprintf(f,"%i\n\n",n);

    //write points
    for(unsigned int i=0;i<n;i++){
        fprintf(f,"1 %+1.8E",points[i][0]);
        fprintf(f," %+1.8E",points[i][1]);
        fprintf(f," %+1.8E\n",points[i][2]);
    }

    n = elements.size();

    fprintf(f,"\n%s\n","[Elements, ARRAY1<STRING>]");
    fprintf(f,"%i\n\n",n);

    //get all the elements in a std::vector
    for (unsigned int i=0; i<n; i++) {
        std::vector<unsigned int> epts = elements[i];
        unsigned int np = epts.size();
        // lowercase to avoid confusion with m3d Tetra=uppercase T
        if (np == 3) {
            fprintf(f,"t"); // triangle
        }
        else if (np == 4){
            fprintf(f,"q"); // quadrangle
        }

        for (unsigned int j= 0; j<np; j++) {
            fprintf(f," %i", epts[j]);
        }

        fprintf(f," 1000.0 0.45 1.0\n");
    }
    fclose(f);

    return true;
}


//-------------------------------------------------------------------
//-------------------------------------------------------------------
bool Services::WriteMeshGetfem(std::string name, const std::shared_ptr<FEMesh> &output){

    /* 2D file */

    vector<Point3D> points = output->getPoints();
    vector<vector<unsigned int> > elements = output->getElements();

    if (elements.empty()) {
        std::cerr << "no output elements\n";
        return false;
    }

    string vol_name = name+".gmf";

    //write the volume mesh
    FILE *f = fopen(vol_name.c_str(),"wt");

    unsigned int n = points.size();

    fprintf(f,"%s\n%s\n\n\n\n","% GETFEM MESH FILE","% GETFEM VERSION 5.1");
    fprintf(f,"%s\n\n","BEGIN POINTS LIST");
    fprintf(f,"%s %i\n\n","POINT COUNT",n);

    //write points
    for(unsigned int i=0;i<n;i++){
        fprintf(f,"  POINT  %i %+1.8E",i,points[i][0]);
        fprintf(f," %+1.8E\n",points[i][1]);
        //fprintf(f," %+1.8E\n",points[i][2]);
    }

    n = elements.size();

    fprintf(f,"\n%s\n\n\n\n","END POINTS LIST");
    fprintf(f,"%s\n\n","BEGIN MESH STRUCTURE DESCRIPTION");

    //get all the elements in a std::vector
    for (unsigned int i=0; i<n; i++) {
        std::vector<unsigned int> epts = elements[i];
        unsigned int np = epts.size();

        //if (np==8 || np==4) {
        fprintf(f,"CONVEX %i    ",i);
        //}


        if (np == 3) {
            fprintf(f,"'GT_PK(2,1)'      ");
        }
        else {
            fprintf(f,"'GT_QK(2,1)'      ");
            unsigned int aux = epts[2];
            epts[2] = epts[3];
            epts[3] = aux;
        }

        for (unsigned int j= 0; j<np; j++) {
            fprintf(f," %i", epts[j]);
        }

        fprintf(f,"\n");
    }

    fprintf(f,"\n%s\n","END MESH STRUCTURE DESCRIPTION");
    fclose(f);

    return true;
}

//-------------------------------------------------------------------
//-------------------------------------------------------------------
// printing "Triangle like" angle Histogram
// see https://www.cs.cmu.edu/~quake/triangle.V.html
bool Services::WriteHistogram(std::string name, const shared_ptr<FEMesh> &output){

    vector<Point3D> points = output->getPoints();
    vector<vector<unsigned int> > elements = output->getElements();

    if (elements.empty()) {
        std::cerr << "no output elements\n";
        return false;
    }

    string hist_name = name+"_histo.txt";
    string hist_name2 = name+"_hplot.txt";

    //write the histo data
    FILE *f = fopen(hist_name.c_str(),"wt");
    FILE *f2 = fopen(hist_name2.c_str(),"wt");

    // printing "Triangle like" angle Histogram
    const vector <double> &minAngles= output->getMinAngles();
    const vector <double> &maxAngles= output->getMaxAngles();
    const array <unsigned int,180> &anglesTriHistogram= output->getAnglesTriHistogram();
    const array <unsigned int,180> &anglesQuadHistogram= output->getAnglesQuadHistogram();
    double minAngleTri=std::numeric_limits<int>::max();
    double minAngleQuad=std::numeric_limits<int>::max();
    double maxAngleTri=std::numeric_limits<int>::min();
    double maxAngleQuad=std::numeric_limits<int>::min();

    // computing min and max angles
    for (unsigned int i=0; i<minAngles.size(); i++) {
        int angle;
        if (elements[i].size()==3) { //triangle
            minAngleTri=min(minAngleTri,minAngles[i]);
            maxAngleTri=max(maxAngleTri,maxAngles[i]);
        }
        else { //quad
            minAngleQuad=min(minAngleQuad,minAngles[i]);
            maxAngleQuad=max(maxAngleQuad,maxAngles[i]);
        }
    }

    // printing "Triangle like" angle Histogram
    fprintf(f,"For Triangles:\n");
    fprintf(f,"Smallest angle:          %g", minAngleTri);
    fprintf(f,"   |  Largest angle:           %g", maxAngleTri);
    fprintf(f,"\nAngle histogram:\n");

    // case interval 10
//    for (unsigned int i=0; i<9; ++i) {
//        fprintf(f,"%3d - %3d degrees:\t %6d | ", i*10, (i+1)*10,anglesTriHistogram[i]);
//        fprintf(f,"%3d - %3d degrees:\t %6d\n", (i+9)*10, (i+10)*10, anglesTriHistogram[i+9]) ;
//    }
    // case interval 1
    for (unsigned int i=0; i<90; ++i) {
        fprintf(f,"%3d - %3d degrees:\t %6d | ", i, (i+1),anglesTriHistogram[i]);
        fprintf(f,"%3d - %3d degrees:\t %6d\n", (i+90), (i+91), anglesTriHistogram[i+90]) ;
    }
    fprintf(f,"For Quads:\n");
    fprintf(f,"Smallest angle:          %g", minAngleQuad);
    fprintf(f,"   |  Largest angle:           %g", maxAngleQuad);
    fprintf(f,"\nAngle histogram:\n");
    // case interval 10
//    for (unsigned int i=0; i<9; ++i) {
//        fprintf(f,"%3d - %3d degrees:\t %6d | ", i*10, (i+1)*10,anglesQuadHistogram[i]);
//        fprintf(f,"%3d - %3d degrees:\t %6d\n", (i+9)*10, (i+10)*10, anglesQuadHistogram[i+9]) ;
//    }
    // case interval 1
    for (unsigned int i=0; i<90; ++i) {
        fprintf(f,"%3d - %3d degrees:\t %6d | ", i, (i+1),anglesQuadHistogram[i]);
        fprintf(f,"%3d - %3d degrees:\t %6d\n", (i+90), (i+91), anglesQuadHistogram[i+90]) ;
    }

    // printing data for Gnuplot
    fprintf(f2,"Angle  Triangles  Quadrangles\n");
    for (unsigned int i=0; i<180; ++i) {
//        fprintf(f2," %d-%d",i*10,(i+1)*10);
// case interval 10        fprintf(f2," %d %3d %3d\n",i*10,anglesTriHistogram[i],anglesQuadHistogram[i]);
        fprintf(f2," %d %3d %3d\n",i,anglesTriHistogram[i],anglesQuadHistogram[i]);
    }
    fprintf(f2," %d %3d %3d\n",180,0,0);

    // closing
    fclose(f);
    fclose(f2);

    return true;
}
}
