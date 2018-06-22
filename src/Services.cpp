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
    tmp.reserve(1);
    if (!ReadMdlMesh(name,tmp)) {
        return false;
    }
    RefinementRegion *rr = new RefinementSurfaceRegion(tmp[0],rrl);
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
                            set<QuadEdge> &edges,
                            vector<unsigned int> &ele_oct_ref,
                            GeometricTransform &gt,
                            unsigned short &minrl,
                            unsigned short &maxrl) {

    char word [256];
    double x,y,z;
    int e1,e2,e3;
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
    std::fscanf(file,"%u",&nl);

    //read each node
    points.reserve(np);

    for(unsigned int i=0;i<np;i++){
        std::fscanf(file,"%s",word);
        x=atof(word);
        std::fscanf(file,"%s",word);
        y=atof(word);
        std::fscanf(file,"%s",word);
        z=atof(word);
        Point3D p (x,y,z);
        points.push_back(p);
    }

    //read edges
    for(unsigned int i=0;i<ne;i++){
        std::fscanf(file,"%i",&e1);
        std::fscanf(file,"%i",&e2);
        std::fscanf(file,"%i",&e3);
        QuadEdge oe(e1,e2);
        unsigned int mid = (unsigned int)e3;
        oe.setMidPoint(mid);
        edges.insert(oe);
    }

    //read the element Quadrant link
    ele_oct_ref.reserve(nl);
    //            unsigned int checksum = 0;
    for (unsigned int i=0; i<no; i++) {
        std::fscanf(file,"%u",&elem);
        for (unsigned int j=0; j<elem; j++) {
            ele_oct_ref.push_back(i);
        }
    }

    //read the Quadrants, its refinement level and
    //the input faces intersected by it.
    Quadrants.reserve(no);
    unsigned int nop = 0, nof = 0, ni=0;
    int orl = 0;

    minrl = 100;
    maxrl = 0;

    //            unsigned int noregular=0;
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
        Quadrant Quadrant (opts,orl);
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

    fclose(file);
    return true;
}

//-------------------------------------------------------------------
//-------------------------------------------------------------------
bool Services::WriteQuadtreeMesh(std::string name, const vector<MeshPoint> &points,
                                 vector<Quadrant> &Quadrants,
                                 set<QuadEdge> &edges,
                                 unsigned int nels,
                                 GeometricTransform &gt){

    //            QuadEdge qe;
    set<QuadEdge>::const_iterator my_edge;

    string vol_name = name+".oct";

    //write the surface mesh
    FILE *f = fopen(vol_name.c_str(),"wt");
    unsigned int np = points.size();
    unsigned int no = Quadrants.size();
    unsigned int ne = edges.size();

    fprintf(f,"%u %u %u %u\n\n", np, ne, no, nels);

    //write points
    for(unsigned int i=0;i<np;i++){
        Point3D p = points[i].getPoint();
        fprintf(f,"%+1.8E  %+1.8E  %+1.8E\n",p[0],p[1],p[2]);
    }
    fprintf(f,"\n");

    //write edges
    for(my_edge=edges.begin();my_edge!=edges.end();my_edge++){
        QuadEdge me = *my_edge;
        fprintf(f,"%i %i %i\n",me[0],me[1],me[2]);
    }
    fprintf(f,"\n");

    //pair sub-elements with Quadrant index.
    //this info is printed per Quadrant and the elements are
    //printed in order in the mesh file so we can compute
    //for each element to which Quadrant it belongs.
    for (unsigned int i=0; i<Quadrants.size(); i++) {
        unsigned int nse = Quadrants[i].getSubElements().size();
        fprintf(f,"%u ",nse);
    }
    fprintf(f,"\n\n");

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
    fprintf(f,"%f %f %f\n",gt.getXAxis(),gt.getYAxis(),gt.getZAxis());

    fclose(f);

    return true;
}

//-------------------------------------------------------------------
//-------------------------------------------------------------------
bool Services::WriteVTK(std::string name, FEMesh &output){

    vector<Point3D> points = output.getPoints();
    vector<vector<unsigned int> > elements = output.getElements();

    if (elements.empty()) {
        std::cout << "no output elements\n";
        return false;
    }

    string vol_name = name+".vtk";

    //write the volume mesh
    FILE *f = fopen(vol_name.c_str(),"wt");

    fprintf(f,"# vtk DataFile Version 2.0\nUnstructured Grid %s\nASCII",name.c_str());
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

    //write refinement levels if computed before (-q option, decoration==true)
    const vector <unsigned short> &reflevels= output.getRefLevels();
    if (reflevels.size()>0) {
        fprintf(f,"\nSCALARS refLevel int 1");
        fprintf(f,"\nLOOKUP_TABLE default");
        for (unsigned int i=0; i<reflevels.size(); i++) {
            if (i%30==0) {fprintf(f,"\n");}
            fprintf(f," %u",reflevels[i]);
        }
    }
    //write minAngles if computed before (-q option, decoration==true)
    const vector <double> &minAngles= output.getMinAngles();
    if (minAngles.size()>0) {
        fprintf(f,"\n\nSCALARS minAngle int 1");
        fprintf(f,"\nLOOKUP_TABLE min");
        for (unsigned int i=0; i<minAngles.size(); i++) {
            if (i%30==0) {fprintf(f,"\n");}
            if (elements[i].size()==3) //triangle
                fprintf(f," %d", (int) (60.-toDegrees(minAngles[i])));
            else //quad
                fprintf(f," %d", (int) (90.-toDegrees(minAngles[i])));
        }
    }

    fclose(f);

    return true;
}

//-------------------------------------------------------------------
//-------------------------------------------------------------------
bool Services::WriteOFF(std::string name, FEMesh &output){

    const vector<Point3D> &points = output.getPoints();
    const vector<vector<unsigned int> > &elements = output.getElements();
    const vector<unsigned int> &colored = output.getColoredCells();

    if (elements.empty()) {
        std::cout << "no output elements\n";
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
bool Services::WriteMixedVolumeMesh(std::string name, FEMesh &output){

    const vector<Point3D> &points = output.getPoints();
    const vector<vector<unsigned int> > &elements = output.getElements();

    if (elements.empty()) {
        std::cout << "no output elements\n";
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
bool Services::WriteOutputMesh(std::string name, FEMesh &output){

    vector<Point3D> points = output.getPoints();
    vector<vector<unsigned int> > elements = output.getElements();

    if (elements.empty()) {
        std::cout << "no output elements\n";
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
bool Services::WriteMeshGetfem(std::string name, FEMesh &output){

    /* 2D file */

    vector<Point3D> points = output.getPoints();
    vector<vector<unsigned int> > elements = output.getElements();

    if (elements.empty()) {
        std::cout << "no output elements\n";
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

    unsigned int aux;

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
}

