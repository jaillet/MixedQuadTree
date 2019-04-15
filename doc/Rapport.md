# Productions

* Analysis scripts : `script/analyse_time.sh̀̀̀`, `script/memory_usage.sh̀̀̀` et `script/massif_analyser.awk`
* Grph plot : `script/plot_time.py` 
* Test program : `build/parallelize_test`

# Comparison 

## Open MP

### Pros
* Easy to use
 
### Cons
* No support of iterators

## Intel TBB

### Pros
* Container thread-safe
* Tasks

### Cons

# Comparaisons de consommation mémoire

## Script for memory analysis

Program `script/memory_usage.sh [N]`

Uses valgrind 


Execute N times the program `mesher_roi` with a.poly as input, and option -s i (with i from 1 to N)
Save the result in `build/memory_usage_mesher_roi_N_DATE`
Files massif.out.\*.i are the detailed results of the memory analysis, i is the level of refinement.
The file memory_usage contains the peak memory usage for each level of refinement.

## Use of deque instead of list

Comparison between vector / list / deque : https://baptiste-wicht.com/posts/2012/12/cpp-benchmark-vector-list-deque.html  

Deque uses slightly more memory.

## Memory usage with list

Command : mesher_roi -p ../data/a.poly -s N

|N  | Peak memory usage|
|---|------------------|
|1  | 115.816 Ko|
|2  | 129.616 Ko|  
|3  | 166.096 Ko|  
|4  | 265.664 Ko|  
|5  | 534.928 Ko|  
|6  | 1.18233 Mo|  
|7  | 2.57091 Mo|  
|8  | 5.66521 Mo|  
|9  | 12.0442 Mo|  
|10 | 24.0443 Mo|  


Command : mesher_roi -p ../data/a.poly -a N

|N  | Peak memory usage|
|---|------------------|
1  | 115.816 Ko
2  | 129.616 Ko
3  | 166.28 Ko
4  | 265.56 Ko
5  | 564.976 Ko
6  | 1.55675 Mo
7  | 5.06608 Mo
8  | 17.9289 Mo
9  | 68.0225 Mo
10 | 263.253 Mo


## Memory usage with deque


Command : mesher_roi -p ../data/a.poly -s N

|N  | Peak memory usage|
|---|------------------|
|1  | 119.712 Ko|
|2  | 134.368 Ko|
|3  | 169.856 Ko|
|4  | 264.68 Ko|
|5  | 543.824 Ko|
|6  | 1.23563 Mo|
|7  | 2.75078 Mo|
|8  | 6.0956 Mo|
|9  | 12.9715 Mo|
|10 | 25.9963 Mo|

Command : mesher_roi -p ../data/a.poly -a N

|N  | Peak memory usage|
|---|------------------|
1  | 119.792 Ko
2  | 134.368 Ko
3  | 169.968 Ko
4  | 264.68 Ko
5  | 570.192 Ko
6  | 1.57748 Mo
7  | 5.12051 Mo
8  | 18.2691 Mo
9  | 68.9012 Mo
10 | 267.182 Mo


# Comparision of execution time

## Script for the analysis

Program `script/analyse_time.sh`

Execute the program `build/parallelize_test` for different number of elements (100 to 10^9), and different number of threads (8 to 1).
Save the result in `script/analyse_time_[DATE]/time_size_[NBELEM]_thread_[NBTRHEAD]`.

## Plot the result

Execute `python3 script/plot_time.py` (pyplot)

3 types of graph :

* elements_times : For each method, a graph that shows the execution time per number of elements, for each number of threads.

* thread_times : For each number of elements, a graph that shows the execution time per number of threads, for each method.

* strong_scaling : Same as thread_times, but with speedup (TimeWith1Thread divided by TimeWithNThread)

# Mesher_roi analysis

## Mesher class

### Member variables

```cpp
vector<MeshPoint> points;
vector<Quadrant> Quadrants;
set<QuadEdge> QuadEdges;
list<RefinementRegion *> regions;
```

### Function generateMesh (or refineMesh)

* generateMesh : Refinement level start at 0.
* refineMesh : Refinement level start at the refinement level of the input.

#### Parameters

	- input : Reference on the Polyline
	- rl : Refinement level
	- name : output name
	- all_reg : list of the refinement region (Surface and/or allRegion and/or Boundary(for refineMesh))
	- decoration : boolean to add decoration to the VTK
	
Additional parameters for refineMesh :

	- roctli : list of starting quadrants
	- gt : geometric transform
	- minrl : minimum refinement level of input quadrants.
	- omaxrl : maximum refinement level of input quadrants.

#### Algorithm

* rotateGridMesh(input, all_reg, gt) : If needed, apply gt on input
* generateGridMesh(input) : Create and store in `Quadrants` the starting Quadrants (from the polyline)
* generateQuadTreeMesh(rl, input, all_reg, name) : Algorithm for mesh generation (For `refineMesh`, call splitQuadrants)

The rest is identical for refineMesh and generateMesh.

### Function generateQuatdreeMesh

#### Parameters
    - rl : refinement level that we want
    - input : The polyline
    - all_reg : list of RefinementRegion
    - name : output name

#### Algo 1 (handle the boundary)

##### Init
        - temp_Quadrants : list of actual Quadrants
        - new_Quadrants : empty list of Quadrants
        - new_pts : empty list of Point3D
        - sv : SplitVisitor
        - i : (output) current level of refinement (start at 0)

##### First level

<pre>
Do  
    Empty <b>new_pts</b>  
    While <b>tmp_Quadrants</b> is not empty 
        Level 2  
    Swap <b>tmp_Quadrants</b> and <b>new_Quadrants</b>
    If <b>new_pts</b> empty
        break
    Add all points in <b>new_pts</b> to <b>points (Mesher attribute)</b> 
    Increment <b>i</b>
While <b>new_pts</b> is not empty 
</pre>

##### Second level

<pre>
Tant que <b>tmp_Quadrants</b> n'est pas vide
    <b>iter</b> <- début d'itérateur de <b>tmp_Quadrants</b>
    
    computeMaxDistance sur <b>iter</b> (Quadrant)
        Lecture <b>points</b> (var classe Mesher)
        Lecture <b>pointindex</b> (var classe Quadrant)
        Lecture <b>sub_elements</b> (var classe Quadrant)
        Ecriture <span style="color:orange"><b>max_dis</b></span> (var classe Quadrant)
        
    <b>to_refine</b> <- (*all_reg.begin())->intersectsQuadrant(points, *iter)
        TODO
    
    S'il n'est pas <b>to_refine</b>
        Ajout de <b>iter</b> dans <b>new_quadrants</b>
    Sinon
        Init ref <b>inter_edges</b> avec <b>iter->intersected_edges</b> (var classe Quadrant)
        Init <b>clipping_coords</b> et set à <b>sv.clipping</b> (SplitVisitor)
        Init <b>split_elements</b> et set à <b>sv.new_eles</b> (SplitVisitor)
        
        Applique <b>sv</b> (SplitVisitor) sur <b>iter</b> (Quadrant)
            Lecture <b>iter.pointindex</b>
            Insertion <b>sv.new_eles</b> (ref vers <b>split_elements</b>)
            Lecture <b>sv.points</b> (ref vers <b>points</b> var classe Mesher)          
            Insertion / Lecture <b>sv.new_pts</b> (ref vers <b>new_pts</b>) 
            Insertion / Suppression / Lecture <b>sv.edges</b> (ref vers <b>QuadEdges</b> var classe Mesher)         
            Insertion <b>sv.clipping</b> (ref vers <b>clipping_coords</b>)
        
        Si <b>inter_edges</b> (ref vers <b>iter->inteersected_edges</b>) vide
            Pour tout <b>split_elements</b>
                Init Quadrant <b>o</b> à partir de <b>split_elements</b>
                Insertion <b>new_Quadrants</b> du quadrant <b>o</b> 
        Sinon
            Pour tout <b>split_elements</b>
                Init Quadrant <b>o</b> à partir de <b>split_elements</b>
                Init iv (IntersectionVisitor)
                Set <b>iv.ply</b> (ref vers <b>input</b>)
                Set <b>iv.edges</b> (ref vers <b>inter_edges</b>)
                Set <b>iv.coords</b> (ref vers <b>clipping_coords</b> du split element)
                
                Applique <b>iv</b> (IntersectionVisitor) sur <b>iter</b> (Quadrant)
                    Insertion / Lecture <b>o.intersected_edges</b>
                    Lecture <b>iv.ply.mVertices</b> (ref vers input.mVertices)
                    Lecture <b>iv.edges</b> (ref vers iter->intersected_edges)
                    Lecture <b>iv.ply.mEdges</b> (ref vers input.mEdges)
                    
                    call intersectsEdge
                        Lecture <b>iv.coords</b> (ref vers <b>clipping_coords</b> du split element)
                        Lecture <b>iv.ply.mVertices</b> (ref vers input.mVertices)
                        Lecture <b>iv.ply.mEdges</b> (ref vers input.mEdges)
                        
                        call computePosition
                            Lecture <b>iv.ply.mEdges</b> (ref vers input.mEdges)
                        
                
                        
                        TODO cpliGeneralCase
                    
                Si retour fonction vrai
                    Insertion <b>new_Quadrants</b> du quadrant <b>o</b>  
                Sinon
                    Apelle fonction isItIn
                        TODO
                                            
                    Si retour fonction vrai
                        Insertion <b>new_Quadrants</b> du quadrant <b>o</b>  
                
    Retire <b>iter</b> de <b>tmp_quadrants</b>  
    
</pre>

#### Algo 2

##### First level

<pre>
For each level until rl  
    Clear <b>new_pts</b>  
    While until <b>tmp_Quadrants</b> is empty (go to Second level)
    Swap <b>tmp_Quadrants</b> and <b>new_Quadrants</b>
    If <b>new_pts</b> empty
        break
    Insert in <b>points (var class Mesher)</b> the content of <b>new_pts</b>
</pre>

##### Second level

This is what we need to parallelize first

<pre>
While until <b>tmp_Quadrants</b> is empty
    <b>iter</b> <- begin iterator of <b>tmp_Quadrants</b>
    
    For each RafinementRegion in <b>all_reg</b>
        <b>region_rl</b> <- getRafinementLevel of RafinementRegion
        If <b>region_rl</b> is lower than <b>i</b> (raffinement level) then continue
        If <b>region_rl</b> is equal or lower than <b>iter</b>  getRafinementLevel then continue
        If <b>reg_iter</b> intersectsQuantrant on points en <b>iter</b> then <b>to_refine</b> <- true
        
    If not <b>to_refine</b>
        Insert <b>iter</b> in <b>new_quadrants</b>
    Else (to refine)
        Init ref <b>inter_edges</b> with <b>iter->intersected_edges</b> (var class Quadrant)
        Init <b>qrl</b> with <b>iter</b> rafinement level
        Init <b>clipping_coords</b> and set <b>sv.clipping</b> (SplitVisitor)
        Init <b>split_elements</b> and set <b>sv.new_eles</b> (SplitVisitor)
        
        Apply <b>sv</b> (SplitVisitor) on <b>iter</b> (Quadrant)
            Read <b>iter.pointindex</b>
            Insert <b>sv.new_eles</b> (ref <b>split_elements</b>)
            Read <b>sv.points</b> (ref <b>points</b> var class Mesher)          
            Insert / Read <b>sv.new_pts</b> (ref <b>new_pts</b>) 
            Insert / Remove / Read <b>sv.edges</b> (ref <b>QuadEdges</b> var class Mesher)         
            Insert <b>sv.clipping</b> (ref <b>clipping_coords</b>)
        
        If <b>inter_edges</b> (ref <b>iter->inteersected_edges</b>) empty
            For each <b>split_elements</b>
                Init Quadrant <b>o</b> based on <b>split_elements</b>
                Insert <b>o</b> quadrant in <b>new_Quadrants</b>
        Else
            For each <b>split_elements</b>
                Init Quadrant <b>o</b> based on <b>split_elements</b>
                Init iv (IntersectionVisitor)
                Set <b>iv.ply</b> (ref <b>input</b>)
                Set <b>iv.edges</b> (ref <b>inter_edges</b>)
                Set <b>iv.coords</b> (ref <b>clipping_coords</b> of split element)
                
                Apply <b>iv</b> (IntersectionVisitor) on <b>iter</b> (Quadrant)
                    Insert / Read <b>o.intersected_edges</b>
                    Read <b>iv.ply.mVertices</b> (ref input.mVertices)
                    Read <b>iv.edges</b> (ref iter->intersected_edges)
                    Read <b>iv.ply.mEdges</b> (ref input.mEdges)
                    
                    Call intersectsEdge
                        Read <b>iv.coords</b> (ref <b>clipping_coords</b> of split element)
                        Read <b>iv.ply.mVertices</b> (ref input.mVertices)
                        Read <b>iv.ply.mEdges</b> (ref input.mEdges)
                        
                        Call computePosition
                            Read <b>iv.ply.mEdges</b> (ref input.mEdges)
                        
                        Call cpliGeneralCase
                            TODO 
                    
                If return true
                    Insert <b>o</b> quadrant in <b>new_Quadrants</b>
                Else
                    Call fonction isItIn
                        TODO
                                            
                    If return true
                        Insert <b>o</b> quadrant in <b>new_Quadrants</b>  
                
    Remove <b>iter</b> from <b>tmp_quadrants</b>  
    
</pre>

##### Concurrent access for second level

Here are the concurrent access for one refinement level :


At the beginning, to know the boolean value of to_refine, virtual function RefinementRegion::intersectsQuadrant(points, iter) is called.
These functions read <b>points</b>, <b>polyline</b>, and Read and can Modify the intersected edge of <b>Quadrant</b>. with function Polyline::getNbFeatures in RefinementboundaryRegion. But this information is not used by other thread (Is this sure ?)

Insert in <b>new_quadrants</b> (output of algo 2)

SplitVisitor (set class var)

Insert / Read <b>sv.new_pts</b> (ref <b>new_pts</b>) 
Insert / Remove / Read <b>sv.edges</b> (ref <b>QuadEdges</b> member variable of class Mesher)  

Remove <b>iter</b> from <b>tmp_quadrants</b> PAUL : where ?

In SplitVisitor, insert new points in <b>new_pts</b> and use the indice of this points to update mid point of edges.  
So, this should be fixed for multi-threading because multiple thread can insert points at the same time.

# Parallel version

A parallel version using tbb and the container concurrent_vector for new_points and a concurrent_unordered_set for QuadEges is implemented.

## Ideas

If the refinement region is a box, or only the quadrants that intersect the polyline, the load balance will not be fair, so maybe add a task only if the quadrant need to be refined. 


# Parallelism with a reduction

To avoid thinking about concurrence issues like inserting new edges, creating new points with good index, we thought about a reduction, ie every thread has private memory, and a single thread join all previous thread. That means that an algorithm need to be implemented to attribute new mesh points indexes.

## Algo

Each thread will have it's own copy of <b>new_pts</b>, <b>new_Quadrants</b> but they will share quadEdges as a concurrent set (because the need it for reading, inserting and removing)


