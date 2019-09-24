## refineMeshTBB.hpp

* function refineMeshParallelTest1TBB : version without reduction using tbb and the container concurrent_vector for new_points and a concurrent_unordered_set (with a custom hash) for QuadEges

* Versions with a reduction
	* function refineMeshReductionTBB : version with reduction using tbb::parallel_reduce
	* function refineMeshCustomReductionTBB : use a vector of RefineMeshReduction objects (class that do the parallel task), and do the join manually.
	* function refineCustomMeshReductionTBBV2 : improve of previous version, with a master join

## refineMeshOpenMP.hpp

* function refineMeshParallelOpenMP : basic parallelization with mutex (critical region) and all datasets shared
* function : implementation of a parallel reduce in openMP, use a custom function to make the reduce
* function make_reduceV1 : first version of reducing
* function make_reduceV2 : second version -> TODO !!  

## splitVisitorOpenMP.cpp / .h

Visitor with mutex (critical region in OpenMP) used in function refineMeshParallelOpenMP

## splitVisitorReductionOpenMP.cpp / .h

The same visitor than the original, but edges is split in 2 dataset : one for reading only (current_edges) and one for new edges (new_edges)

