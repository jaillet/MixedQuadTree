cmake_minimum_required(VERSION 2.6)
project(mesher_roi)

#set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ltbb -std=c++11 -fopenmp")

if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  # Update if necessary
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -pedantic")
endif()

#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")

set(SOURCE_FILES Main
        Point3D Polyline PolyEdge  Quadrant QuadEdge
	QuadTransition QuadSurfTemplate 
        Mesher GridMesher MeshPoint
        RefinementRegion RefinementCubeRegion RefinementSurfaceRegion
        RefinementInputSurfaceRegion RefinementAllRegion RefinementBoundaryRegion
        GeometricTransform Services
        Visitors/Visitor Visitors/SplitVisitor Visitors/EdgeVisitor 
	Visitors/IntersectionsVisitor
        Visitors/OneIrregularVisitor Visitors/PointMovedVisitor
        Visitors/TransitionPatternVisitor Visitors/SurfaceTemplatesVisitor
        Visitors/RemoveSubElementsVisitor
        ParallelizeTest/CloneMesher/
        )

set(SOURCE_MEANTIME MeanTimePerformance/Main)
set(SOURCE_PARALLELIZE ParallelizeTest/ParallelizeTest.cpp ParallelizeTest/TimeDecorator.hpp)

add_executable(mesher_roi ${SOURCE_FILES})
add_executable(meanTime ${SOURCE_MEANTIME})
add_executable(parallelize_test ${SOURCE_PARALLELIZE})

# obsolete 3D files to be removed

# obsolete 3D files already been removed
#TransitionPoint TransitionCube
#FaceDiagonal Prism Octant OctreeEdge
#PatternFactory TransitionTemplate
#SurfTemplate1 SurfTemplate2 SurfTemplate3 SurfTemplate4
#SurfTemplate5 SurfTemplate6 SurfTemplate7 TriMesh HexRotation
#SurfPrism SurfPyramid SurfHexahedron SurfTetrahedron SurfEdge SurfTriangle SurfEdgeContainer SurfPair
