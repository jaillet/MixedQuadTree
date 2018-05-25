cmake_minimum_required(VERSION 2.6)
project(mesher_roi)

#set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")

set(SOURCE_FILES Main MeshPoint Point3D TriMeshPolyline PolyEdge  Quadrant QuadEdge
	QuadTransition QuadSurfTemplate 
	HexRotation
        SurfTemplate1 SurfTemplate2 SurfTemplate3 SurfTemplate4
        SurfTemplate5 SurfTemplate6 SurfTemplate7
        Mesher GridMesher TransitionPoint PatternFactory TransitionTemplate
        TransitionCube SurfPair SurfEdgeContainer
        RefinementRegion RefinementCubeRegion RefinementSurfaceRegion
        RefinementInputSurfaceRegion RefinementAllRegion
        FaceDiagonal SurfPrism SurfPyramid SurfHexahedron SurfTetrahedron 
         GeometricTransform
        Visitors/Visitor Visitors/SplitVisitor Visitors/EdgeVisitor 
	Visitors/IntersectionsVisitor
        Visitors/OneIrregularVisitor Visitors/PointMovedVisitor
        Visitors/TransitionPatternVisitor Visitors/SurfaceTemplatesVisitor
        Visitors/RemoveSubElementsVisitor

        SurfEdge SurfTriangle )


add_executable(mesher_roi ${SOURCE_FILES})
