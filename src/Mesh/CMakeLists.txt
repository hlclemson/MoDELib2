target_sources(${PROJECT_NAME} PRIVATE 
		${CMAKE_CURRENT_SOURCE_DIR}/BarycentricTraits.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/GmshReader.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/MeshPlane.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/MeshLine.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/PlanarMeshFace.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/MeshBoundarySegment.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/Simplex.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/MeshRegion.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/MeshRegionBoundary.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/triangle.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/SimplicialMesh.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/TriangularMesh.cpp
			)

