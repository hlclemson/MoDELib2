target_sources(${PROJECT_NAME} PRIVATE 
	${CMAKE_CURRENT_SOURCE_DIR}/GlidePlane.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GlidePlaneFactory.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/PeriodicGlidePlaneFactory.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/PeriodicGlidePlane.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GlidePlaneBase.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/NoiseTraits.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GlidePlaneNoiseBase.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/GlidePlaneNoise.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/AnalyticalSolidSolutionNoise.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/MDSolidSolutionNoise.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/MDStackingFaultNoise.cpp
			)
