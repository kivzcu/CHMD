cmake_minimum_required(VERSION 2.8.11)
project(FibreLibrary)

#Add support for Edit & Continue fo MSVC
IF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /ZI")
ENDIF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")

set(HDRS_FibreLibrary
				hungarian.h
				memdebug.h
				targetver.h
				vtkLHPMuscleFibresAnalysisFilter.h
				vtkLHPMuscleFibresMath.h
				vtkLHPMuscleFibresResample.h
				vtkLHPSetGet.h
				vtkMAFPolyDataCutOutFilterEx.h
				vtkMapperActor.h
				vtkMuscleDecomposer.h
				vtkMuscleDecomposerSplineInterpolator.h
				vtkMuscleDecomposerUtils.h
				vtkLHPPolyLinesBuilder.h
				vtkCatmullRomSpline.h
				vtkParametricCatmullRomSpline.h
				vtkParametricCatmullRomSurface.h
				vtkTriangleMeshLaplacianDeformation.h
				csparse.h
				vtkLHPSetGet.h
)

if (CHMD_BUILD_PRIVATE_APPS)
  list(APPEND  HDRS_FibreLibrary
				vtkCreateMTFiber.h
				vtkSortMTFiber.h
				vtkDirectionMTFiber.h
				vtkRemovePoints.h
				vtkDataLoader.h
				vtkInfo.h
				vtkNearestNeighbours.h
				vtkLinesRestorer.h
				vtkStdMcr.h
  )
endif()

set(SRCS_FibreLibrary
				hungarian.c
				vtkLHPMuscleFibresAnalysisFilter.cpp
				vtkLHPMuscleFibresMath.cpp
				vtkLHPMuscleFibresResample.cpp
				vtkMAFPolyDataCutOutFilterEx.cxx
				vtkMapperActor.cpp
				vtkMuscleDecomposer.cpp
				vtkMuscleDecomposerSplineInterpolator.cpp
				vtkMuscleDecomposerUtils.cpp
				vtkLHPPolyLinesBuilder.cpp
				vtkCatmullRomSpline.cxx
				vtkParametricCatmullRomSpline.cxx
				vtkParametricCatmullRomSurface.cxx
				vtkTriangleMeshLaplacianDeformation.cpp
				csparse.c
)

if (CHMD_BUILD_PRIVATE_APPS)
  list(APPEND  SRCS_FibreLibrary
				vtkCreateMTFiber.cxx
				vtkSortMTFiber.cxx
				vtkDirectionMTFiber.cxx
				vtkRemovePoints.cxx
				vtkDataLoader.cxx
				vtkInfo.cxx
				vtkLinesRestorer.cxx
				vtkNearestNeighbours.cxx
  )
endif()

if (CHMD_BUILD_PRIVATE_APPS)
  add_subdirectory(HaoLi)
endif()

add_library(FibreLibrary ${SRCS_FibreLibrary} ${HDRS_FibreLibrary})