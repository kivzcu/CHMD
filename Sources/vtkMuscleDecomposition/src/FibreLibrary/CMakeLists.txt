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
				vtkCreateMTFiber.h
				vtkSortMTFiber.h
				vtkDirectionMTFiber.h
				vtkRemovePoints.h
				vtkDataLoader.h
				vtkInfo.h
				vtkLHPMuscleFibresAnalysisFilter.h
				vtkLHPMuscleFibresMath.h
				vtkLHPMuscleFibresResample.h
				vtkLHPSetGet.h
				vtkMAFPolyDataCutOutFilterEx.h
				vtkMapperActor.h
				vtkMuscleDecomposer.h
				vtkMuscleDecomposerSplineInterpolator.h
				vtkMuscleDecomposerUtils.h
				vtkStdMcr.h
				vtkLHPPolyLinesBuilder.h
				vtkCatmullRomSpline.h
				vtkParametricCatmullRomSpline.h
				vtkParametricCatmullRomSurface.h
				vtkNearestNeighbours.h
				vtkLinesRestorer.h
				vtkTriangleMeshLaplacianDeformation.h
				csparse.h
)

set(SRCS_FibreLibrary
				hungarian.c
				vtkCreateMTFiber.cxx
				vtkSortMTFiber.cxx
				vtkDirectionMTFiber.cxx
				vtkRemovePoints.cxx
				vtkDataLoader.cxx
				vtkInfo.cxx
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
				vtkNearestNeighbours.cxx
				vtkLinesRestorer.cxx
				vtkTriangleMeshLaplacianDeformation.cpp
				csparse.c
)

add_subdirectory(HaoLi)
add_library(FibreLibrary ${SRCS_FibreLibrary} ${HDRS_FibreLibrary})

