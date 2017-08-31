cmake_minimum_required(VERSION 2.8.11)

project(vtkMuscleDecomposition)

if(CMAKE_COMPILER_IS_GNUCXX)
	message(STATUS "GCC detected, adding compile flags")
	set(CMAKE_CXX_FLAGS "-std=c++11 -Wno-write-strings -Wno-format-security -m64 -g")
endif(CMAKE_COMPILER_IS_GNUCXX)

#Add support for Edit & Continue fo MSVC
IF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /ZI")
  SET(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} /SAFESEH:NO")
ENDIF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")


include_directories(${FibreLibrary_SOURCE_DIR})

set(SRCS_vtkMuscleDecomposition
				demoApp.cpp
				stdafx.cpp
)

add_executable(vtkMuscleDecomposition ${SRCS_vtkMuscleDecomposition})

target_link_libraries(vtkMuscleDecomposition LINK_PUBLIC FibreLibrary ${VTK_LIBRARIES})
