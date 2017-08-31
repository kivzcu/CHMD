cmake_minimum_required(VERSION 2.8.11)
project(MuscleDecomposition)

if(CMAKE_COMPILER_IS_GNUCXX)
	message(STATUS "GCC detected, adding compile flags")
	set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -Wextra -Wno-delete-non-virtual-dtor -Wno-unknown-pragmas -Wno-variadic-macros -Wno-long-long -Wno-format-security -Wno-write-strings -m64 -O3")
endif(CMAKE_COMPILER_IS_GNUCXX)

if(EXISTS "../src/vtk_dir")
	set(VTK_DIR "vtk_dir")
endif()

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
include_directories(SYSTEM ${VTK_INCLUDE_DIRS})

if (CHMD_BUILD_PRIVATE_APPS)
  find_package(Ceres REQUIRED)
  include_directories(SYSTEM ${CERES_INCLUDE_DIRS})
endif()

add_subdirectory(FibreLibrary)
if (CHMD_BUILD_PRIVATE_APPS)
  add_subdirectory(FibreOperation)
endif()
add_subdirectory(vtkMuscleDecomposition)

if (CHMD_BUILD_PRIVATE_APPS)
  add_subdirectory(ValidationDataGenerator)
endif()
