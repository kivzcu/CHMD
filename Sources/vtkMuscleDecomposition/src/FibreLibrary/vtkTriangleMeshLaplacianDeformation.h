#pragma once

#ifndef __vtkLHPMuscleFibresResample_h
#define __vtkTriangleMeshLaplacianDeformation_h

#include "vtkObject.h"
#if VTK_MAJOR_VERSION < 5
#include "vtkPolyDataToPolyDataFilter.h"
#define vtkPolyDataAlgorithm vtkPolyDataToPolyDataFilter
#else
#include "vtkPolyDataAlgorithm.h"
#endif

#include <map>

class vtkIdList;
class vtkPoints;

/** Laplacian deformation of triangular meshes. */
class vtkTriangleMeshLaplacianDeformation : public vtkPolyDataAlgorithm
{
public:
	vtkTypeMacro(vtkTriangleMeshLaplacianDeformation, vtkPolyDataAlgorithm);

	/** Constructs a new analysis filter that calculates nothing by the default. */
	static vtkTriangleMeshLaplacianDeformation *New();	
	
	/** Gets/Sets the indices of the fixed points */
	vtkGetObjectMacro(FixedPointsIds, vtkIdList);
	vtkSetObjectMacro(FixedPointsIds, vtkIdList);

	/** Gets/Sets the indices of the moved points */
	vtkGetObjectMacro(MovedPointsIds, vtkIdList);
	vtkSetObjectMacro(MovedPointsIds, vtkIdList);
	
	/** Gets/Sets the coordinates of the moved points */
	vtkGetObjectMacro(MovedPointsCoords, vtkPoints);
	vtkSetObjectMacro(MovedPointsCoords, vtkPoints);
	

protected:
	vtkTriangleMeshLaplacianDeformation();
	~vtkTriangleMeshLaplacianDeformation();

protected:
#if VTK_MAJOR_VERSION < 5
	/**
	This method is the one that should be used by subclasses, right now the
	default implementation is to call the backwards compatibility method */
	/*virtual*/void ExecuteData(vtkDataObject *output);
#else
	// Description:	
	// This is the method of the algorithm in which the algorithm should fill in the output ports
	/*virtual*/ int RequestData(vtkInformation* request,
		vtkInformationVector** inputVector,	vtkInformationVector* outputVector);
#endif

protected:
	struct VCoord
	{
		double xyz[3];
	};
	

	typedef std::map< vtkIdType, double> RowMap;

	/** Constructs sparse Laplacian matrix for the input triangular mesh*/
	void ConstructMatrixL(vtkPolyData* inputPoly, std::vector< RowMap >& matrixL) const;

	/** Constructs the vector B = A^T*b where b is n+m+o x 1 vector of right side 
	N.B. this method assumes matrixL is correctly filled by ConstructMatrixL*/
	void ConstructVectorB(vtkPolyData* inputPoly, const std::vector< RowMap >& matrixL, std::vector< VCoord >& vectorB) const;

	/** Constructs sparse matrix Q = A^T*A
	N.B. this method assumes matrixL is correctly filled by ConstructMatrixL */
	void ConstructMatrixAtA(vtkPolyData* inputPoly, const std::vector< RowMap >& matrixL, std::vector< RowMap >& matrixQ) const;
	
protected:
	vtkIdList* FixedPointsIds;
	vtkIdList* MovedPointsIds;
	vtkPoints* MovedPointsCoords;
	
};

#if VTK_MAJOR_VERSION < 5
#undef vtkPolyDataAlgorithm
#endif

#endif //vtkTriangleMeshLaplacianDeformation
