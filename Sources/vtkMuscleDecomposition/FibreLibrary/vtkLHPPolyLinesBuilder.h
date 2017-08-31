#pragma once

#ifndef __vtkLHPPolyLinesBuilder_h
#define __vtkLHPPolyLinesBuilder_h

#include "vtkObject.h"
#if VTK_MAJOR_VERSION < 5
#include "vtkPolyDataToPolyDataFilter.h"
#define vtkPolyDataAlgorithm vtkPolyDataToPolyDataFilter
#else
#include "vtkPolyDataAlgorithm.h"
#endif


class vtkPoints;
class vtkCellArray;

/** 
Constructs a polydata in which each cell represents one polyline.

If the input polydata do not have lines specified, i.e., it contains point sets,
the output polydata contains just one polyline formed by vertices with indices 
0, 1, ... N-1, where N is the number of points.

If the input polydata contains multiple lines whose indices of end-poins matches,
these lines are combined to form a single polyline. For example, cells (0,1), (1,2),
(2,3),(4,5),(5,6) will result in two polylines (0,1,2,3) and (4,5,6).

If the input polydata contains cells of different type, error is raised.

Data arrays associated with the points are passed, data arrays associated with cells
are passed only if there is no change, in which case everything is bypassed.
*/
class vtkLHPPolyLinesBuilder : public vtkPolyDataAlgorithm
{
public:
	vtkTypeMacro(vtkLHPPolyLinesBuilder, vtkPolyDataAlgorithm);

	/** Constructs a new analysis filter that calculates nothing by the default. */
	static vtkLHPPolyLinesBuilder *New();
			
protected:
	vtkLHPPolyLinesBuilder();
	~vtkLHPPolyLinesBuilder();

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
	/** Creates a new polyline containing numPoints points. */
	vtkCellArray* CreatePolyLine(int numPoints);

	/** Merges the passed polylines together to create a minimal set of polylines. 
	If the input polylines are already the minimal set, the method returns it.
	N.B. the caller must call Delete() method on the return cell array in both cases. */
	vtkCellArray* MergePolyLines(vtkCellArray* inputLines);

protected:	

};

#if VTK_MAJOR_VERSION < 5
#undef vtkPolyDataAlgorithm
#endif

#endif //vtkLHPPolyLinesBuilder
