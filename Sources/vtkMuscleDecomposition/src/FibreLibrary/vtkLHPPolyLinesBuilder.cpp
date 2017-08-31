/*=========================================================================
Program: Musculoskeletal Modeling (VPHOP WP10)
Module: vtkLHPPolyLinesBuilder.cpp

Authors: Josef Kohout
==========================================================================
Copyright (c) 2014 University of West Bohemia (www.zcu.cz)
See the COPYINGS file for license details
=========================================================================
*/
#include "vtkLHPPolyLinesBuilder.h"
#include "vtkObjectFactory.h"
#include "vtkPolyData.h"
#include "vtkCellArray.h"
#include "vtkMath.h"
#include "vtkLHPMuscleFibresMath.h"
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <unordered_map>
#include <list>

vtkStandardNewMacro(vtkLHPPolyLinesBuilder);


vtkLHPPolyLinesBuilder::vtkLHPPolyLinesBuilder()
{
	
}


vtkLHPPolyLinesBuilder::~vtkLHPPolyLinesBuilder()
{
	
}

#if VTK_MAJOR_VERSION < 5
#define VTK6_RESULT_OK
#define VTK6_RESULT_FAIL
/**
Executes the data operation.

@param [in,out]	output	If non-null, the output.
*/
/*virtual*/ void vtkLHPPolyLinesBuilder::ExecuteData(vtkDataObject *output)
{
	vtkPolyData* inputPoly = GetInput();
	vtkPolyData* outputPoly = vtkPolyData::SafeDownCast(output);
#else

#define VTK6_RESULT_OK		1
#define VTK6_RESULT_FAIL	0

// This is the method of the algorithm in which the algorithm should fill in the output ports
/*virtual*/ int vtkLHPPolyLinesBuilder::RequestData(vtkInformation* request,
	vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
	vtkPolyData* inputPoly = vtkPolyData::GetData(inputVector[0]);
	vtkPolyData* outputPoly = vtkPolyData::GetData(outputVector);
#endif

	//check whether output is valid	
	if (inputPoly == NULL || inputPoly->GetPoints() == NULL ||
		inputPoly->GetPoints()->GetNumberOfPoints() == 0)
	{
		vtkErrorMacro(<< "Invalid input for vtkLHPPolyLinesBuilder.");
		return VTK6_RESULT_FAIL;   //we have no valid input
	}
	
	if (outputPoly == NULL)
	{
		vtkErrorMacro(<< "Invalid output for vtkLHPPolyLinesBuilder.");
		return VTK6_RESULT_FAIL;   //we have no valid output
	}

	//pass points
	outputPoly->SetPoints(inputPoly->GetPoints());

	//pass points data
	outputPoly->GetPointData()->PassData(inputPoly->GetPointData());
		
	vtkCellArray* lines = NULL;
	vtkCellArray* inputLines = inputPoly->GetLines();
	if (inputLines == NULL || inputLines->GetNumberOfCells() == 0) {
		lines = CreatePolyLine(outputPoly->GetNumberOfPoints());
	}
	else {
		lines = MergePolyLines(inputLines);
	}

	outputPoly->SetLines(lines);
	lines->Delete();

	//if nothing has changed, pass cell data
	if (inputLines == lines)
		outputPoly->GetCellData()->PassData(inputPoly->GetCellData());

	return VTK6_RESULT_OK;
}

//Creates a new polyline containing numPoints points.
vtkCellArray* vtkLHPPolyLinesBuilder::CreatePolyLine(int numPoints)
{
	//if there are no cells, points should form a single fibre
	vtkCellArray* lines = vtkCellArray::New();	
	lines->InsertNextCell(numPoints);

	for (int i = 0; i < numPoints; i++) {
		lines->InsertCellPoint(i);
	}

	return lines;
}

//Merges the passed polylines together to create a minimal set of polylines.
//If the input polylines are already the minimal set, the method returns it.
//N.B.the caller must call Delete() method on the return cell array in both cases
vtkCellArray* vtkLHPPolyLinesBuilder::MergePolyLines(vtkCellArray* inputLines)
{
	int nCells = inputLines->GetNumberOfCells();
	vtkIdType* pIdsPtr = inputLines->GetPointer();
	
	struct POLYLINE
	{		
		int nPoints;	//number of points forming the polyline
		std::list< vtkIdType* > segments;	//list of segments forming the polyline
	};

	POLYLINE* pols = new POLYLINE[nCells];	//there will be at most as many polylines as in the input
	int nPolylines = 0;	

	//maps end-points to polylines
	std::unordered_map<vtkIdType, POLYLINE*> mapBeg;
	std::unordered_map<vtkIdType, POLYLINE*> mapEnd;

	POLYLINE* pol;	//current polyline
	for (int i = 0; i < nCells; i++)
	{	
		//if we denote the current polyline B, there are three options: 
		//there is polyline A such that A.begin = B.end => prepend B to A (B->A)
		//there is polyline A such that A.end = B.begin => append B to A  (A->B)
		//otherwise B is a new independent polyline
		
		auto itBeg = mapBeg.find(pIdsPtr[*pIdsPtr]);
		if (itBeg != mapBeg.end())
		{
			//case 1: prepend B to A
			pol = itBeg->second;
			mapBeg.erase(itBeg);
			
			pol->segments.push_front(pIdsPtr);						
			mapBeg[pIdsPtr[1]] = pol;
		}
		else
		{
			auto itEnd = mapEnd.find(pIdsPtr[1]);
			if (itEnd != mapEnd.end())
			{
				//case 2: we will append the current polyline 
				pol = itEnd->second;
				mapEnd.erase(itEnd);
			}
			else
			{
				//case 3: start a new POLYLINE
				pol = &pols[nPolylines];
				nPolylines++;

				pol->nPoints = 1;

				//update the beginning
				mapBeg[pIdsPtr[1]] = pol;
			}
			
			pol->segments.push_back(pIdsPtr);			
			mapEnd[pIdsPtr[*pIdsPtr]] = pol;
		}	
		
		pol->nPoints += (*pIdsPtr) - 1;	//increase the number of points
		pIdsPtr += (*pIdsPtr) + 1;
	}

	vtkCellArray* ret = NULL;
	if (nPolylines == nCells)
	{
		//nothing to merge done
		ret = inputLines;
		ret->Register(this);
	}
	else
	{
		//we need to merge cells
		ret = vtkCellArray::New();
		
		//first get the number of points we will have
		int nPoints = 0;
		for (int i = 0; i < nPolylines; i++){
			nPoints += pols[i].nPoints;
		}

		vtkIdTypeArray* pIds = vtkIdTypeArray::New();
		pIdsPtr = pIds->WritePointer(0, nPoints + nPolylines);

		for (int i = 0; i < nPolylines; i++)
		{
			*pIdsPtr = pols[i].nPoints; //write number of points per cell
			pIdsPtr++;

			*pIdsPtr = (*pols[i].segments.cbegin())[1];	//first point
			pIdsPtr++;

			//copy the data from every segment
			for (auto it = pols[i].segments.cbegin(); it != pols[i].segments.cend(); it++)
			{				
				vtkIdType* src = *it;	//do not copy the first (shared) point
				memcpy(pIdsPtr, &src[2], (*src - 1) * sizeof(vtkIdType));

				pIdsPtr += (*src - 1);
			}
		}		

		ret->SetCells(nPolylines, pIds);
		pIds->Delete();
	}

	delete[] pols;
	return ret;
}
