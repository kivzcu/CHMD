/*========================================================================= 
Program: Multimod Application Framework RELOADED 
Module: $RCSfile: vtkMAFPolyDataCutOutFilterEx.cxx,v $ 
Language: C++ 
Date: $Date: 2012-03-19 12:45:22 $ 
Version: $Revision: 1.1.2.2 $ 
Authors: Josef Kohout
========================================================================== 
Copyright (c) 2008 University of Bedfordshire (www.beds.ac.uk)
Copyright (c) 2011University of West Bohemia (www.zcu.cz)
See the COPYINGS file for license details 
=========================================================================
*/

//----------------------------------------------------------------------------
// Include:
//----------------------------------------------------------------------------
#include "vtkMAFPolyDataCutOutFilterEx.h"
//#include "vtkMAFVisualDebugger.h"
#include "vtkObjectFactory.h"
#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkCellArray.h"
#include "vtkMath.h"
#include "vtkSmartPointer.h"
#include <float.h>
#include <queue>
#include <stack>
#include <map>

vtkStandardNewMacro(vtkMAFPolyDataCutOutFilterEx);

//#include "mafMemDbg.h"
//#include "mafDbg.h"

#include "memdebug.h"

//crtdbg.h unix wrapper
#ifdef __unix__
	#include <assert.h>
	#define _ASSERT(...) assert(__VA_ARGS__)
	#define _ASSERTE(...) assert(__VA_ARGS__)
#endif

//----------------------------------------------------------------------------
vtkMAFPolyDataCutOutFilterEx::vtkMAFPolyDataCutOutFilterEx()
	//----------------------------------------------------------------------------
{
	CuttingPolyline = vtkPolyData::New();
	PointTolerance = 0.05;	//5% of edge length
	CutOutSide = 0;

	OutputCuttingPolyline = NULL;
	OutputClippedPart = NULL;

	this->Debug = 0;
	
	m_LastExecuteTimeStamp = 0;
	m_nCurrentMark = 0;

#ifdef DIJKSTRA_CUT
	m_nCurrentDijkstraMark = 0;
#endif
}

//----------------------------------------------------------------------------
vtkMAFPolyDataCutOutFilterEx::~vtkMAFPolyDataCutOutFilterEx()
	//----------------------------------------------------------------------------
{
	if (CuttingPolyline != NULL)
	{
		CuttingPolyline->Delete();
		CuttingPolyline = NULL;
	}

	if (OutputClippedPart != NULL)
	{
		OutputClippedPart->Delete();
		OutputClippedPart = NULL;
	}
	
	if (OutputCuttingPolyline != NULL)
	{
		OutputCuttingPolyline->Delete();
		OutputCuttingPolyline = NULL;
	}
}

//------------------------------------------------------------------------
//Return this object's modified time.
/*virtual*/ vtkMTimeType vtkMAFPolyDataCutOutFilterEx::GetMTime()
	//------------------------------------------------------------------------
{
	unsigned long mtime = Superclass::GetMTime();  
	unsigned long t1 = CuttingPolyline->GetMTime();

	return (t1 > mtime) ? t1 : mtime;        
}

#if VTK_MAJOR_VERSION < 5
#define VTK6_RESULT_OK
#define VTK6_RESULT_FAIL
/**
Executes the data operation.

@param [in,out]	output	If non-null, the output.
*/
/*virtual*/ void vtkMAFPolyDataCutOutFilterEx::ExecuteData(vtkDataObject *output)
{
	vtkPolyData* inputPoly = GetInput();
	vtkPolyData* outputPoly = vtkPolyData::SafeDownCast(output);
#else

#define VTK6_RESULT_OK		1
#define VTK6_RESULT_FAIL	0

// This is the method of the algorithm in which the algorithm should fill in the output ports
/*virtual*/ int vtkMAFPolyDataCutOutFilterEx::RequestData(vtkInformation* request,
	vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
	vtkPolyData* inputPoly = vtkPolyData::GetData(inputVector[0]);
	vtkPolyData* outputPoly = vtkPolyData::GetData(outputVector);
#endif

	unsigned long mtime = GetMTime();
	if (mtime <= m_LastExecuteTimeStamp)
		return VTK6_RESULT_OK; //no change

	//check whether output is valid
	if (inputPoly == NULL || inputPoly->GetPoints() == NULL || 
		inputPoly->GetPoints()->GetNumberOfPoints() == 0)    
	{
		vtkErrorMacro(<< "Invalid input for vtkMAFPolyDataCutOutFilterEx.");
		return VTK6_RESULT_FAIL;   //we have no valid input
	}

	if (this->CuttingPolyline == NULL || this->CuttingPolyline->GetPoints() == NULL || 
		this->CuttingPolyline->GetPoints()->GetNumberOfPoints() == 0)  
	{
		vtkErrorMacro(<< "Invalid CuttingPolyline for vtkMAFPolyDataCutOutFilterEx.");
		return VTK6_RESULT_FAIL;   //we have no valid input
	}

	if (outputPoly == NULL)
	{
		vtkErrorMacro(<< "Invalid output for vtkMAFPolyDataCutOutFilterEx.");
		return VTK6_RESULT_FAIL;   //we have no valid output
	}

	vtkIdList* pInsIds = vtkIdList::New();

	//performs the cutting
	InitMesh(inputPoly);

	//insert cutting points            
	if (InsertCutPoints(pInsIds)){				
		MarkTriangles(pInsIds);	//mark triangles to be saved
	}
	else
	{    
		m_nCurrentMark++;   //to force DoneMesh to produce empty dataset   
	}    
	
	DoneMesh(outputPoly);  

	//do we have additional outputs?
	//check, if the clipped part should be also preserved
	if (this->OutputClippedPart != NULL)
	{
		//invert marks
		int nOldCutside = this->CutOutSide;
		this->CutOutSide = 1 - this->CutOutSide;	//revert the side

		MarkTriangles(pInsIds);

		this->CutOutSide = nOldCutside;		

		DoneMesh(this->OutputClippedPart);
	}
	
	//do we need to save also cutting polyline?
	if (this->OutputCuttingPolyline != 0) {
		DoneCuttingPolyline(this->OutputCuttingPolyline, pInsIds);		
	}
	
	pInsIds->Delete();
	m_LastExecuteTimeStamp = mtime;   //new output available
	return VTK6_RESULT_OK;
}

//------------------------------------------------------------------------
//Creates internal structures for mesh
/*virtual*/ void vtkMAFPolyDataCutOutFilterEx::InitMesh(vtkPolyData* input)
	//------------------------------------------------------------------------
{	 
	int nNumOfTriangle = input->GetNumberOfCells();
	int nNumOfVertex = input->GetNumberOfPoints();

	m_Vertices.resize(nNumOfVertex);
	m_Triangles.resize(nNumOfTriangle);

	//load vertices
	for (int i = 0;i < nNumOfVertex; i++)
	{
		VERTEX v;    
		input->GetPoint(i, v.dCoord);
		v.nMark = -1;
#ifdef DIJKSTRA_CUT
		v.Dijstra.nDMark = 0;
#endif

		m_Vertices[i] = v;
	}

	int nEdges = 0;

	//load triangles   
	input->BuildCells();
	input->BuildLinks();
	vtkIdList* pIdList = vtkIdList::New();
	for(int i = 0; i < nNumOfTriangle; i++)
	{
		vtkIdType nPts, *pPts;
		input->GetCellPoints(i, nPts, pPts);

		if(!(nPts == 3)) continue;

		TRIANGLE t;
		for (int j = 0; j < 3; j++)
		{
			t.aVertex[j] = pPts[j];

			//get neighbours
			input->GetCellEdgeNeighbors(i, pPts[j], pPts[(j + 1) % 3], pIdList); 
			if (pIdList->GetNumberOfIds() == 0)
				t.pnb[j] = -1;  //isolated triangle
			else
				t.pnb[j] = pIdList->GetId(0);         
		}

		t.nMark = -1;

		m_Triangles[i] = t;
	}   

	pIdList->Delete();
	m_nCurrentMark = 0;
}

//------------------------------------------------------------------------
//Releases internal structures for mesh
/*virtual*/ void vtkMAFPolyDataCutOutFilterEx::ReleaseMesh()
//------------------------------------------------------------------------
{
	//release memory
	m_Vertices.clear();
	m_Triangles.clear();
}

//------------------------------------------------------------------------
//Saves the polyline in pIdsIns into the output
/*virtual*/ void vtkMAFPolyDataCutOutFilterEx::DoneCuttingPolyline(vtkPolyData* output, vtkIdList* pIdsIns)
	//------------------------------------------------------------------------
{
	//initialize output
	output->Initialize();

	//find the number of vertices
	int nVertices = (int)pIdsIns->GetNumberOfIds();
	vtkIdType* pIdPtr = pIdsIns->GetPointer(0);

	//save vertices
	vtkPoints* points = vtkPoints::New();
	points->SetNumberOfPoints(nVertices);
	for (int i = 0; i < nVertices; i++)
	{
		VERTEX& v = m_Vertices[*pIdPtr];		
		points->SetPoint(i, v.dCoord);
		pIdPtr++;
	}

	vtkCellArray* cells = vtkCellArray::New();	
	if (nVertices >= 2)
	{
		vtkIdType cell[2] = {nVertices - 1, 0};
		cells->InsertNextCell(2, cell);		

		for (int i = 0; i < nVertices; i++)
		{
			cell[0] = cell[1];
			cell[1] = i;

			cells->InsertNextCell(2, cell);		
		}
	}


	output->SetPoints(points);
	output->SetLines(cells);
	points->Delete();
	cells->Delete();
}

//------------------------------------------------------------------------
//Saves marked parts of mesh (or the whole mesh, if bExtractMarkedPartOnly is false) into output
/*virtual*/ void vtkMAFPolyDataCutOutFilterEx::DoneMesh(vtkPolyData* output, bool bExtractMarkedPartOnly)
	//------------------------------------------------------------------------
{
	//initialize output
	output->Initialize();

	int nVertices = (int)m_Vertices.size();  
	int nTriangles = (int)m_Triangles.size();    

	vtkPoints* points = vtkPoints::New();
	vtkCellArray* cells = vtkCellArray::New();
	
	if (!bExtractMarkedPartOnly)
	{
		//extract everything - it will be a bit easier
		points->SetNumberOfPoints(nVertices);
		for (int i = 0; i < nVertices; i++)
		{
			VERTEX& v = m_Vertices[i];		
			points->SetPoint(i, v.dCoord);		
		}
		cells->Allocate(nTriangles);
		for (int i = 0; i < nTriangles; i++)
		{
			TRIANGLE& t = m_Triangles[i];

			vtkIdType tri[3];
			for (int j = 0; j < 3; j++){
				tri[j] = t.aVertex[j];
			}

			cells->InsertNextCell(3, tri);			
		}
	}
	else
	{
		//OK, this will be a bit more difficult
		//find the number of vertices	
		int* pVIds = new int[nVertices];    //vertex ids
		int nValidVertices = 0;

		for (int i = 0; i < nVertices; i++)
		{
			VERTEX& v = m_Vertices[i];
			if (v.nMark != m_nCurrentMark)
				pVIds[i] = -1;
			else    
				pVIds[i] = nValidVertices++; 
		}

		_ASSERT(nValidVertices != 0 && nValidVertices != nVertices);

		//save vertices

		points->SetNumberOfPoints(nValidVertices);
		for (int i = 0; i < nVertices; i++)
		{
			VERTEX& v = m_Vertices[i];
			if (v.nMark == m_nCurrentMark) {  //set point
				points->SetPoint(pVIds[i], v.dCoord);
			}
		}

		//find the number of triangles	
		int nValidTriangles = 0;
		for (int i = 0; i < nTriangles; i++)
		{
			TRIANGLE& t = m_Triangles[i];
			if (t.nMark != m_nCurrentMark)          
				nValidTriangles++; 
		}

		cells->Allocate(nValidTriangles);
		for (int i = 0; i < nTriangles; i++)
		{
			TRIANGLE& t = m_Triangles[i];
			if (t.nMark == m_nCurrentMark)
			{
				vtkIdType tri[3];
				for (int j = 0; j < 3; j++){
					tri[j] = pVIds[t.aVertex[j]];
				}

				cells->InsertNextCell(3, tri);
			}
		}
		delete[] pVIds;
	}

	output->SetPoints(points);
	output->SetPolys(cells);
	points->Delete();
	cells->Delete();	
}

//------------------------------------------------------------------------
//Sorts all points from the cutting polyline to form a closed contour.
//Returns false, if the polyline is invalid; otherwise true and in pIds
//are ordered indices of points.
bool vtkMAFPolyDataCutOutFilterEx::SortPolylinePoints(vtkIdList* pIds)
	//------------------------------------------------------------------------
{
	int nPoints = this->CuttingPolyline->GetNumberOfPoints();  
	vtkIdType* pIdsPtr = pIds->WritePointer(0, nPoints);
	pIds->SetNumberOfIds(nPoints);

	this->CuttingPolyline->BuildCells();  //initialize cells
	int nCells = this->CuttingPolyline->GetNumberOfCells();  

	bool* CellProcessed = new bool[nCells];
	memset(CellProcessed, 0, sizeof(bool)*nCells);

	bool bOK = true;
	vtkIdType nPts, *pPts;
	double x_prev[3], x_next[3];
	int nCurrCell, nLastCell = -1;
	int nPrevCellPoint = -1, nTotalPoints = 0;

	for (int i = 0; i < nCells && bOK == true; i++)
	{ 
		//find the next cell to be processed (it append the previous one)    
		bool bFound = false;
		for (int j = 0; j < nCells; j++)
		{
			nCurrCell = (nLastCell + 1 + j) % nCells;
			if (!CellProcessed[nCurrCell])
			{
				this->CuttingPolyline->GetCellPoints(nCurrCell, nPts, pPts);
				if (nPrevCellPoint < 0 || nPrevCellPoint == pPts[0]) {
					bFound = true; break;
				}
			}      
		}

		if (!bFound) 
		{  
			//cannot continue, no continuing cell, as a last resort, check coordinates
			//the contour might be geometrically closed, even if it is not topologically
			//closed, when some point is duplicated      
			this->CuttingPolyline->GetPoint(nPrevCellPoint, x_prev);

			for (int j = 0; j < nCells; j++)
			{
				nCurrCell = (nLastCell + 1 + j) % nCells;
				if (!CellProcessed[nCurrCell])
				{
					this->CuttingPolyline->GetCellPoints(nCurrCell, nPts, pPts);
					this->CuttingPolyline->GetPoint(pPts[0], x_next);

					if (x_prev[0] == x_next[0] && x_prev[1] == x_next[1] && x_prev[2] == x_next[2]) {
						bFound = true; break;
					}
				}      
			}

			if (!bFound)
			{ //not found, even if we checked this
				bOK = false; 
				break;
			}
		}

		//copy indices
		memcpy(&pIdsPtr[nTotalPoints], pPts, sizeof(vtkIdType)*(nPts - 1)); //-1 because the last point == the first of the next cell
		nTotalPoints += nPts - 1;

		nPrevCellPoint = pPts[nPts - 1];
		CellProcessed[nLastCell = nCurrCell] = true;
	} //end for every cell

	delete[] CellProcessed;

	if (!bOK)
	{    
		vtkErrorMacro(<< "Contour not continuous nor oriented.");
		return false;
	}

	//check, if the contour is closed
	if (pIdsPtr[0] != nPrevCellPoint)
	{
		//check coordinates as last resort
		this->CuttingPolyline->GetPoint(nPrevCellPoint, x_prev);
		this->CuttingPolyline->GetPoint(pIdsPtr[0], x_next);

		if (x_prev[0] != x_next[0] || x_prev[1] != x_next[1] || x_prev[2] != x_next[2])
		{
			vtkErrorMacro(<< "Contour not closed.");
			return false;
		}
	}

	pIds->SetNumberOfIds(nTotalPoints);
	return true;
}

//------------------------------------------------------------------------
//Inserts all points from the cutting polyline into mesh.
//The routine returns indices of inserted points. Indices are ordered to form 
//oriented polyline. If something goes wrong (e.g., the input polyline is 
//incorrect, the routine returns false; true otherwise. The algorithm is as
//follows. If the point to be inserted corresponds to already existing mesh
//point, it is not inserted again but the index of mesh point is stored,
//otherwise the point is added to the list of vertices. If the point lies
//inside a triangle, this triangle is subdivided into three new triangles
//(one is replaced, two newly constructed). If the point lies on an edge,
//both triangles sharing this edge are subdivided into two new triangles
//(one is replaced, one newly constructed).
bool vtkMAFPolyDataCutOutFilterEx::InsertCutPoints(vtkIdList* pIds)
	//------------------------------------------------------------------------
{ 
	vtkSmartPointer< vtkIdList > pContourIds = vtkIdList::New();
	pContourIds->UnRegister(this);  //decrease its reference by one => refcount = 1

	if (!SortPolylinePoints(pContourIds.GetPointer()))   
		return false; //polyline is wrong, cannot be ordered  

	int nPoints = pContourIds->GetNumberOfIds();
	vtkIdType* pContourIdsPtr = pContourIds->GetPointer(0);

	pIds->SetNumberOfIds(nPoints);        //allocate pIds
	vtkIdType* pIdsPtr = pIds->WritePointer(0, nPoints);


	VERTEX vertex;
	vertex.nMark = m_nCurrentMark;  
#ifdef DIJKSTRA_CUT
	vertex.Dijstra.nDMark = 0;
#endif
	int nVertices = (int)m_Vertices.size();

	//process every point from every cell  
	int nStatus, nTriId = -1;
	int nPointsProcessed = 0;
	
	//vtkSmartPointer < vtkMAFVisualDebugger > vd = vtkMAFVisualDebugger::New();
	//vd->UnRegister(this);	//to remove our reference


	//process all points from this cell
	for (int j = 0; j < nPoints; j++)
	{
		if (this->Debug != 0){
			//Debug_Visualize_Progress2(vd.GetPointer(), pContourIds.GetPointer(), j, pIdsPtr, nPointsProcessed);
		}

		//locate the point starting from nTriId
		double* pcoords = this->CuttingPolyline->GetPoint(pContourIdsPtr[j]); 
		nTriId = LocatePoint(pcoords, nTriId, nStatus, vertex.dCoord);
		if (nStatus == ptlcOutside)
		{
			vtkErrorMacro(<< "Polyline does not lie on the surface (LocatePoint)."); 
			return false;
		}

		TRIANGLE& t = m_Triangles[nTriId];
		if ((nStatus & ptlcVertex) == ptlcVertex)
		{
			//point is already existing vertex
			pIdsPtr[nPointsProcessed] = t.aVertex[nStatus & ~ptlcVertex];

			//if the previous point was detected to be the same vertex
			//skip the insertion to avoid having in pIds both points, 
			//to avoid situations like this: ...,105,106,10,10,107,...
			if (nPointsProcessed > 0 && pIdsPtr[nPointsProcessed] == pIdsPtr[nPointsProcessed - 1])
				nPointsProcessed--;         
		}
		else
		{
			//we need to add a new point
			pIdsPtr[nPointsProcessed] = nVertices++;
			m_Vertices.push_back(vertex);

			if (nStatus == ptlcInside)
				SubdivideTriangle(nTriId, pIdsPtr[nPointsProcessed]);
			else
				Subdivide2Triangles(nTriId, nStatus & ~ptlcEdge, pIdsPtr[nPointsProcessed]);
		}		

		nPointsProcessed++;   //next point was processed
	} //end for every point 

	//if landmarks are too close to each other
	if (nPointsProcessed > 1 && pIdsPtr[0] == pIdsPtr[nPointsProcessed - 1])
		nPointsProcessed--;

	//everything is OK, finally
	pIds->SetNumberOfIds(nPointsProcessed);

	if (this->Debug != 0){
			//Debug_Visualize_Progress2(vd.GetPointer(), pContourIds.GetPointer(), -1, pIdsPtr, nPointsProcessed);
	}
	
	//now make sure that every two adjacent points in pIds are connected by an edge in the mesh
	CreateCutPointsNonIntersectingContour(pIds, nTriId);		
	return true;
}

//------------------------------------------------------------------------
//Finds the triangle containing vertex with nVertexIndex index.
//The method performs bread-first search starting with the triangle identified 
//by nStartTriId and returns id of the located triangle and in iVertexPos also 
//which vertex of the triangle is the located one.
int vtkMAFPolyDataCutOutFilterEx::FindTriangle(int nStartTriId, int nVertexIndex, int& iVertexPos)
	//------------------------------------------------------------------------
{	
	//find the first triangle
	m_nCurrentMark++; //search with new mark

	std::queue< int > queue;
	queue.push(nStartTriId);

	while (!queue.empty())
	{
		int nCurTriId = queue.front();
		queue.pop();

		TRIANGLE& tCur = m_Triangles[nCurTriId];
		if (tCur.nMark == m_nCurrentMark)
			continue; //already checked in previous steps
		
		for (iVertexPos = 0; iVertexPos < 3; iVertexPos++) {
			if (tCur.aVertex[iVertexPos] == nVertexIndex) {
				return nCurTriId;
			}
		}
		
		//we need to continue
		for (int i = 0; i < 3; i++)
		{
			if (tCur.pnb[i] >= 0)
				queue.push(tCur.pnb[i]);
		}

		tCur.nMark = m_nCurrentMark;
	}

	return -1;
}

//------------------------------------------------------------------------
//Finds the first triangle in the fan around the vertex identified by triangle id and vertex pos. 
//Input parameters are those returned by the previous call of FindFirstTriangleInFan
//or FindNextTriangleInFan. When the method returns nCurTriId contains id of the 
//returned triangle and iVertexPos the position of vertex around which the 
//search was initiated (by FindFirstTriangleInFan). It supports open fans.
//N.B. the loop is completed when the method returns the same pointer as  FindFirstTriangleInFan. 
vtkMAFPolyDataCutOutFilterEx::TRIANGLE* 
	vtkMAFPolyDataCutOutFilterEx::FindNextTriangleInOpenFan(const TRIANGLE* pSimplex, int& nCurTriId, int& iVertexPos)
	//------------------------------------------------------------------------
{
	if (pSimplex->pnb[iVertexPos] < 0)
	{
		//we are at the end => go to the end of the other side
		const TRIANGLE* pSimplex2 = pSimplex;
		int iVertexPos2 = (iVertexPos + 2) % 3;
		while (pSimplex2->pnb[iVertexPos2] >= 0)
		{
			TRIANGLE* pNbSimplex2 = &m_Triangles[nCurTriId = pSimplex2->pnb[iVertexPos2]];
			iVertexPos2 = GetFarVertexPosition(pSimplex2, pNbSimplex2, iVertexPos2);
			pSimplex2 = pNbSimplex2;
		}	

		iVertexPos = (iVertexPos2 + 1) % 3;
		return const_cast<TRIANGLE*>(pSimplex2);
	}
	else
	{
		TRIANGLE* pNbSimplex = &m_Triangles[nCurTriId = pSimplex->pnb[iVertexPos]];
		int iNbVxPos = GetFarVertexPosition(pSimplex, pNbSimplex, iVertexPos);				 		
		iVertexPos = (iNbVxPos + 2) % 3;
		return pNbSimplex;
	}
}

#ifndef DIJKSTRA_CUT
//------------------------------------------------------------------------
//Finds the triangle that contains the vertex nIndex.
//The search is done in the fan of triangles around the vertex identified by nCurTriId and iVertexPos.
//Both nCurTriId and iVertexPos are update to contain the located vertex (or the last tested one) 
//The method returns NULL, if triangle has not been located.
vtkMAFPolyDataCutOutFilterEx::TRIANGLE* 
	vtkMAFPolyDataCutOutFilterEx::FindTriangleInFanV(int nIndex, int& nCurTriId, int &iVertexPos)
//------------------------------------------------------------------------
{
	TRIANGLE* pStartSimplex = FindFirstTriangleInFan(nCurTriId, iVertexPos);
	TRIANGLE* pSimplex = pStartSimplex;
	do
	{
		//check, if we have nNextIndex in the other two vertices
		int iV1 = (iVertexPos + 1) % 3;
		if (nIndex == pSimplex->aVertex[iV1]) {
			iVertexPos = iV1; 
			return pSimplex;
		}

		//the other will be tested in the next triangle

		pSimplex = FindNextTriangleInFan(pSimplex, nCurTriId, iVertexPos);
	} while (pStartSimplex != pSimplex);

	return NULL;
}

//------------------------------------------------------------------------
//Finds the mesh edge intersected (in projection) by the line going towards the vertex nIndex.
//	The method processes triangles in the fan of triangles around the vertex FanVertIndex identified by 
//	nCurTriId and iVertexPos parameters. N.B. these parameters are updated as the method proceeds.
//	When the method finds the triangle with an edge intersected by the line FanVertIndex - nIndex projected onto the plane 
//	of this triangle, it determines the status of intersection (ptlcEdge, ptlcVertex) and computes the 
//	coordinates of intersection.  The method may return NULL, if due to some numeric inaccuracy 
//	something went terribly wrong, otherwise the triangle containing x_out point
vtkMAFPolyDataCutOutFilterEx::TRIANGLE* 
	vtkMAFPolyDataCutOutFilterEx::FindTriangleInFanE(
	int nIndex, int& nCurTriId, int &iVertexPos,
	int& nStatus, double* x_out)
//------------------------------------------------------------------------
{
	const static double eps = 1e-16;			
	TRIANGLE* pStartSimplex = FindFirstTriangleInFan(nCurTriId, iVertexPos);
	
	//create line 
	double* P = m_Vertices[pStartSimplex->aVertex[iVertexPos]].dCoord;
	double* Q = m_Vertices[nIndex].dCoord;					

	double l_dir[3];
	for (int k = 0; k < 3; k++) {
		l_dir[k] = Q[k] - P[k];
	}

	double dblBestAdjT = DBL_MAX, dblBestT = 0.0;	//optimal dblBestAdjT should be 0.0
	int nBestTriId = -1, iBestVertexPos = -1;
	TRIANGLE* pBestSimplex = NULL;

	TRIANGLE* pSimplex = pStartSimplex;
	do
	{
		//construct the plane for the current triangle		
		int iV1 = (iVertexPos + 1) % 3;
		int iV2 = (iVertexPos + 2) % 3;

		double* A = m_Vertices[pSimplex->aVertex[iV1]].dCoord;
		double* B = m_Vertices[pSimplex->aVertex[iV2]].dCoord;
		double* C = m_Vertices[pSimplex->aVertex[iVertexPos]].dCoord;

		//compute normal of this plane
		double u[3], v[3], n[3];
		for (int j = 0; j < 3; j++)
		{
			u[j] = B[j] - A[j];
			v[j] = C[j] - A[j];      
		}

		vtkMath::Cross(u, v, n);
		if (vtkMath::Normalize(n) <= eps) 				
			goto nextLoop;	//ignore this triangle, it is degenerated into a line

		//we have a valid plane, now project l_dir onto this plane
		//according to P.Schneider: Geometric Tools for Computer Graphics, pg. 665
		//projection of vector l_dir onto plane A*n + d = 0 can be computed as
		//l_dir_proj = l_dir  - (l_dir *n)*n assuming that vector n is normalized
		double l_dir_proj[3];
		double dblVal = vtkMath::Dot(l_dir, n);
		for (int j = 0; j < 3; j++) {
			l_dir_proj[j] = (l_dir[j] - dblVal*n[j]);
		}

		//does P + l_dir_proj*t intersect line AB?
		dblVal = vtkMath::Dot(u, l_dir_proj);
		if (fabs(dblVal - 1.0) <= eps)
			goto nextLoop;	//apparently u is parallel with l_dir_proj

		//there is some intersection => compute it
		//according to http://mathworld.wolfram.com/Line-LineIntersection.html,
		//intersection of lines A + t*u and P + s*dir => t = ((P-A) x dir)*(u x dir))/ |(u x dir)|^2
		for (int j = 0; j < 3; j++) {
			v[j] = P[j] - A[j];
		}

		double c1[3], c2[3];
		vtkMath::Cross(u, l_dir_proj, c1);
		vtkMath::Cross(v, l_dir_proj, c2);

		//OK, line AB is intersected by line l_dir_proj at t
		double t = vtkMath::Dot(c1, c2) / 
			(c1[0]*c1[0] + c1[1]*c1[1] + c1[2]*c1[2]);
		
		double tAdj = fabs(t - 0.5);
		if (tAdj < dblBestAdjT)
		{
			//this intersection is better than the other
			//however, we need to check, if  it is also intersected by the half-line
			//calculate also s, since s must be > 0
			int iMaxSlope = -1;
			double dblMaxSlope = 0.0;
			for (int j = 0; j < 3; j++) {
				double absslope = fabs(l_dir_proj[j]);
				if (absslope > dblMaxSlope) {
					iMaxSlope = j; dblMaxSlope = absslope;
				}
			}

			double s = (t*u[iMaxSlope] - v[iMaxSlope]) / l_dir_proj[iMaxSlope];
			if (s <= 0.0) {
					goto nextLoop;	//the intersection is not in the direction toward the point Q
			}

			dblBestAdjT = tAdj; dblBestT = t; 
			nBestTriId = nCurTriId;
			iBestVertexPos = iVertexPos;
			pBestSimplex = pSimplex;
		}

nextLoop:
		pSimplex = FindNextTriangleInFan(pSimplex, nCurTriId, iVertexPos);
	} while(pStartSimplex != pSimplex);

	//we are ready
	if (pBestSimplex == NULL)
		return NULL;	//FATAL ERROR, we could not find it, some restart is needed

	//determine what is our next point
	nCurTriId = nBestTriId;
	iVertexPos = iBestVertexPos;

	int iV1 = (iVertexPos + 1) % 3;
	int iV2 = (iVertexPos + 2) % 3;

	double* A = m_Vertices[pBestSimplex->aVertex[iV1]].dCoord;
	double* B = m_Vertices[pBestSimplex->aVertex[iV2]].dCoord;
				
	double u[3];
	for (int j = 0; j < 3; j++){
		u[j] = B[j] - A[j];			      
	}

	double tolerance = this->PointTolerance  / vtkMath::Norm(u);

	//we are on the valid edge, great
	if (dblBestT <= tolerance)
	{
		//the next point is A
		for (int j = 0; j < 3; j++) {
			x_out[j] = A[j];
		}

		nStatus = ptlcVertex | (iVertexPos = iV1);					
	} 
	else if (dblBestT >= 1 - tolerance)
	{
		//the next point will be B
		for (int j = 0; j < 3; j++) {
			x_out[j] = B[j];
		}

		nStatus = ptlcVertex | (iVertexPos = iV2);
	}
	else
	{
		//compute x_out
		for (int j = 0; j < 3; j++) {
			x_out[j] = A[j] + dblBestT*u[j];
		}

		nStatus = ptlcEdge | iV1;	//iVertexPos must remain unchanged!!!
	}

	return pBestSimplex;
}


//------------------------------------------------------------------------
//Finds the mesh vertex closest to the ray going towards the vertex nIndex.
//	The method processes triangles in the fan of triangles around the vertex FanVertIndex identified by 
//	nCurTriId and iVertexPos parameters. N.B. these parameters are updated as the method proceeds.
//When the method exits, nCurTriId and iVertexPos contain the located vertex
//The method never returns NULL
vtkMAFPolyDataCutOutFilterEx::TRIANGLE* 
	vtkMAFPolyDataCutOutFilterEx::FindTriangleInFanL(
	int nIndex, int& nCurTriId, int &iVertexPos)
//------------------------------------------------------------------------
{
	TRIANGLE* pStartSimplex = FindFirstTriangleInFan(nCurTriId, iVertexPos);

	//create ray
	double* P = m_Vertices[pStartSimplex->aVertex[iVertexPos]].dCoord;
	double* Q = m_Vertices[nIndex].dCoord;					

	double l_dir[3];
	for (int k = 0; k < 3; k++) {
		l_dir[k] = Q[k] - P[k];
	}

	vtkMath::Normalize(l_dir);	//make it unit
	
	double dblMinDist = DBL_MAX;
	int iBestTriId = -1, iBestVertexPos = -1;

	TRIANGLE* pSimplex = pStartSimplex;
	do
	{
		//check, if we have nNextIndex in the other two vertices
		int iV1 = (iVertexPos + 1) % 3;	
		//if (m_Vertices[pSimplex->aVertex[iV1]].nMark != m_nCurrentMark)	//avoid infinite loops
		{
			double* x = m_Vertices[pSimplex->aVertex[iV1]].dCoord;

			//the closest point on the ray can be computed as:
			//P + (u*(X - P))*u, where u is normalized vector l_dir
			double v[3];
			for (int k = 0; k < 3; k++){			
				v[k] = x[k] - P[k];
			}

			double w = vtkMath::Dot(l_dir, v);		
			if (w >= 0)
			{		
				//closest point lies on the line supporting the edge
				double x_cl[3];
				for (int k = 0; k < 3; k++) {
					x_cl[k] = P[k] + w*l_dir[k];
				}


				double dblDist = vtkMath::Distance2BetweenPoints(x_cl, x);	
				if (dblDist < dblMinDist) {
					iBestTriId = nCurTriId;
					iBestVertexPos = iV1;
					dblMinDist = dblDist;
				}
			}
		}		

		pSimplex = FindNextTriangleInFan(pSimplex, nCurTriId, iVertexPos);
	} while (pStartSimplex != pSimplex);
	
	_VERIFY_RETVAL(iBestTriId >= 0, NULL)		
	return FindFirstTriangleInFan(nCurTriId = iBestTriId, iVertexPos = iBestVertexPos);
}


//------------------------------------------------------------------------
/** Finds all triangles sharing the vertex FanVertIndex identified by nCurTriId and iVertexPos parameters.
pOutput is allocated in the method and the caller is responsible for its release when the data is no longer needed.
The method returns the number of retrieved triangles.
N.B. triangles are ordered. */
int vtkMAFPolyDataCutOutFilterEx::GetTriangleFan(int nCurTriId, int iVertexPos, TRIANGLE**& pOutput)
	//------------------------------------------------------------------------
{
	int nTriangles = 0;
	TRIANGLE* pStartSimplex = FindFirstTriangleInFan(nCurTriId, iVertexPos);
	TRIANGLE* pSimplex = pStartSimplex;
	do 
	{
		nTriangles++;
		pSimplex = FindNextTriangleInFan(pSimplex, nCurTriId, iVertexPos);
	} 
	while (pStartSimplex != pSimplex);
	
	pOutput = new TRIANGLE*[nTriangles];
	for (int i = 0; i < nTriangles; i++) 
	{		
		pOutput[0] = pSimplex;
		pSimplex = FindNextTriangleInFan(pSimplex, nCurTriId, iVertexPos);
	} 
		
	return nTriangles;
}
#else

//------------------------------------------------------------------------
// Moves vertex with the given nIndex in the direction to head of Dijkstra priority queue until the proper place is found.
//Returns the new index of vertex at head of the queue. 
void vtkMAFPolyDataCutOutFilterEx::DijkstraPQDown(int nIndex, int* pPQ, int nPQSize)
	//------------------------------------------------------------------------
{
	double dblCurPri = m_Vertices[nIndex].Dijstra.dSum;
	int nPQIndex = m_Vertices[nIndex].Dijstra.iPQIndex;

	while (2*nPQIndex + 1 < nPQSize)
	{	
		int iChild = 0;
		int nPQChildIndex[2];		
		double dblChildPri[2];
			
		nPQChildIndex[0] = 2*nPQIndex + 1;
		dblChildPri[0] = m_Vertices[pPQ[nPQChildIndex[0]]].Dijstra.dSum;
		
		nPQChildIndex[1] = nPQChildIndex[0] + 1;				
		if (nPQChildIndex[1] < nPQSize)
		{
			dblChildPri[1] = m_Vertices[pPQ[nPQChildIndex[1]]].Dijstra.dSum;
			if (dblChildPri[1] < dblChildPri[0]) {
				iChild++;
			}
		}

		if (dblCurPri <= dblChildPri[iChild])
			break;	//we are ready

		//swap items
		DijkstraPQSwap(nPQIndex, nPQChildIndex[iChild], pPQ);		
		nPQIndex = nPQChildIndex[iChild];
	}
}

//------------------------------------------------------------------------
//Moves vertex with the given nIndex in the direction to tail of Dijkstra priority queue until the proper place is found.
//Returns the new index of vertex at head of the queue.
void vtkMAFPolyDataCutOutFilterEx::DijkstraPQUp(int nIndex, int* pPQ)
	//------------------------------------------------------------------------
{
	double dblCurPri = m_Vertices[nIndex].Dijstra.dSum;
	int nIndexPQ = m_Vertices[nIndex].Dijstra.iPQIndex;

	while (nIndexPQ > 0)
	{
		int nParIndexPQ = (nIndexPQ + 1) / 2 - 1;
		double dblParPri = m_Vertices[pPQ[nParIndexPQ]].Dijstra.dSum;
		
		if (dblParPri <= dblCurPri)
			break;	//we are ready

		//swap nIndexPQ and nParIndexPQ
		DijkstraPQSwap(nIndexPQ, nParIndexPQ, pPQ);

		nIndexPQ = nParIndexPQ;
	}	
}

//------------------------------------------------------------------------
//Deletes the minimum from the PQ. 
//Returns the new index of vertex at head of the queue.
int vtkMAFPolyDataCutOutFilterEx::DijkstraPQDeleteMin(int* pPQ, int& nPQSize)
//------------------------------------------------------------------------
{
	int nRetIndex = pPQ[0];
	if (0 != --nPQSize)
	{
		//swap the first and last item of PQ
		DijkstraPQSwap(0, nPQSize, pPQ);	
		DijkstraPQDown(pPQ[0], pPQ, nPQSize);
	}
	return nRetIndex;
}

//------------------------------------------------------------------------
//Inserts the new vertex into the PQ.
//N.B. the value of the vertex must be already specified.
//Returns the new index of vertex at head of the queue. 
void vtkMAFPolyDataCutOutFilterEx::DijkstraPQInsert(int nIndex, int* pPQ, int& nPQSize)
//------------------------------------------------------------------------
{
	pPQ[nPQSize] = nIndex;
	m_Vertices[nIndex].Dijstra.iPQIndex = nPQSize;	
	nPQSize++;

	DijkstraPQUp(nIndex, pPQ);	
}

//------------------------------------------------------------------------
//Finds the shortest path from vertex nFromIndex to nToIndex.
//All indices on the path except for the last (nToIndex) are inserted into pContourIds and vertices are marked.
//nTriId is a helper ID to a valid triangle close to the first point to speed-up the process. 
//Returns ID to a valid triangle close to nToIndex to speed-up the further process.
int vtkMAFPolyDataCutOutFilterEx::CreateDijskraCutPath(int nFromIndex, int nToIndex, vtkIdList* pContourIds, int nTriId, int* pPQBuf)
//------------------------------------------------------------------------
{
	// Initializations
	int* pPQ = pPQBuf;
	if (pPQ == NULL){
		pPQ = new int[m_Vertices.size()];
	}

	m_nCurrentDijkstraMark++;	//next dijkstra mark
	
	VERTEX& vini = m_Vertices[nFromIndex];
	if (vini.Dijstra.nDMark == 0)
	{
		//find the first triangle containing the first point (it has not been found in previous iterations)
		vini.Dijstra.iVertexPos = 0;
		vini.Dijstra.nTriId = FindTriangle(nTriId, nFromIndex, vini.Dijstra.iVertexPos);
	}	
	
	vini.Dijstra.dSum = 0.0;
	//vini.Dijstra.iNextPQ = vini.Dijstra.iPrevPQ = -1;
	vini.Dijstra.nDMark = m_nCurrentDijkstraMark;	//already in PQ
	vini.Dijstra.nPrevVertexId = -1;
	
	int nPQSize = 0;
	DijkstraPQInsert(nFromIndex, pPQ, nPQSize);
	
	int nToIndexDegree = 0;	
	while (nPQSize > 0)
	{
		int uIndex = DijkstraPQDeleteMin(pPQ, nPQSize);
		if (uIndex == nToIndex){	
				break;	//we have found it => we are finished			
		}

		//get the vertex on top
		VERTEX& u = m_Vertices[uIndex];		
		u.Dijstra.nDMark = -m_nCurrentDijkstraMark;	//mark the vertex to be already processed

		///search for each neighbor of the current vertex
		int nCurTriId = u.Dijstra.nTriId;
		int iVertexPos = u.Dijstra.iVertexPos;
		
		TRIANGLE* pStartSimplex = FindFirstTriangleInFan(nCurTriId, iVertexPos);
		TRIANGLE* pSimplex = pStartSimplex;
		do 
		{
			int vIndex = pSimplex->aVertex[(iVertexPos + 1) % 3];	//get next vertex
			VERTEX& v = m_Vertices[vIndex];

			//if the vertex v has not been yet processed
			if (v.Dijstra.nDMark != -m_nCurrentDijkstraMark)
			{
				//compute distance between u and v
				double dist_uv = 0.0;
				for (int i = 0; i < 3; i++) 
				{
					double d = v.dCoord[i] - u.dCoord[i];
					dist_uv += d*d;
				}
				
				double alt = u.Dijstra.dSum + sqrt(dist_uv);
				if (v.Dijstra.nDMark == m_nCurrentDijkstraMark)
				{
					//v is already in PQ, so check, if we have found shorter path					
					if (alt < v.Dijstra.dSum) 
					{
						//we have a shorter path => alter it
						v.Dijstra.dSum = alt;
						v.Dijstra.nPrevVertexId = uIndex;

						DijkstraPQUp(vIndex, pPQ);
					}
				}
				else
				{
					//the vertex is not in PQ, so let us add it there
					v.Dijstra.dSum = alt;
					v.Dijstra.nTriId = nCurTriId;
					v.Dijstra.iVertexPos = (iVertexPos + 1) % 3;
					v.Dijstra.nPrevVertexId = uIndex;
					v.Dijstra.nDMark = m_nCurrentDijkstraMark;

					DijkstraPQInsert(vIndex, pPQ, nPQSize);
				}			
			}

			//get the next vertex
			pSimplex = FindNextTriangleInFan(pSimplex, nCurTriId, iVertexPos);
		} 
		while (pStartSimplex != pSimplex);
	}

	//we are ready with Dijkstra, so store now the path
	int nPathLen = 0, nCurIndex = nToIndex;		
	while (m_Vertices[nCurIndex].Dijstra.nPrevVertexId >= 0) {
		nPathLen++; nCurIndex = m_Vertices[nCurIndex].Dijstra.nPrevVertexId;
		pContourIds->InsertNextId(-1);	//reserve value for ID
	}

	int nIdPos = pContourIds->GetNumberOfIds();	
	nCurIndex = nToIndex;	

	while (m_Vertices[nCurIndex].Dijstra.nPrevVertexId >= 0) 
	{
		nCurIndex = m_Vertices[nCurIndex].Dijstra.nPrevVertexId;
		pContourIds->SetId(--nIdPos, nCurIndex);
	}


	if (pPQBuf == NULL){
		delete[] pPQ;
	}
	return m_Vertices[nToIndex].Dijstra.nTriId;
}
#endif

//------------------------------------------------------------------------
//Updates the contour in pIds list in such a way that points form a connected chain (connected by
//mesh edges) and there are no duplicities in the contour (this is fixed). 
//nTriId is a helper ID to a valid triangle close to the first point to speed-up the process.
void vtkMAFPolyDataCutOutFilterEx::CreateCutPointsNonIntersectingContour(vtkIdList* pIds, int nTriId)
	//------------------------------------------------------------------------
{	
	int nLastNumberOfIds = pIds->GetNumberOfIds(); 
#ifndef DIJKSTRA_CUT
	while (true)
	{		
#endif
		CreateCutPointsContour(pIds, nTriId);
#ifndef DIJKSTRA_CUT
		int nNewNumberOfIds = pIds->GetNumberOfIds(); 
		if (nNewNumberOfIds == nLastNumberOfIds)
			break;	//we have not changing list

		nLastNumberOfIds = nNewNumberOfIds;
	} 
#endif

	//vtkSmartPointer < vtkMAFVisualDebugger > vd = vtkMAFVisualDebugger::New();
	//vd->UnRegister(this);	//to remove our reference
	
	//BES: 13.3.2012 - The input data is really nasty!!!
	//It is necessary to do the following: if there is a crossing, split it there into two contours 
	//and reorient one of them so that the path does not cross but only touches the other,
	//after that vertices where the contour touches itself must be split and the contour corrected

	vtkIdList* pRetList = vtkIdList::New();

	//check for duplicities
	int nPoints = pIds->GetNumberOfIds();
	vtkIdType* pIdsPtr = pIds->GetPointer(0);		
	pRetList->InsertNextId(pIdsPtr[0]);

	typedef std::map< int, int > VertMap;	//Key = vertex index, Value = position in pIds
	VertMap dupl;	//to detect duplicated vertices
	dupl[pIdsPtr[0]] = 0;	

	for (int i = 1; i < nPoints; i++) 
	{
		VertMap::const_iterator it = dupl.find(pIdsPtr[i]);
		if (it == dupl.end())
		{
			pRetList->InsertNextId(pIdsPtr[i]);
			dupl[pIdsPtr[i]] = i;
		}
		else
		{
			if (this->Debug != 0) {
				//Debug_Visualize_Progress3(vd.GetPointer(), pIds, i, pRetList);
			}

			//we have here the duplicated part of the contour
			//get all adjacent indices, so we could fix it
			int iPos = it->second;	//iPos should be always > 0
			vtkIdType ptInOut[4] = { pIdsPtr[iPos - 1], pIdsPtr[iPos + 1], 
															pRetList->GetId(pRetList->GetNumberOfIds() - 1), //pIdsPtr[i  - 1] may not contain the correct ID
															pIdsPtr[(i + 1) % nPoints] };

			//Now, search the fan around the current (common) vertex pIdsPtr[i] starting from outPt[0]
			//to outPt[1] and check, if any vertex from inPt is present in the chain
			//If none or both are present, then contour touches itself at pIdsPtr[i],
			//otherwise, the contour intersects itself at pIdsPtr[i] (this must be corrected)
#pragma region REORIENTATION
			int iVertexPos = 0;
			int nCurTriId = FindTriangle(nTriId, pIdsPtr[i], iVertexPos);

			//RELEASE NOTE: this assumes closed fans!
			int ptOrders[4] = {-1, -1, -1, -1};
			int nPtFound = 0;

			int nTriangles = 0;
			TRIANGLE* pStartSimplex = FindFirstTriangleInFan(nCurTriId, iVertexPos);
			TRIANGLE* pSimplex = pStartSimplex;
			do 
			{				
				int iVertex = (pSimplex->aVertex[(iVertexPos + 1) % 3]);
				for (int k = 0; k < 4; k++) {
					if (iVertex == ptInOut[k]) {
						ptOrders[k] = nPtFound++; 
					}
				}				

				nTriangles++;
				pSimplex = FindNextTriangleInFan(pSimplex, nCurTriId, iVertexPos);
			} 
			while (pStartSimplex != pSimplex);

			//reorder ptOrders, so the InOut chain starts at 0
			int nAdd = 4 - ptOrders[0];
			for (int k = 0; k < 4; k++) {
				ptOrders[k] = (ptOrders[k] + nAdd) % 4;
			}

			if (ptOrders[1] == 2)
			{
				//So, we the contour intersects itself => we need to take the contour from ptInOut[1] to ptInOut[2]
				//and invert it, so the contour will touch itself but not intersect
				//there are at least two vertices in the considered contour [because otherwise ptOrders[1] would be 1]
				InvertCurve(pRetList->GetPointer(0), iPos + 1, pRetList->GetNumberOfIds() - 1);

				int tmp = ptOrders[1];
				ptOrders[1] = ptOrders[2];
				ptOrders[2] = tmp;

				//update ptInOut
				tmp = ptInOut[1];
				ptInOut[1] = ptInOut[2];
				ptInOut[2] = tmp;
			}
#pragma endregion

#pragma region VERTEX SPLIT
			//At this moment, we have contour that touches itself at vertex pIdsPtr[i]
			//we need to process triangles between OUT vertices in triangle ORDER =>
			//swap our OUT vertices, if their order go against the triangle order
			int iSwpRetList = pRetList->GetNumberOfIds();
			if (ptOrders[2] > ptOrders[3])
			{
				vtkIdType tmp = ptInOut[2];
				ptInOut[2] = ptInOut[3];
				ptInOut[3] = tmp;
			}

			//get the first triangle, where everything starts, i.e., which contains ptInOut[2]
			do 
			{
				int iVertex = (pSimplex->aVertex[(iVertexPos + 1) % 3]);
				if (iVertex == ptInOut[2])
					break;

				pSimplex = FindNextTriangleInFan(pSimplex, nCurTriId, iVertexPos);
			} 
			while (pStartSimplex != pSimplex);

			//pSimplex now contains reference to the triangle whose iVertexPos + 1 == OUT but it
			//the first triangle to be split is its neighbor
			pStartSimplex =  pSimplex = FindNextTriangleInFan(pSimplex, nCurTriId, iVertexPos);
			int nCurTriId2 = nCurTriId, iVertexPos2 = iVertexPos;	//store info

			//now cycle until the other OUT point is found			
			int nValidTri = -1;						
			for (int k = 0; k < nTriangles; k++)
			{
				int iVertex = (pSimplex->aVertex[(iVertexPos + 1) % 3]);
				if (iVertex == ptInOut[3]) {
					nValidTri = k;	break; //nValidTri denotes number of triangles between both OUT vertices
				}

				_ASSERTE(iVertex != ptInOut[0] && iVertex != ptInOut[1]);
				pSimplex = FindNextTriangleInFan(pSimplex, nCurTriId, iVertexPos);
			} 			

			if (nValidTri == 0)
			{
				//singular case: both OUT vertices are in the same triangle, so all that we need
				//to do, is to skip pIdsPtr[i] and continue, i.e., this branch will be empty
			}
			else
			{
				//for each triangle, add  vertex on its more distant edge
				double* A = m_Vertices[pIdsPtr[i]].dCoord;

				VERTEX vertex;
				vertex.nMark = 0;
#ifdef DIJKSTRA_CUT
				vertex.Dijstra.nDMark = 0;
#endif
				nCurTriId = nCurTriId2; iVertexPos = iVertexPos2;
				for (int iCurTri = 0; iCurTri < nValidTri; iCurTri++)
				{
					//edge (iVerexPos + 1) % 3 will be split
					int iEdgePos = iVertexPos; //(iVertexPos + 2) % 3;
					double* B = m_Vertices[m_Triangles[nCurTriId].aVertex[(iEdgePos + 1) % 3]].dCoord;					
					for (int k = 0; k < 3; k++) {
						vertex.dCoord[k] = B[k] - A[k];
					}

					for (int k = 0; k < 3; k++) {
						vertex.dCoord[k] = A[k] + vertex.dCoord[k]*
#if defined(_DEBUG)
						(this->Debug == 0 ? this->PointTolerance : 0.5);	//0.5, so it is visible
#else
						this->PointTolerance;
#endif
					}

					//insert the new vertex
					int nCurIndex = (int)m_Vertices.size();
					m_Vertices.push_back(vertex);					

					Subdivide2Triangles(nCurTriId, iEdgePos, nCurIndex);
					pRetList->InsertNextId(nCurIndex);

					//after subdivision, our next triangle has changed but we know exactly where it is
					nCurTriId = m_Triangles[nCurTriId].pnb[1];
					iVertexPos = 0;	//the new point was put at the pos 
				}

				if (ptOrders[2] > ptOrders[3]) {
					//we need to invert the chain of just created vertices
					InvertCurve(pRetList->GetPointer(0), iSwpRetList, pRetList->GetNumberOfIds() - 1);
				}
			} //end [there are some triangles to be subdivided]

			if (this->Debug != 0) {
				//Debug_Visualize_Progress3(vd.GetPointer(), pIds, i, pRetList);
			}
#pragma endregion
		}
	}

	if (this->Debug != 0) {
		//Debug_Visualize_Progress3(vd.GetPointer(), pIds, 0, pRetList);
	}

	//copy contours to pIds	
	pIds->DeepCopy(pRetList);
	pRetList->Delete();	
}

//------------------------------------------------------------------------
//Inverts the curve stored in pInvertPtr at indices from iLeft to iRight (inclusively)
void vtkMAFPolyDataCutOutFilterEx::InvertCurve(vtkIdType* pInvertPtr, int iLeft, int iRight)
{

	while (iLeft < iRight)
	{
		vtkIdType tmp = pInvertPtr[iLeft];
		pInvertPtr[iLeft] = pInvertPtr[iRight];
		pInvertPtr[iRight] = tmp;

		iLeft++; iRight--;
	}
}

//------------------------------------------------------------------------
//Updates the contour in pIds list in such a way that points forms a connected chain (connected by
//	mesh edges)  connected but does not guarantee that the contour will be manifold.
//	nTriId is a helper ID to a valid triangle close to the first point to speed-up the process. 
//	The method is called by  CreateCutPointsNonIntersectingContour.
void vtkMAFPolyDataCutOutFilterEx::CreateCutPointsContour(vtkIdList* pIds, int nTriId)
	//------------------------------------------------------------------------
{	
	//vtkSmartPointer < vtkMAFVisualDebugger > vd = vtkMAFVisualDebugger::New();
	//vd->UnRegister(this);	//to remove our reference


	VERTEX vertex;
	vertex.nMark = 0;
#ifdef DIJKSTRA_CUT
	vertex.Dijstra.nDMark = 0;
#endif

	int nVertices = (int)m_Vertices.size();

	vtkSmartPointer< vtkIdList > pContourIds = vtkIdList::New();
	pContourIds->UnRegister(this);  //decrease its reference by one => refcount = 1

	int nPoints = pIds->GetNumberOfIds();
	vtkIdType* pIdsPtr = pIds->GetPointer(0);	

#ifdef DIJKSTRA_CUT
	int* pPQ = new int[nVertices];

	int nCurTriId = nTriId;
#else		

	//find the first triangle containing the first point
	int iVertexPos = 0;
	int nCurTriId = FindTriangle(nTriId, pIdsPtr[0], iVertexPos);
#endif

	m_nCurrentMark++;
	for (int i = 0; i < nPoints; i++)
	{
		//find the path from the current index to the next index
		int nCurIndex = pIdsPtr[i];
		int nNextIndex = pIdsPtr[(i + 1) % nPoints];		

#ifdef DIJKSTRA_CUT
		if (this->Debug != 0) {
				//Debug_Visualize_Progress(vd.GetPointer(), pIds, pContourIds.GetPointer(), nCurIndex, nNextIndex);
		}

		nCurTriId = CreateDijskraCutPath(nCurIndex, nNextIndex, pContourIds.GetPointer(), nCurTriId, pPQ);
#else
		while (true)
		{			
			//record nCurIndex
			pContourIds->InsertNextId(nCurIndex);
			m_Vertices[nCurIndex].nMark = m_nCurrentMark;

			if (this->Debug != 0) {
				Debug_Visualize_Progress(vd.GetPointer(), pIds, pContourIds.GetPointer(), nCurIndex, nNextIndex);
			}

			//search the triangle fan around nCurIndex (P) to find the triangle containing nNextIndex				
			if (NULL != FindTriangleInFanV(nNextIndex, nCurTriId, iVertexPos))
				break; //we have found it, so we may proceed with the next point

			//OK, so we have to find triangle that is intersected by the line from nCurIndex -> nNextIndex
			int nStatus;						
			TRIANGLE* pSimplex = FindTriangleInFanE(
#ifndef CUTOUT_IN_ONE_DIRECTION
				nNextIndex, 
#else
				direction,
#endif
				nCurTriId, iVertexPos, nStatus, vertex.dCoord);
			if (pSimplex != NULL)
			{									
				//point is already existing vertex => just continue with this vertex
				if ((nStatus & ptlcVertex) == ptlcVertex)
					nCurIndex = pSimplex->aVertex[iVertexPos];		
				else
				{
					//we need to add a new point								
					m_Vertices.push_back(vertex);
					nCurIndex = nVertices++;

					Subdivide2Triangles(nCurTriId, nStatus & ~ptlcEdge, nCurIndex);
					iVertexPos = 2;	//the new point was put at the pos 2
				}
			}
			else
			{
				//oops, due to numerical problems we have not found the triangle where to continue
				//we will continue with the vertex of triangle that is closest to the current ray 
				pSimplex = FindTriangleInFanL(
#ifndef CUTOUT_IN_ONE_DIRECTION
				nNextIndex, 
#else
				direction,
#endif
				nCurTriId, iVertexPos);	
				_VERIFY_CMD(pSimplex != NULL, goto finish);

				nCurIndex = pSimplex->aVertex[iVertexPos];
			}
		} //end while (true)
#endif
	} //end for
#ifndef DIJKSTRA_CUT
finish:
#else
	delete[] pPQ;
#endif
	//copy contours to pIds	
	pIds->DeepCopy(pContourIds.GetPointer());
}

//------------------------------------------------------------------------
//Locates the given point in the current mesh.
//It searches for the triangle that contains the given point x in the mesh.
//The indentifier of this triangle is returned and the mutual position of 
//the point and the triangle is returned in nStatus, the coordinates
//to be stored for the point (might be shifted in comparison to x) are
//returned in x_out variable. When the method is called for the first time, 
//nTriId must be -1 (nStatus is ignored) and the algorithm performs 
//a brute-force search (using VTK to speed up the process). The successive
//calls have nTriId set to the identifier returned from the previous call
//and nStatus unchanged. Only the triangles around the last processed
//points are checked.	The method supposes that the caller modifies
//the triangulation regarding to the returned values.
//N.B. This method never returns ERROR.
int vtkMAFPolyDataCutOutFilterEx::LocatePoint(double* x, int nTriId, 
	int& nStatus, double* x_out)
	//------------------------------------------------------------------------
{	
	double tol = 1e-16;

	if (nTriId < 0)
	{    
		//locate the initial cell (we use VTK, as it has some internal structures
		//to avoid brute-force search)
		int subId;    
		double x_par[3], weights[3];		
		while (tol <= 1e-8 &&
			(nTriId = vtkPolyData::SafeDownCast(GetInput())->FindCell(x, NULL, NULL, -1, tol, subId, x_par, weights)) < 0)
		{
			tol *= 10;	//increase the tolerance
		}

		if (nTriId < 0)
			nTriId = 0;	//we will need to perform brute search
	}

	m_nCurrentMark++; //search with new mark

	std::queue< int > queue;
	queue.push(nTriId);

	//here is stored the closest triangle, even, if it is beyond the tolerance	
	double dblMinDist = DBL_MAX;
	int nMinDistTriId = -1;

	int iStage = 0;
	while (true)	//perform at least one stage
	{
		while (!queue.empty())
		{
			int nCurTriId = queue.front();		
			queue.pop();

			TRIANGLE& tCur = m_Triangles[nCurTriId];
			if (tCur.nMark == m_nCurrentMark)
				continue; //already checked in previous steps  		

			double* A = m_Vertices[tCur.aVertex[0]].dCoord;
			double* B = m_Vertices[tCur.aVertex[1]].dCoord;
			double* C = m_Vertices[tCur.aVertex[2]].dCoord;

			//compute normal of this plane
			double u[3], v[3], n[3];
			for (int j = 0; j < 3; j++)
			{
				u[j] = B[j] - A[j];
				v[j] = C[j] - A[j];      
			}

			vtkMath::Cross(u, v, n);
			if (vtkMath::Normalize(n) >= tol)     
			{
				//we have a valid plane, now compute the closest point
				//according to P.Schneider: Geometric Tools for Computer Graphics, pg. 665
				//projection of vector xA onto plane A*n + d = 0 can be computed as
				//xA_proj = xA - (xA*n)*n assuming that vector n is normalized, thus
				//the closest point x_cl = A + xA_proj
				double x_cl[3], xA[3];
				for (int j = 0; j < 3; j++) {          
					xA[j] = x[j] - A[j];
				}

				double dblVal = vtkMath::Dot(xA, n);
				for (int j = 0; j < 3; j++) {
					x_cl[j] = A[j] + (xA[j] - dblVal*n[j]);
				}

				//in stage 2 only one triangle is selected and the best edge or vertex is chosen
				if (iStage != 2)
				{
					//in the early stage 0, we filter out planes too far from x
					//but we are trying to find out the triangle for whose x_cl is the closest
					//the stage 1 is performed, if stage 0 failed completely, i.e., we have no idea about
					//the closest triangle, so we need to repeat the whole process without filtering planes too far from x
					double dblDist2 = iStage == 0 ? vtkMath::Distance2BetweenPoints(x, x_cl) : 0.0;
					if (dblDist2 <= this->PointTolerance)
					{
						//point lies on the plane, use an imprecise (using tolerances) test
						//to decide whether the point lies on edges (or in vertices)
						double dblDistABC;
						GetPointInTriangleStatus(x_cl, A, B, C, nStatus, x_out, &dblDistABC);
						if (nStatus != ptlcOutside)  //edges or vertex        
							return nCurTriId;

						//we are outside the triangle
						if (dblDistABC < dblMinDist) {
							dblMinDist = dblDistABC; 
							nMinDistTriId = nCurTriId;
						}
					}
				}
				else
				{			
					GetOutsidePointInTriangleStatus(x_cl, A, B, C, nStatus, x_out, NULL);
					return nCurTriId;	//we have found it!
				}				      				
			}

			//we need to continue
			for (int i = 0; i < 3; i++)
			{
				if (tCur.pnb[i] >= 0)
					queue.push(tCur.pnb[i]);
			}

			tCur.nMark = m_nCurrentMark;
		}

		//we did not find it in this stage, so use the  nMinDistTriId as the last resort		
		m_nCurrentMark++; 
		iStage++;
		if (nMinDistTriId < 0)
			queue.push(nTriId);
		else
		{
			queue.push(nMinDistTriId);
			iStage++;	//so we are at stage 2
		}		
	} //end while(true)
}

//------------------------------------------------------------------------
//Determine the mutual position of the given point x and the triangle ABC.
//The outcome of the test is returned in nStatus. The point x is assumed to lie on the
//plane defined by ABC and having ptlcOutside  status from GetPointInTriangleStatus.
//The method always succeeds returning in x_out the coordinates of the point that lies in the triangle ABC 
//and is the closest to x and the distance from x to x_out in dist variable. The method
//never returns nStatus to be ptlcOutside or ptlcInside.
//N.B. x_out and dist may be NULL, if the caller does not need the 
//extra information passed in these variables.
void vtkMAFPolyDataCutOutFilterEx::GetOutsidePointInTriangleStatus(double* x, double* A, double* B, double* C, 
	int& nStatus, double* x_out, double* dist )
	//------------------------------------------------------------------------
{
	double* ABC[3] = {A, B, C};

	//so it belongs to one of edges of the triangle (or its vertices)
	double dblEdgeDist;
	int iEdge = GetClosestEdge(x, ABC, &dblEdgeDist, x_out);
	if (iEdge >= 0)    
	{
		nStatus = ptlcEdge | iEdge;  
		if (dist != NULL)
			*dist = dblEdgeDist;
	}
	else
	{
		//it must be one of vertices
		double dblVertDist;
		int iVert = GetClosestVertex(x, ABC, &dblVertDist);    

		//OK, we are at vertex
		if (x_out != NULL) {  //if the closest point required, copy coordinates
			for (int j = 0; j < 3; j++){
				x_out[j] = ABC[iVert][j];
			}
		}

		nStatus = ptlcVertex | iVert;  
		if (dist != NULL)
			*dist = dblVertDist;
	}	
}

//------------------------------------------------------------------------
//Determine the mutual position of the given point x and the triangle ABC.
//The outcome of the test is returned in nStatus. Due to possible numerical inaccuracy, 
//	the test works with PointTolerance (set	this variable to zero, for accurate testing) and, 
//	therefore, it returns in x_out the coordinates of the point that lies in the triangle ABC 
//	and is the closest to x and the distance from x to x_out in dist variable. In a case 
//	that the point x is definitely outside the triangle, x_out and dist are filled but 
//	nStatus contains ptlcOutside status.
//	N.B. x_out and dist may be NULL, if the caller does not need the 
//	extra information passed in these variables
void vtkMAFPolyDataCutOutFilterEx::GetPointInTriangleStatus(double* x, 
	double* A, double* B, double* C, int& nStatus, 
	double* x_out, double* dist)
	//------------------------------------------------------------------------
{  
	double tol2 = this->PointTolerance*this->PointTolerance;
	double* ABC[3] = {A, B, C};

	//BES: 27.1.2009 - First, check for vertices
	double dblVertDist;
	int iVert = GetClosestVertex(x, ABC, &dblVertDist);    
	if (dblVertDist < tol2)
	{
		//OK, we are at vertex
		if (x_out != NULL) {  //if the closest point required, copy coordinates
			for (int j = 0; j < 3; j++){
				x_out[j] = ABC[iVert][j];
			}
		}

		nStatus = ptlcVertex | iVert;  
		if (dist != NULL)
			*dist = dblVertDist;

		return;
	}  

	//now, find the closest edge
	double dblEdgeDist;
	int iEdge = GetClosestEdge(x, ABC, &dblEdgeDist, x_out);
	if (iEdge >= 0 && dblEdgeDist < tol2)    
	{
		nStatus = ptlcEdge | iEdge;  
		if (dist != NULL)
			*dist = dblEdgeDist;
	}
	else
	{  
		if (PointInTriangle(x, A, B, C))
		{
			nStatus = ptlcInside;

			if (x_out != NULL) 
			{  
				//if the closest point required, copy coordinates
				for (int j = 0; j < 3; j++){
					x_out[j] = x[j];
				}
			}

			if (dist != NULL)
				*dist = 0.0;	//we are in triangle
		}
		else
		{
			nStatus = ptlcOutside; 
			if (dblEdgeDist < dblVertDist) 
			{
				//x_out is already filled
				if (dist != NULL)
					*dist = dblEdgeDist;	
			}
			else 
			{
				if (x_out != NULL) {  //if the closest point required, copy coordinates
					for (int j = 0; j < 3; j++){
						x_out[j] = ABC[iVert][j];
					}
				}

				if (dist != NULL)
					*dist = dblVertDist;	//we are in triangle
			}			
		}		
	}
}

//------------------------------------------------------------------------
//Determines which edge of the triangle ABC is the closest to the given point x.
//If x_out is not NULL, the projection of x onto the closest edge is stored in x_out.
//If dist is not NULL, the method returns in dist the distance from the closest point to x. 
//Returns -1, if there is no closest edge => the point is outside
int vtkMAFPolyDataCutOutFilterEx::GetClosestEdge(double *x, double *ABC[3], double* dist, double *x_out)
	//------------------------------------------------------------------------
{
	const static double eps_num = 1e-15;
	
	int iEdge = -1;
	double dblMinDist = DBL_MAX;
	for (int i = 0; i < 3; i++)
	{
		//the closest point on the edge can be computed as:
		//A + (u*(X - A))*u, where u is normalized vector B-A
		double u[3], v[3];
		for (int k = 0; k < 3; k++)
		{
			u[k] = ABC[(i + 1) % 3][k] - ABC[i][k];
			v[k] = x[k] - ABC[i][k];
		}

		double edge_len = vtkMath::Normalize(u);
		if (edge_len < eps_num)
			continue; //this is the collapsed edge! 

		double w = vtkMath::Dot(u, v);		

		//closest point lies on the line supporting the edge
		double x_cl[3];
		for (int k = 0; k < 3; k++) {
			x_cl[k] = ABC[i][k] + w*u[k];
		}

		double dblDist = DBL_MAX;
		if (w < -eps_num)	//eps_num is here for tolerance
		{
			//we are before the edge and so far there is no better edge
			if (iEdge < 0) {
				dblDist = vtkMath::Distance2BetweenPoints(ABC[i], x);	
			}
		}
		else if (w > edge_len + eps_num)
		{
			//we are after the edge and so far there is no better edge
			if (iEdge < 0) {
				dblDist = vtkMath::Distance2BetweenPoints(ABC[(i + 1) % 3], x);	
			}
		}
		else		
		{
			dblDist = vtkMath::Distance2BetweenPoints(x_cl, x);	
			if (dblDist < dblMinDist) {
				iEdge = i;
			}
		}

		//common update of x_out and dist
		if (dblDist < dblMinDist)
		{
			if (x_out != NULL)
			{
				//if the closest point required, copy coordinates
				for (int j = 0; j < 3; j++){
					x_out[j] = x_cl[j];
				}
			}

			dblMinDist = dblDist;				
		}		
	} //end for i (all edges)

	if (dist != NULL)
		*dist = dblMinDist;

	return iEdge;
}

//------------------------------------------------------------------------
//Determine which of points A, B, C is the closest to the given point x and returns its index. 
//If dist is not NULL, it also returns the closest distance.
int vtkMAFPolyDataCutOutFilterEx::GetClosestVertex(double* x, double** ABC, double* dist)
	//------------------------------------------------------------------------
{	
	double dblMinDist = vtkMath::Distance2BetweenPoints(x, ABC[0]);
	int iRetIndex = 0;

	for (int i = 1; i < 3; i++)
	{
		double dblDist = vtkMath::Distance2BetweenPoints(x, ABC[i]);
		if (dblDist < dblMinDist)
		{
			dblMinDist = dblDist;
			iRetIndex = i;
		}
	}

	if (dist != NULL)
		*dist = dblMinDist;

	return iRetIndex;
}

//------------------------------------------------------------------------
//Returns true, if the given point x lies inside (or on edge) of the give triangle, false, otherwise.
//N.B. x must lie on the plane defined by points ABC
//Barycentric coordinates are used to get accurate results
bool vtkMAFPolyDataCutOutFilterEx::PointInTriangle(double *x, double *A, double *B, double *C)
	//------------------------------------------------------------------------
{
	//See: http://www.blackpawn.com/texts/pointinpoly/default.html
	// Compute vectors
	double v0[3], v1[3], v2[3];
	for (int i = 0; i < 3; i++)
	{
		v0[i] = C[i] - A[i];
		v1[i] = B[i] - A[i];
		v2[i] = x[i] - A[i];
	}

	// Compute dot products
	double dot00 = vtkMath::Dot(v0, v0);
	double dot01 = vtkMath::Dot(v0, v1);
	double dot02 = vtkMath::Dot(v0, v2);
	double dot11 = vtkMath::Dot(v1, v1);
	double dot12 = vtkMath::Dot(v1, v2);

	// Compute barycentric coordinates
	double invDenom = 1.0 / (dot00 * dot11 - dot01 * dot01);
	double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
	double v = (dot00 * dot12 - dot01 * dot02) * invDenom;

	// Check if point is in triangle
	return ((u >= 0) && (v >= 0) && (u + v <= 1));
}

//------------------------------------------------------------------------
//Subdivides the given triangle into three new with common vertex nPtId. 
//Connectivity is updated as necessary.
void vtkMAFPolyDataCutOutFilterEx::SubdivideTriangle(int nTriId, int nPtId)
	//------------------------------------------------------------------------
{
	TRIANGLE tmpTri[2];
	TRIANGLE& t = m_Triangles[nTriId];
	for (int i = 0; i < 2; i++)
	{
		tmpTri[i].aVertex[0] = t.aVertex[i + 1];
		tmpTri[i].aVertex[1] = t.aVertex[(i + 2) % 3];
		tmpTri[i].aVertex[2] = nPtId;
		tmpTri[i].nMark = m_nCurrentMark;
	}

	t.aVertex[2] = nPtId;

	//update the connectivity
	int nTriangles = (int)m_Triangles.size();
	UpdateNeighbour(t.pnb[1], nTriId, nTriangles);
	UpdateNeighbour(t.pnb[2], nTriId, nTriangles + 1);

	tmpTri[0].pnb[0] = t.pnb[1];
	tmpTri[0].pnb[1] = nTriangles + 1;
	tmpTri[0].pnb[2] = nTriId;

	tmpTri[1].pnb[0] = t.pnb[2];
	tmpTri[1].pnb[1] = nTriId;
	tmpTri[1].pnb[2] = nTriangles;    

	t.pnb[1] = nTriangles;
	t.pnb[2] = nTriangles + 1;  

	//insert new triangles
	m_Triangles.push_back(tmpTri[0]);
	m_Triangles.push_back(tmpTri[1]);
}

//------------------------------------------------------------------------
//Subdivides the given triangle and its neighbor sharing edge iEdge 
//into two new triangles with common vertex nPtId. 
//Connectivity is updated as necessary.
void vtkMAFPolyDataCutOutFilterEx::Subdivide2Triangles(int nTriId, int iEdge, int nPtId)
	//------------------------------------------------------------------------
{ 
	TRIANGLE tmpTri[2];
	tmpTri[0].nMark = m_nCurrentMark;
	tmpTri[1].nMark = m_nCurrentMark;

	TRIANGLE* pTris[2];
	int nTriIds[2], iEdges[2];

	nTriIds[0] = nTriId;
	pTris[0] = &m_Triangles[nTriId];
	nTriIds[1] = pTris[0]->pnb[iEdge];    
	pTris[1] = (nTriIds[1] >= 0) ? &m_Triangles[nTriIds[1]] : NULL;

	iEdges[0] = iEdge;
	iEdges[1] = -1;

	//find the edge in the other triangle (if it exists)
	if (pTris[1] != NULL)
	{
		int ptHack = pTris[0]->aVertex[(iEdge + 1) % 3];
		for (int i = 0; i < 3; i++)
		{
			if (pTris[1]->aVertex[i] == ptHack) {
				iEdges[1] = i; break; //edge found
			}
		}

		if(!(iEdges[1] != -1)) return;
	}

	//subdivide both triangles
	int nTriangles = (int)m_Triangles.size();
	for (int iTri = 0; iTri < 2; iTri++)
	{ 
		if (pTris[iTri] == NULL)
			continue; //this triangle does not exist

		TRIANGLE& t1 = *pTris[iTri];
		TRIANGLE& t2 = tmpTri[iTri];
		iEdge = iEdges[iTri];
		nTriId = nTriIds[iTri];    

		int iv[3];
		for (int i = 0; i < 3; i++){
			iv[i] = t1.aVertex[(iEdge + i)  % 3];
		}

		t2.aVertex[0] = iv[1];
		t2.aVertex[1] = iv[2];
		t2.aVertex[2] = nPtId;

		t1.aVertex[0] = iv[2];
		t1.aVertex[1] = iv[0];
		t1.aVertex[2] = nPtId;

		//update connectivity
		t2.pnb[0] = t1.pnb[(iEdge + 1) % 3];
		t2.pnb[1] = nTriId;     //t1
		t2.pnb[2] = nTriIds[(iTri + 1) % 2];
		t1.pnb[0] = t1.pnb[(iEdge + 2) % 3];    
		t1.pnb[1] = nTriangles + 1 - iTri;
		t1.pnb[2] = nTriangles + iTri; //t2

		UpdateNeighbour(t2.pnb[0], nTriId, nTriangles + iTri);
	}

	//insert new triangles  
	if (pTris[1] == NULL)
	{
		pTris[0]->pnb[1] = -1;  //RELEASE NOTE: this call must be before push_back
		m_Triangles.push_back(tmpTri[0]);
	}
	else
	{
		m_Triangles.push_back(tmpTri[0]);
		m_Triangles.push_back(tmpTri[1]);
	}    
}

//------------------------------------------------------------------------
//Mark triangles to the right of the polyline.
void vtkMAFPolyDataCutOutFilterEx::MarkTriangles(vtkIdList* pIds)
	//------------------------------------------------------------------------
{    
	vtkIdType* pIdsPtr = pIds->GetPointer(0);
	int nPoints = pIds->GetNumberOfIds();

	//prepare fast cut locating routine  
	int nVertices = (int)m_Vertices.size();
	VERTEX_CUT_RING* pCutInfo = new VERTEX_CUT_RING[nVertices];
	memset(pCutInfo, 0, sizeof(VERTEX_CUT_RING)*nVertices);     

	//get the number of edges
	for (int i = 0; i < nPoints; i++)
	{
		//get next edge
		int iA = pIdsPtr[i];
		int iB = pIdsPtr[(i + 1) % nPoints];

		pCutInfo[iA].nAdjVerts++;
		pCutInfo[iB].nAdjVerts++;    
	}

	//allocated edges
	EDGE_DIR* pEdges = new EDGE_DIR[2*nPoints];
	int nTotalEntries = 0;

	for (int i = 0; i < nPoints; i++)
	{
		int iA = pIdsPtr[i];
		pCutInfo[iA].pAdjVerts = &pEdges[nTotalEntries];
		nTotalEntries += pCutInfo[iA].nAdjVerts;
		pCutInfo[iA].nAdjVerts = 0;
	}

	//and store edges
	for (int i = 0; i < nPoints; i++)
	{   
		//get next edge
		int iA = pIdsPtr[i];
		int iB = pIdsPtr[(i + 1) % nPoints];        

		pCutInfo[iA].pAdjVerts[pCutInfo[iA].nAdjVerts].nVertId = iB;
		pCutInfo[iA].pAdjVerts[pCutInfo[iA].nAdjVerts++].nOrient = 0;

		pCutInfo[iB].pAdjVerts[pCutInfo[iB].nAdjVerts].nVertId = iA;
		pCutInfo[iB].pAdjVerts[pCutInfo[iB].nAdjVerts++].nOrient = 1;    
	}


	bool bOrientationDetermined = false;
	bool bOrientationValid = true;     

	int nCurTriId = 0;
	std::queue< int > queue;
	queue.push(nCurTriId);

	m_nCurrentMark++;
	while (!queue.empty())
	{
		nCurTriId = queue.front();
		queue.pop();

		TRIANGLE& tCur = m_Triangles[nCurTriId]; 
		if (tCur.nMark == m_nCurrentMark)
			continue; //already processed, skip it

		for (int i = 0; i < 3; i++)
		{
			if (tCur.pnb[i] < 0)
				continue; //there is no neighbor, OK

			//check, if we can cross the current edge
			//(it cannot be crossed, if it is cut edge
			bool bOK = true;      
			int iA = tCur.aVertex[i];
			int iB = tCur.aVertex[(i + 1) % 3];      
			for (int j = 0; j < pCutInfo[iA].nAdjVerts; j++)
			{
				if (iB == pCutInfo[iA].pAdjVerts[j].nVertId) 
				{                              
					if (!bOrientationDetermined)
					{
						//BES 28.1.2009 - new method for orientation check 

						//The idea is to determine whether the far point (not lying on
						//this edge) of the current triangle lies to the left of edge
						//iA, iB or to the right measured in the plane defined by the
						//triangle. The orientation flag can be then easily determined
						//according to CutSide parameter
						int iC = tCur.aVertex[(i + 2) % 3]; 

						//get points
						double* A = m_Vertices[iA].dCoord;
						double* B = m_Vertices[iB].dCoord;
						double* C = m_Vertices[iC].dCoord;

						double u[3], v[3], w[3], n[3];                        
						for (int k = 0; k < 3; k++){
							u[k] = B[k] - A[k]; v[k] = C[k] - A[k];              
						}

						vtkMath::Cross(u, v, w);
						if (vtkMath::Norm(w) >= 1e-15) //eps (degenerated triangles)
						{
							//regular triangle, u and w define plane, check, if C is left
							//i.e., in the positive half space
							vtkMath::Cross(w, u, n);
							double dblDot = vtkMath::Dot(v, n);

							//if the cut edge A,B is oriented from A to B, i.e., nOrient == 0, then              
							bOrientationValid = (dblDot > 0.0) == (this->CutOutSide == 0);

							//otherwise we need to invert it
							if (pCutInfo[iA].pAdjVerts[j].nOrient != 0)
								bOrientationValid = !bOrientationValid;

							bOrientationDetermined = true;
						}            
					}

					bOK = false;           
					break;
				}
			}

			if (bOK)  //if we are allowed to cross the edge
				queue.push(tCur.pnb[i]);      
		}    

		tCur.nMark = m_nCurrentMark;  //marked
	} //end while (all triangles in queue)


	//mark every point
	int nLastMark = m_nCurrentMark;
	if (!bOrientationValid) {
		m_nCurrentMark++;
	}


	//we have now marked all triangles to be removed
	//=> perform inversion of marking
	int nTriangles = (int)m_Triangles.size();
	for (int i = 0; i < nTriangles; i++)
	{
		TRIANGLE& t = m_Triangles[i];
		if ((t.nMark == nLastMark) == bOrientationValid)
		{
			t.nMark = m_nCurrentMark;
			for (int j = 0; j < 3; j++)
			{
				VERTEX& v = m_Vertices[t.aVertex[j]];
				v.nMark = m_nCurrentMark;
			}
		}
	}

	delete[] pCutInfo;
	delete[] pEdges;
}
