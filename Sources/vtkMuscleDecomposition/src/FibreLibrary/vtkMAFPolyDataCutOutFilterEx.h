/*========================================================================= 
Program: Multimod Application Framework RELOADED 
Module: $RCSfile: vtkMAFPolyDataCutOutFilterEx.h,v $ 
Language: C++ 
Date: $Date: 2012-03-19 12:45:22 $ 
Version: $Revision: 1.1.2.2 $ 
Authors: Josef Kohout
========================================================================== 
Copyright (c) 2008 University of Bedfordshire (www.beds.ac.uk), 
Copyright (c) 2011University of West Bohemia (www.zcu.cz)
See the COPYINGS file for license details 
=========================================================================
vtkMAFPolyDataCutOutFilterEx cuts the input polydata surface
by the cut defined by the input polyline (see below) splitting
triangles affected by the cut into two polygons (which are
later triangulated). The connectivity between vertices of these
two polygons is severed (thus it is physical cut). By cutting,
the surface, therefore, splits into two parts (it is assumed
that the polyline is closed). One part (to the CutOutSide from the 
polyline) is discarded (cut out).

The filter has some strong assumptions on the input polyline:
1) polyline must be continuous (i.e., it must have just one 
component; it is one "contour" only) and closed
2) polyline edges must be oriented (i.e., for instance cells 
can be given as 01,12,23,30 but not 01,32,03,12 - which is
typical output from vtkCutter). The algorithm cuts out the part
of surface that lies to the left or right (see CutOutSide) from the polyline .
3) all points of polyline must lie on the input surface (within 
specified m_AbsCellTolerance and m_AbsPointTolerance tolerances)
4) two neighbouring points of the polyline (i.e., they form 
an edge) must lie in the same triangle of surface (i.e., on
the same plane)

Another assumption is that the input polydata contains of
triangles only and these triangles are consistently oriented.

Typical use of this filter is to cut out unwanted parts of 
surface that are identified by a polyline obtained from vtkCutter
based filters (beware that vtkCutter outputs violate 
assumptions 1 and 2) or filters for finding the shortest
path (in this case all conditions are usually fullfiled)


Changes
-------
29.10.2009 - the contour points must be either vertices of the mesh 
or they must lie on edges, SubdivideTriangle is no longer supported.
Contour is no longer checked, if it lies on the surface. It simply
must lie there and the closest points are picked. 

17.10.2011 - SubdivideTriangle reintroduced, the contour points may be 
distant (the path between these points is automatically constructed, so it can be cut)*/

#ifndef __vtkMAFPolyDataCutOutFilterEx_h
#define __vtkMAFPolyDataCutOutFilterEx_h

//----------------------------------------------------------------------------
// Include:
//----------------------------------------------------------------------------
#pragma once

#pragma warning(push)
#pragma warning(disable:4996)

#include "vtkObject.h"
#if VTK_MAJOR_VERSION < 5
#include "vtkPolyDataToPolyDataFilter.h"
#define vtkPolyDataAlgorithm vtkPolyDataToPolyDataFilter
#else
#include "vtkPolyDataAlgorithm.h"
#if VTK_MAJOR_VERSION < 7
typedef unsigned long vtkMTimeType;
#endif
#endif

#pragma warning(pop)

#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include <vector>

//#include "vtkLHPConfigure.h"
//#include "mafDllMacros.h"

#define DIJKSTRA_CUT //produces noisy cut but works OK even in cases of folded surface (unlike the other method)

//----------------------------------------------------------------------------
// forward declarations :
//----------------------------------------------------------------------------
class vtkIdList;
class vtkCamera;
class vtkRenderer;
class vtkActor;
class vtkMAFVisualDebugger;

//----------------------------------------------------------------------------
// vtkMAFPolyDataCutOutFilterEx class 
//----------------------------------------------------------------------------


class VTK_EXPORT vtkMAFPolyDataCutOutFilterEx : public vtkPolyDataAlgorithm
{
protected:
#pragma region Nested Structures and Classes	
	typedef struct TRIANGLE
	{
		int   aVertex[3];   //indices of vertices
		int   pnb[3];       //indices of neighbouring triangles
		int   nMark;        //mark for various things
	} TRIANGLE;

	typedef struct VERTEX
	{
		double  dCoord[3];  //coordinates
		int     nMark;      //mark for various things
#ifdef DIJKSTRA_CUT
		struct {
			int nDMark;				//mark for Dijskra
			int nTriId;				//index of triangle containing the vertex
			int iVertexPos;		//position of the vertex in the triangle

			double dSum;				//distance cost
			int nPrevVertexId;	//previous vertex Id

//			int iNextPQ;				//next vertex with lower priority (or -1, if this is the last one)
//			int iPrevPQ;					//previous vertex with higher priority
			int iPQIndex;			//index to position of this vertex to PQ
		} Dijstra;
#endif
	} VERTEX;

	typedef struct EDGE_DIR
	{
		int nVertId;   //adjacent vertex
		int nOrient;   //orientation of the edge
	} EDGE_DIR;

	typedef struct VERTEX_CUT_RING
	{
		int nAdjVerts;        //number of stored adj. vertices (0, if the point is normal, 2 if it is cut)
		EDGE_DIR* pAdjVerts;  //adjacent vertices (with highlighted orientation)
		//int pAdjVerts[2];   //indices of adjacent cut vertices

		//debug
		//bool bDetected[2];
	} VERTEX_CUT_RING;

	typedef struct POINTLOC
	{
		int nBestTriId;		//the currently closest triangle to the point
		int nEdgeId;			//edge closest to the point
		bool bOnEdge;			//true, if the closest point lies directly on the edge
		//n.b., in this case dblTime must be <0..1>

		double dblTime;		//position of the closest point on the line supported by the edge
		double dblDist;		//the distance between the closest point on the EDGE and the given point

	} POINTLOC;

	//enum with possible mutual position of the point and the simplex
	//suitable up to d (dimension) = 15
	//NOTE: HIWORD denotes the kind of position, LOWORD the exact ID of edge, etc.  
	enum PointLocation
	{
		ptlcInside	= 0x00000,	  //Point lies inside the simplex, e.g., inside the triangle or tetrahedron
		ptlcVertex	= 0x10000,	  //Point lies in vertex m_Value & ~Vertex
		ptlcEdge	= 0x20000,	    //Point lies on edge m_Value & ~Edge      
		ptlcOutside = 0x80000000,	//Point lies outside			    
	};

#pragma endregion Nested Structures and Classes
	
public:
	static vtkMAFPolyDataCutOutFilterEx *New();

	vtkTypeMacro(vtkMAFPolyDataCutOutFilterEx, vtkPolyDataAlgorithm);

protected:
	vtkMAFPolyDataCutOutFilterEx();           
	virtual ~vtkMAFPolyDataCutOutFilterEx();

protected:
	vtkPolyData* CuttingPolyline;         //<polyline used as cutter
	int CutOutSide;                       ///<0 (by default), if the the mesh left from polyline is to be clipped, 1 otherwise

	vtkPolyData* OutputCuttingPolyline;		///<NULL by default, output polyline
	vtkPolyData* OutputClippedPart;				///<NULL by default, the discarded (clipped) part of the mesh
	
	int Debug;														//<1, if the debug mode is ON (displays cutting progress)	

	double PointTolerance;                //<the minimal relative distance between two points (on a unit edge), closer points are considered to be one 
	
#pragma warning(suppress: 4251)					//This is to suppress the warning that we are unable to access these members from a client, if this class is exported from a DLL, since it is protected and in fact will not be accessed from the client 
	std::vector< VERTEX > m_Vertices;     //<vertices of the surface
#pragma warning(suppress: 4251)
	std::vector< TRIANGLE > m_Triangles;  //<triangles of the surface  

	unsigned long m_LastExecuteTimeStamp; //<timestamp for last succesful execute
	int m_nCurrentMark;                   //<current mark to be used  

#ifdef DIJKSTRA_CUT
	int m_nCurrentDijkstraMark;						//<current mark to be used for Dijsktra	
#endif

public:
	/** Gets the current cutting polyline. */
	vtkGetObjectMacro(CuttingPolyline, vtkPolyData);

	/** Sets a new cutting polyline. */
	vtkSetObjectMacro(CuttingPolyline, vtkPolyData);
	
	/** Gets the output cutting polyline. 
	N.B. by default, this is NULL unless non-NULL reference was passed to SetOutputCuttingPolyline
	prior to the call of Update method (then the object is updated)*/
	vtkGetObjectMacro(OutputCuttingPolyline, vtkPolyData);

	/** Sets a new object where output cutting polyline should be stored. */
	vtkSetObjectMacro(OutputCuttingPolyline, vtkPolyData);

	/** Gets the clipped part. 
	N.B. by default, this is NULL unless non-NULL reference was passed to SetOutputClippedPart 
	prior to the call of Update method (then the object is updated). */
	vtkGetObjectMacro(OutputClippedPart, vtkPolyData);

	/** Sets a new object where clipped part should be stored. */
	vtkSetObjectMacro(OutputClippedPart, vtkPolyData);

	/** Gets the current minimal distance for two different points. 
	The value is given in a fraction of unit edge. */
	vtkGetMacro(PointTolerance, double);

	/** Sets a new minimal distance for two different points. 
	Points closer than this tolerance are considered to be one. 
	The value is given in a fraction of unit edge. */
	vtkSetMacro(PointTolerance, double);

	/** Gets the current settings of which part of mesh to cut out.   
	0 = mesh to the left from the cutting polyline is to be discarded,
	1 = mesh to the right from the cutting polyline is to be discarded.*/
	vtkGetMacro(CutOutSide, int);

	/** Sets which part of mesh to cut out. 
	0 = mesh to the left from the cutting polyline is to be discarded,
	1 = mesh to the right from the cutting polyline is to be discarded.*/
	vtkSetMacro(CutOutSide, int);

	/** Gets the current settings of which part of mesh to cut out.   
	0 = mesh to the left from the cutting polyline is to be discarded,
	1 = mesh to the right from the cutting polyline is to be discarded.*/
	vtkGetMacro(Debug, int);

	/** Sets which part of mesh to cut out. 
	0 = mesh to the left from the cutting polyline is to be discarded,
	1 = mesh to the right from the cutting polyline is to be discarded.*/
	vtkSetMacro(Debug, int);

public:
	/** Return this object's modified time. */  
	/*virtual*/ vtkMTimeType GetMTime();

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
		vtkInformationVector** inputVector, vtkInformationVector* outputVector);
#endif

	/** Creates internal structures for mesh */
	virtual void InitMesh(vtkPolyData* input);

	/** Saves marked parts of mesh (or the whole mesh, if bExtractMarkedPartOnly is false) into output.*/
	virtual void DoneMesh(vtkPolyData* output, bool bExtractMarkedPartOnly = true);

	/** Saves the polyline in pIdsIns into the output*/
	virtual void DoneCuttingPolyline(vtkPolyData* output, vtkIdList* pIdsIns);

	/** releases internal structures for mesh*/
	virtual void ReleaseMesh();

protected:
	/** Inserts all points from the cutting polyline into mesh.
	The routine returns indices of inserted points. Indices are ordered to form 
	oriented polyline. If something goes wrong (e.g., the input polyline is 
	incorrect, the routine returns false; true otherwise. The algorithm is as
	follows. If the point to be inserted corresponds to already existing mesh
	point, it is not inserted again but the index of mesh point is stored,
	otherwise the point is added to the list of vertices. If the point lies
	inside a triangle, this triangle is subdivided into three new triangles
	(one is replaced, two newly constructed). If the point lies on an edge,
	both triangles sharing this edge are subdivided into two new triangles
	(one is replaced, one newly constructed).  */
	bool InsertCutPoints(vtkIdList* pIds);

	/** Sorts all points from the cutting polyline to form a closed contour.
	Returns false, if the polyline is invalid; otherwise true and in pIds
	are ordered indices of points. */
	bool SortPolylinePoints(vtkIdList* pIds);  

	/** Locates the given point in the current mesh.
	It searches for the triangle that contains the given point x in the mesh.
	The indentifier of this triangle is returned and the mutual position of 
	the point and the triangle is returned in nStatus, the coordinates
	to be stored for the point (might be shifted in comparison to x) are
	returned in x_out variable. When the method is called for the first time, 
	nTriId must be -1 (nStatus is ignored) and the algorithm performs 
	a brute-force search (using VTK to speed up the process). The successive
	calls have nTriId set to the identifier returned from the previous call
	and nStatus unchanged. Only the triangles around the last processed
	points are checked.	The method supposes that the caller modifies
	the triangulation regarding to the returned values.
	N.B. This method never returns ERROR. */
	int LocatePoint(double* x, int nTriId, int& nStatus, double* x_out);

	/**
	returns the position of the not shared vertex V of the simplex pSimplexAdj
	that shares the edge nEdge with the given simplex pSimplex	*/
	inline int GetFarVertexPosition(const TRIANGLE* pSimplex, const TRIANGLE* pSimplexAdj, int nEdge)
	{
		//This implementation is a slightly faster than traditional for
		unsigned int nEdgeID = pSimplex->aVertex[nEdge] + pSimplex->aVertex[(nEdge + 1) % 3];	
		if (pSimplexAdj->aVertex[0] + pSimplexAdj->aVertex[1] == nEdgeID)	
			return 2;

		if (pSimplexAdj->aVertex[0] + pSimplexAdj->aVertex[2] == nEdgeID)	
			return 1;

		return 0;
	}

	/** Determine the mutual position of the given point x and the triangle ABC.
	The outcome of the test is returned in nStatus. Due to possible numerical inaccuracy, 
	the test works with PointTolerance (set	this variable to zero, for accurate testing) and, 
	therefore, it returns in x_out the coordinates of the point that lies in the triangle ABC 
	and is the closest to x and the distance from x to x_out in dist variable. In a case 
	that the point x is definitely outside the triangle, x_out and dist are filled but 
	nStatus contains ptlcOutside status.
	N.B. x_out and dist may be NULL, if the caller does not need the 
	extra information passed in these variables.*/
	void GetPointInTriangleStatus(double* x, double* A, double* B, double* C, 
		int& nStatus, double* x_out, double* dist );

	/** Determine the mutual position of the given point x and the triangle ABC.
	The outcome of the test is returned in nStatus. The point x is assumed to lie on the
	plane defined by ABC and having ptlcOutside  status from GetPointInTriangleStatus.
	The method always succeeds returning in x_out the coordinates of the point that lies in the triangle ABC 
	and is the closest to x and the distance from x to x_out in dist variable. The method
	never returns nStatus to be ptlcOutside or ptlcInside.
	N.B. x_out and dist may be NULL, if the caller does not need the 
	extra information passed in these variables.*/
	void GetOutsidePointInTriangleStatus(double* x, double* A, double* B, double* C, 
		int& nStatus, double* x_out, double* dist );

	/** Determine which of points A, B, C is the closest to the given point x and returns its index. 
	If dist is not NULL, it also returns the closest distance.*/
	int GetClosestVertex(double* x, double** ABC, double* dist);
	
	/** Determines which edge of the triangle ABC is the closest to the given point x.
	If x_out is not NULL, the projection of x onto the closest edge is stored in x_out.
	If dist is not NULL, the method returns in dist the distance from the closest point to x. */
	int GetClosestEdge(double *x, double *ABC[3], double* dist, double *x_out);

	/** Returns true, if the given point x lies inside (or on edge) of the give triangle, false, otherwise.
	N.B. x must lie on the plane defined by points ABC
	Barycentric coordinates are used to get accurate results*/
	bool PointInTriangle(double *x, double *A, double *B, double *C);

	/** Subdivides the given triangle into three new with common vertex nPtId. 
	Connectivity is updated as necessary. */
	void SubdivideTriangle(int nTriId, int nPtId);

	/** 
	Subdivides the given triangle and its neighbor sharing edge iEdge 
	into two new triangles with common vertex nPtId. 
	Connectivity is updated as necessary. */
	void Subdivide2Triangles(int nTriId, int iEdge, int nPtId);

	/** Updates the link nOldNeib in triangle nTriId by nNewNeib.  */
	inline void UpdateNeighbour(int nTriId, int nOldNeib, int nNewNeib);


	/** Finds the triangle containing vertex with nVertexIndex index.
	The method performs bread-first search starting with the triangle identified 
	by nStartTriId and returns id of the located triangle and in iVertexPos also 
	which vertex of the triangle is the located one.*/
	int FindTriangle(int nStartTriId, int nVertexIndex, int& iVertexPos);

	/** Finds the first triangle in the fan around the vertex identified by triangle id and vertex pos. 
	N.B. use FindNextTriangleInFan to traverse to the next triangle. */
	inline TRIANGLE* FindFirstTriangleInFan(int& nCurTriId, int& iVertexPos);

	/** Finds the first triangle in the fan around the vertex identified by triangle id and vertex pos. 
	Input parameters are those returned by the previous call of FindFirstTriangleInFan
	or FindNextTriangleInFan. When the method returns nCurTriId contains id of the 
	returned triangle and iVertexPos the position of vertex around which the 
	search was initiated (by FindFirstTriangleInFan).
	N.B. the loop is completed when the method returns the same pointer as  FindFirstTriangleInFan.. */
	inline TRIANGLE* FindNextTriangleInFan(const TRIANGLE* pSimplex, int& nCurTriId, int& iVertexPos);

private:
	/** Finds the first triangle in the fan around the vertex identified by triangle id and vertex pos. 
	Input parameters are those returned by the previous call of FindFirstTriangleInFan
	or FindNextTriangleInFan. When the method returns nCurTriId contains id of the 
	returned triangle and iVertexPos the position of vertex around which the 
	search was initiated (by FindFirstTriangleInFan). It supports open fans.
	N.B. the loop is completed when the method returns the same pointer as  FindFirstTriangleInFan. */
	inline TRIANGLE* FindNextTriangleInOpenFan(const TRIANGLE* pSimplex, int& nCurTriId, int& iVertexPos);

protected:
#ifdef DIJKSTRA_CUT
	/** Finds the shortest path from vertex nFromIndex to nToIndex.
	All indices on the path except for the last (nToIndex) are inserted into pContourIds and vertices are marked.
	nTriId is a helper ID to a valid triangle close to the first point to speed-up the process. 
	Returns ID to a valid triangle close to nToIndex to speed-up the further process.*/
	int CreateDijskraCutPath(int nFromIndex, int nToIndex, vtkIdList* pContourIds, int nTriId = 0, int* pPQBuf = NULL);
	
	/** Deletes the minimum from the PQ. 
	Returns the new index of vertex at head of the queue. */
	int DijkstraPQDeleteMin(int* pPQ, int& nPQSize); // nPQHeadIndex);

	/** Inserts the new vertex into the PQ.
	N.B. the value of the vertex must be already specified.
	Returns the new index of vertex at head of the queue. */
	void DijkstraPQInsert(int nIndex, int* pPQ, int& nPQSize);

	/** Moves vertex with the given nIndex in the direction to head of Dijkstra priority queue until the proper place is found.
	Returns the new index of vertex at head of the queue. */
	void DijkstraPQDown(int nIndex, int* pPQ, int nPQSize);

	/** Moves vertex with the given nIndex in the direction to tail of Dijkstra priority queue until the proper place is found.
	Returns the new index of vertex at head of the queue. */
	void DijkstraPQUp(int nIndex, int* pPQ);

	inline void DijkstraPQSwap(int iPQ, int jPQ, int* pPQ) 
	{
		int nTmp = pPQ[iPQ];
		pPQ[iPQ] = pPQ[jPQ];
		pPQ[jPQ] = nTmp;

		m_Vertices[pPQ[iPQ]].Dijstra.iPQIndex = iPQ;
		m_Vertices[pPQ[jPQ]].Dijstra.iPQIndex = jPQ;
	}	
#else
	/** Finds the triangle that contains the vertex nIndex.
	The search is done in the fan of triangles around the vertex identified by nCurTriId and iVertexPos.
	Both nCurTriId and iVertexPos are update to contain the located vertex (or the last tested one) 
	The method returns NULL, if triangle has not been located. */
	TRIANGLE* FindTriangleInFanV(int nIndex, int& nCurTriId, int &iVertexPos);

	/** Finds the mesh edge intersected (in projection) by the line going from the vertex FanVertIndex towards the vertex nIndex.
	The method processes triangles in the fan of triangles around the vertex FanVertIndex identified by 
	nCurTriId and iVertexPos parameters. N.B. these parameters are updated as the method proceeds.
	When the method finds the triangle with an edge intersected by the line FanVertIndex - nIndex projected onto the plane 
	of this triangle, it determines the status of intersection (ptlcEdge, ptlcVertex) and computes the 
	coordinates of intersection.  The method may return NULL, if due to some numeric inaccuracy 
	something went terribly wrong, otherwise the triangle containing x_out point.  */
	TRIANGLE* FindTriangleInFanE(int nIndex, int& nCurTriId, int &iVertexPos, int& nStatus, double* x_out);

	/** Finds the mesh vertex closest to the ray going towards the vertex nIndex (or in the given direction).
	The method processes triangles in the fan of triangles around the vertex FanVertIndex identified by 
	nCurTriId and iVertexPos parameters. N.B. these parameters are updated as the method proceeds.
	When the method exits, nCurTriId and iVertexPos contain the located vertex
	The method never returns NULL.*/
	TRIANGLE* FindTriangleInFanL(int nIndex, int& nCurTriId, int &iVertexPos);

	/** Finds all triangles sharing the vertex FanVertIndex identified by nCurTriId and iVertexPos parameters.
	pOutput is allocated in the method and the caller is responsible for its release when the data is no longer needed.
	The method returns the number of retrieved triangles.
	N.B. triangles are ordered. */
	int GetTriangleFan(int nCurTriId, int iVertexPos, TRIANGLE**& pOutput);
#endif
	/** 
	Updates the contour in pIds list in such a way that points form a connected chain (connected by
	mesh edges) and there are no duplicities in the contour (this is fixed). 
	nTriId is a helper ID to a valid triangle close to the first point to speed-up the process. */
	void CreateCutPointsNonIntersectingContour(vtkIdList* pIds, int nTriId = 0);

	/** Updates the contour in pIds list in such a way that points forms a connected chain (connected by
	mesh edges)  connected but does not guarantee that the contour will be manifold.
	nTriId is a helper ID to a valid triangle close to the first point to speed-up the process. 
	The method is called by  CreateCutPointsNonIntersectingContour. */
	void CreateCutPointsContour(vtkIdList* pIds, int nTriId = 0);

	/** Inverts the curve stored in pInvertPtr at indices from iLeft to iRight (inclusively)*/
	void InvertCurve(vtkIdType* pInvertPtr, int iLeft, int iRight);

	/** Mark triangles to the right of the polyline. */
	void MarkTriangles(vtkIdList* pIds);

#pragma region DEBUG Visualization
	/** Debug visualization of CreateCutPointsContour.	
	Visualize the mesh with highlighted points to be inserted (pInsertedPts)
	and already connected points (pCountourPts). ptIdFrom and ptIdTo are
	indices of points currently being connected. Called only, if Debug != 0 */
	void Debug_Visualize_Progress(vtkMAFVisualDebugger* vd, const vtkIdList* pInsertedPts, 
		const vtkIdList* pContourPts, int ptIdFrom, int ptIdTo);

	/** Debug visualization of InsertCutPoints.
	Visualize the mesh with highlighted points already processed (pProcessed, nProcessed)
	and input poly-line contour sorted as given in pInputContourPts. iCurId is the point of
	the input contour to be inserted in the current iteration.	Called only, if Debug != 0 */
	void Debug_Visualize_Progress2(vtkMAFVisualDebugger* vd, const vtkIdList* pInputContourPts, int iCurId, 
		const vtkIdType* pProcessed, int nProcessed);

	/** Debug visualization of correcting the contour.
	Visualize the mesh with highlighted intersecting contour pInputContourPts 
	and the  current non-intersecting contour  pProcessed. iCurId is the point of pInputContourPts
	where the contour is intersecting (or touching) itself.	Called only, if Debug != 0 */
	void Debug_Visualize_Progress3(vtkMAFVisualDebugger* vd, const vtkIdList* pInputContourPts, int iCurId, const vtkIdList* pProcessed);
#pragma endregion	
private:
	vtkMAFPolyDataCutOutFilterEx(const vtkMAFPolyDataCutOutFilterEx&);  // Not implemented.
	void operator = (const vtkMAFPolyDataCutOutFilterEx&);  // Not implemented.  
};

#pragma region INLINES
//------------------------------------------------------------------------
//Updates the link nOldNeib in triangle nTriId by nNewNeib.
inline void vtkMAFPolyDataCutOutFilterEx::UpdateNeighbour(int nTriId, int nOldNeib, int nNewNeib)
	//------------------------------------------------------------------------
{
	if (nTriId < 0)
		return;

	TRIANGLE& t = m_Triangles[nTriId];
	for (int i = 0; i < 3; i++)
	{
		if (t.pnb[i] == nOldNeib) {
			t.pnb[i] = nNewNeib; return;
		}
	}

	//error
}

//------------------------------------------------------------------------
/** Finds the first triangle in the fan around the vertex identified by triangle id and vertex pos. 
N.B. use FindNextTriangleInFan to traverse to the next triangle. */
inline vtkMAFPolyDataCutOutFilterEx::TRIANGLE* 
	vtkMAFPolyDataCutOutFilterEx::FindFirstTriangleInFan(int& nCurTriId, int& iVertexPos)
	//------------------------------------------------------------------------
{
	return &m_Triangles[nCurTriId];
}

//------------------------------------------------------------------------
/** Finds the first triangle in the fan around the vertex identified by triangle id and vertex pos. 
Input parameters are those returned by the previous call of FindFirstTriangleInFan
or FindNextTriangleInFan. When the method returns nCurTriId contains id of the 
returned triangle and iVertexPos the position of vertex around which the 
search was initiated (by FindFirstTriangleInFan).
N.B. the loop is completed when the method returns the same pointer as  FindFirstTriangleInFan. */
inline vtkMAFPolyDataCutOutFilterEx::TRIANGLE* 
	vtkMAFPolyDataCutOutFilterEx::FindNextTriangleInFan(const TRIANGLE* pSimplex, int& nCurTriId, int& iVertexPos)
	//------------------------------------------------------------------------
{
	if (pSimplex->pnb[iVertexPos] < 0)
		return FindNextTriangleInOpenFan(pSimplex, nCurTriId, iVertexPos);
	else
	{
		TRIANGLE* pNbSimplex = &m_Triangles[nCurTriId = pSimplex->pnb[iVertexPos]];
		int iNbVxPos = GetFarVertexPosition(pSimplex, pNbSimplex, iVertexPos);				 		
		iVertexPos = (iNbVxPos + 2) % 3;
		return pNbSimplex;
	}
}

#pragma endregion

#if VTK_MAJOR_VERSION < 5
#undef vtkPolyDataAlgorithm
#endif
#endif