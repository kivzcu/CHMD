/*=========================================================================
Module:    $RCSfile: vtkMuscleDecomposer.cxx $
Language:  C++
Date:      $Date: 2011-05-03 21:10 $
Version:   $Revision: 0.1.0.0 $
Author:    David Cholt
Notes:
=========================================================================*/

// Blender: Forward: Y Forward, Up: Z Up, sel only, objects as objects

//----------------------------------------------------------------------------
// Include:
//----------------------------------------------------------------------------
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable: 4996)	// warning C4996: 'strcpy': This function or variable may be unsafe.
#endif
#include "vtkMuscleDecomposer.h"
#include "vtkObjectFactory.h"

#include "vtkMuscleDecomposerSplineInterpolator.h"
#include "vtkLHPPolyLinesBuilder.h"
#include "vtkLHPMuscleFibresResample.h"
#include "vtkPolyDataWriter.h"
#include "vtkPointData.h"

#include <assert.h>

//_RPT0 linux wrapper
#ifdef __unix__
#define _RPT0(...) fprintf(stderr,"ASSERT 0")
#define _RPT1(...) fprintf(stderr,"ASSERT 1")
#endif

//----------------------------------------------------------------------------
#ifdef _MSC_VER
#pragma warning(pop)
#endif

#define debug // comment for release!!

void my_assert(bool expr)
{
	if (expr)
		return;

	_RPT0(_CRT_WARN, "ASSERT");
}

#define _CHECK_POINT(x) //my_assert(_finite((x)[0]) && _finite((x)[1]) && _finite((x)[2]))

vtkStandardNewMacro(vtkMuscleDecomposer);

//// =========================================================================
//// Nested classes
//// =========================================================================

//// FIBER

//----------------------------------------------------------------------------
vtkMuscleDecomposer::Fiber::Fiber()
//----------------------------------------------------------------------------
{
	DLayer = 0;
}

//----------------------------------------------------------------------------
vtkMuscleDecomposer::Fiber::~Fiber()
//----------------------------------------------------------------------------
{
	PointIDs.clear();
}

//// POINT

//----------------------------------------------------------------------------
vtkMuscleDecomposer::Point::Point(double x, double y, double z)
//----------------------------------------------------------------------------
{
	DCoord[0] = x; //x
	DCoord[1] = y; //y
	DCoord[2] = z; //z

	AttachmentConnectionID = -1;
	Interpolated = false;
	DistToNext = 0;
	// by default, all points are anchors
	LeftAnchorID = -1;
	RightAnchorID = -1;
	Splits = false;

	MuscleTendonInfo = -1;
}

//----------------------------------------------------------------------------
vtkMuscleDecomposer::Point::~Point()
//----------------------------------------------------------------------------
{
}

//----------------------------------------------------------------------------
// Computes coordinates of point on interpolated curve, using given parameter.
// Substracting the segment start t from point's t is neccessary for some
// interpolation methods.
void vtkMuscleDecomposer::Point::ComputeInterpolatedCurveCoord(double t, double* result, bool Subtract)
//----------------------------------------------------------------------------
{
	//double* result = new double[3];
	if (!this->Interpolated)
	{
		result[0] = result[1] = result[2] = 0;
		return;
	}
	double h = t;
	if (Subtract)
	{
		h = (t - this->t);
	}
	double u = 0;

	for (int axis = 0; axis < 3; axis++)
	{
		//u = a[axis] * h * h * h + b[axis] * h * h + c[axis] * h + d[axis];
		u = d[axis] + h*(c[axis] + h* (b[axis] + h* (a[axis]))); // horner's schema
		result[axis] = u;
	}
}

//----------------------------------------------------------------------------
// Assigns normal to this Point
//----------------------------------------------------------------------------
void vtkMuscleDecomposer::Point::SetNormal(double* normal)
//----------------------------------------------------------------------------
{
	this->DNormal[0] = normal[0];
	this->DNormal[1] = normal[1];
	this->DNormal[2] = normal[2];
}

//// =========================================================================
//// vtkMuscleDecomposer class
//// =========================================================================

//----------------------------------------------------------------------------
// vtkMuscleDecomposer Constructor
vtkMuscleDecomposer::vtkMuscleDecomposer()
//----------------------------------------------------------------------------
{
	InterpolationMethod = SPLINE_CATMULL_ROM;
	InterpolationSubdivision = 10;
	SurfaceSubdivision = 5;
	ArtifactEpsilon = 10.0;
	DoFlipNormals = true;
	DoOffsetFibers = true;
	DoForceFlipNormals = false;
	DoSpreadFibers = true;
	DoConnectFibersToAA = true;
	FiberCount = 0;
	PennationTreshold = vtkMath::Pi() * 4 / 5;
	tree = vtkOBBTree::New();
	tree->SetTolerance(0.0);

	this->NormalMesh = NULL;

#if VTK_MAJOR_VERSION >= 6
	//The filter is implemented skipping the default Input -> Output path, so actually it is a source
	this->SetNumberOfInputPorts(0);
#endif
}

//----------------------------------------------------------------------------
// vtkMuscleDecomposer Destructor
vtkMuscleDecomposer::~vtkMuscleDecomposer()
//----------------------------------------------------------------------------
{
	ClearMesh();

	for (int side = 0; side < 2; side++)
	{
		for (int i = 0; i < (int)AttachmentsSurfaces[side].size(); i++)
		{
			if (AttachmentsSurfaces[side][i])
			{
				AttachmentsSurfaces[side][i]->Delete();
				AttachmentsSurfaces[side][i] = NULL;
			}
		}
		AttachmentsSurfaces[side].clear();
	}

	int nFibs = (int)ControlFibers.size();
	for (int i = 0; i < nFibs; i++)
	{
		ControlFibers[i]->UnRegister(this);
	}

	if (this->MuscleMesh != NULL)
	{
		this->MuscleMesh->Delete();
		this->MuscleMesh = NULL;
	}

	if (this->NormalMesh != NULL)
	{
		this->NormalMesh->Delete();
		this->NormalMesh = NULL;
	}

	tree->Delete();
	tree = NULL;
}

//----------------------------------------------------------------------------
// Adds single control fiber polyline data. Fiber data has to be sorted
// (from one end to another)and adding of whole fibers has to be in order from
// "left to right" also
void vtkMuscleDecomposer::AddFiberData(vtkPolyData *Data)
//----------------------------------------------------------------------------
{
	ControlFibers.push_back(Data);
	Data->Register(this);
}

//----------------------------------------------------------------------------
// Add Attachment data. If Bone data = NULL, then Data represents allready
// cut off attachment area surface. If Bone != NULL, then Data represent cut
// boundary polyline. That must be continuous, ordered and must lie on bone surface
void vtkMuscleDecomposer::AddAttachmentData(vtkPolyData *Data, vtkPolyData *Bone, int Position)
//----------------------------------------------------------------------------
{
	vtkPolyData* area = vtkPolyData::New();

	if (Bone == NULL) // if Bone is null, we expect allready cut off surface of the area
	{
		area->ShallowCopy(Data);
		AttachmentsSurfaces[Position].push_back(area);
	}
	else // If not, we have to cut it out
	{
		Log("  Extracting attachment area...");
		GetAttachmentAreaSurface(Data->GetPoints(), Bone, area);
		if (area->GetNumberOfPoints() > 0)
		{
			AttachmentsSurfaces[Position].push_back(area);
		}
	}
}

//----------------------------------------------------------------------------
//Creates the surface polydata that represents the attachment area on the muscle.
//Attachment area points (projpts) must lie on the muscle surface (surface), i.e., use ProjectAttachmentArea method prior to this one.
void vtkMuscleDecomposer::GetAttachmentAreaSurface(const vtkPoints* projpts, const vtkPolyData* surface, vtkPolyData* output)
//----------------------------------------------------------------------------
{
	// Adapted from vtkMuscleDecompositionAdv (ASM)
	//create cells for projpts
	int N = (int)const_cast<vtkPoints*>(projpts)->GetNumberOfPoints();
	vtkIdType ptIds[2] = { 0, 1 };

	vtkCellArray* cells = vtkCellArray::New();
	for (int i = 1; i < N; i++)
	{
		cells->InsertNextCell(2, ptIds);
		ptIds[0]++; ptIds[1]++;
	}

	ptIds[0] = N - 1;
	ptIds[1] = 0;
	cells->InsertNextCell(2, ptIds);

	//now we have here a contour
	vtkPolyData* contour = vtkPolyData::New();
	contour->SetPoints(const_cast<vtkPoints*>(projpts));
	contour->SetLines(cells);

	cells->Delete();

	//create cutter
	vtkMAFPolyDataCutOutFilterEx* cutter = vtkMAFPolyDataCutOutFilterEx::New();
#if VTK_MAJOR_VERSION < 5
	cutter->SetInput(const_cast<vtkPolyData*>(surface));
#else
	cutter->SetInputData(const_cast<vtkPolyData*>(surface));
#endif
	cutter->SetCuttingPolyline(contour);

	vtkPolyData* clippedPart = vtkPolyData::New();
	cutter->SetOutputClippedPart(clippedPart);

	//cutter->SetDebug((this->DebugMode & dbgVisualizeAttachmentConstruction) == dbgVisualizeAttachmentConstruction);
	cutter->Update();

	vtkPolyData* outputL = cutter->GetOutput();
	vtkPolyData* outputR = cutter->GetOutputClippedPart();

	//the smaller of both parts is the desired attachment area (surface)
	int numL = (int)outputL->GetNumberOfPoints();
	int numR = (int)outputR->GetNumberOfPoints();
	//BES 15.3.2012 - cutter may fail, if the input is corrupted, producing an empty set, which later causes problems
	//hence, here is a fix: we accept smaller part unless it is empty
	if (numL == 0)
		output->ShallowCopy(outputR);
	else if (numR == 0)
		output->ShallowCopy(outputL);
	else
		output->ShallowCopy((numL < numR ? outputL : outputR));

	contour->Delete();
	clippedPart->Delete();
	cutter->Delete();
	/*
	vtkSmartPointer< vtkSTLWriter > writer;
	writer->SetFileName("surface.stl");
	writer->SetInput(const_cast<vtkPolyData*>(surface));
	writer->Update();

	writer->SetFileName("attachment1.stl");
	writer->SetInput(output);
	writer->Update();
	*/
}

//----------------------------------------------------------------------------
// Adds mesh data for the muscle to be decomposed
void vtkMuscleDecomposer::AddMeshData(vtkPolyData *Mesh)
//----------------------------------------------------------------------------
{
	this->MuscleMesh = Mesh;	//TODO: remove Mesh and make it VTK consistent
	this->MuscleMesh->Register(this);

	this->SetInputData(0, Mesh);

	tree->SetDataSet(this->MuscleMesh);
	RebuildLocator = true; // The locator of the OBB tree will have to be updated
}

//// =========================================================================
//// Filter execution methods (the main entry points)
//// =========================================================================

//----------------------------------------------------------------------------
// Initializes all input data and performs data check (rough check, bigger problems
// arise as the decomposition goes along
bool vtkMuscleDecomposer::InitData()
//----------------------------------------------------------------------------
{
	Log("Initializing and checking input data...");
	OutputMesh = this->GetOutput(); // Link Output
	InitFibers();
	Log("  Number of fiber polylines: ", (int)FiberCount);
	if (FiberCount < 2)
	{
		Log("  Not enough fibers!");
		return false;
	}
	Log("  Number of vertices: ", (int)(PointCloud.size()));
	Log("  Fiber data seems OK.");

	if (MuscleMesh == NULL)
	{
		Log("  Muscle mesh not specified!");
		return false;
	}
	Log("  Muscle mesh seems OK.");

	if (AttachmentsSurfaces[0].size() == 0 || AttachmentsSurfaces[0].size() == 0)
	{
		Log("  Attachment data not specified. Muscle will not be connected to bones.");
	}
	else
	{
		Log("  Some attachment data found.");
	}

	Log("Data initialized OK.\n");
	return true;
}

//----------------------------------------------------------------------------
// Executes muscle decompositor
#if VTK_MAJOR_VERSION < 5
#define VTK6_RESULT_OK
#define VTK6_RESULT_FAIL
/**
Executes the data operation.

@param [in,out]	output	If non-null, the output.
*/
/*virtual*/ void vtkMuscleDecomposer::ExecuteData(vtkDataObject *output)
{
	//vtkPolyData* inputPoly = GetInput();
	//vtkPolyData* outputPoly = vtkPolyData::SafeDownCast(output);
#else

#define VTK6_RESULT_OK		1
#define VTK6_RESULT_FAIL	0

// This is the method of the algorithm in which the algorithm should fill in the output ports
/*virtual*/ int vtkMuscleDecomposer::RequestData(vtkInformation* request,
	vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
	//vtkPolyData* inputPoly = vtkPolyData::GetData(inputVector[0]);
	//vtkPolyData* outputPoly = vtkPolyData::GetData(outputVector);
#endif
	Log("Source executed.\n");
	if (!InitData()) {
		Log("Some data is incorrect. Stopping.");
		return VTK6_RESULT_FAIL;
	}
	BuildMesh();
	DoneMesh();
	ClearMesh();
	Log("Execution successful.");
	return VTK6_RESULT_OK;
}

//----------------------------------------------------------------------------
// Performs the decomposition and builds the output polylines
void vtkMuscleDecomposer::BuildMesh()
//----------------------------------------------------------------------------
{
	// It's called BuildMesh for nostalgic reasons from simpler times ;)
	bool check = true;

	times.push_back((long)MuscleMesh->GetNumberOfCells());
	start = clock();

	Log("Decomposing muscle...");
	int TargetTendonCount = (int)
		(SurfaceSubdivision > 0 ? (SurfaceSubdivision * (Fibers.size()) + 1) : Fibers.size()) //fibres in one layer
		*	(VolumeSubdivision + 1);  // volume fibers
	if (VolumeSubdivision > 0)
		TargetTendonCount += (VolumeSubdivision - 1) * 2; // side fibers

	PrepareTendons(TargetTendonCount);
	//RemoveArtifacts(ArtifactEpsilon, Fibers);
	times.push_back(clock() - start);

	// interpolate the top layer of fibers (input fibres) constructing SurfaceSubdivision fibres
	// between every pair of input fibres
	// we interpolate only top surface as the layers are interpolated automatically after sampling
	InterpolateSurface(SurfaceSubdivision, InterpolationMethod, Fibers);
	times.push_back(clock() - start);

	//estimate normals in the points of all the fibres
	ComputeNormals();
	times.push_back(clock() - start);
	//AlignMesh(); --- DOES NOT WORK

	//compute the thickness of the muscle in the points and spread the fibres if needed	
	ComputeThickness();
	times.push_back(clock() - start);
	Log("");
	
	//determine the scaling factor from the maximal thickness at fibres points
	double maxThickness = 0.0;

	std::vector<double> tScale;
	tScale.resize(Fibers.size());

	for (int i = 0; i < (int)Fibers.size(); i++)
	{
		//get the maximum thickness		
		for (int j = 0; j < (int)Fibers[i]->PointIDs.size(); j++) 
		{
			double thickness = PointCloud[Fibers[i]->PointIDs[j]]->Thickness;
			if (thickness > tScale[i]) {
				tScale[i] = thickness;
			}
		}

		if (tScale[i] > maxThickness) {
			maxThickness = tScale[i];
		}		
	}

	maxThickness = 0.8 * maxThickness;	//80% is allowed compression without removing some fibres
	for (int i = 0; i < (int)Fibers.size(); i++)
	{
		if (tScale[i] >= maxThickness)
			tScale[i] = 1.0;		//no multiplication
		else
			tScale[i] = maxThickness / tScale[i];
	}
	
	std::vector<Fiber*> layer;
	double step = 1.0f / VolumeSubdivision;
	// layer by layer, we reconstruct the volume
	for (double i = step; i < 1; i += step)
	{
		BuildBottomFibers(i, layer, step, tScale);
		if (DoOffsetFibers)
		{
			OffsetFibers(layer);
		}

		BottomFibers.insert(BottomFibers.end(), layer.begin(), layer.end());
		layer.clear();
	}

	if (VolumeSubdivision > 0)
	{
		for (int i = 0; i < (int)Fibers.size(); i++){			
			tScale[i] = 1.0;		//no multiplication
		}

		BuildBottomFibers(1, layer, 0, tScale);
		BottomFibers.insert(BottomFibers.end(), layer.begin(), layer.end());
		BuildSideLayers(VolumeSubdivision, InterpolationMethod, Fibers, layer);

		layer.clear();
	}
	times.push_back(clock() - start);

	if (this->GetDoConnectFibersToAA())
	{
		Log("\n  Connecting fibres...");
		ConnectAttachments();
		if (VolumeSubdivision > 0)
		{
			GatherAllFibers();
			Log("  Checking connection...");
			for (int i = 0; i < 10 && check; i++)
			{
				check = false;
				if (!CheckConnection(ATTACHMENT_START))
				{
					check = true;
					ShortenFibersByLayers(0.1, ATTACHMENT_START);
				}
				if (!CheckConnection(ATTACHMENT_END))
				{
					check = true;
					ShortenFibersByLayers(0.1, ATTACHMENT_END);
				}
				if (check)
				{
					Log("    Reconnecting...");
					if (i == 9)
					{
						Log("    Warning: This is too much reconnecting for the intput data to be\n    correct, I recon. Recheck the data please.");
					}
					ConnectAttachments();
				}
			}
		}

		times.push_back(clock() - start);
	}

	FiberCount = 0; // Now FiberCount represents the number of output fibers
	BuildAndOutputConnectedData(Fibers, "top fiber", true);
	BuildAndOutputConnectedData(BottomFibers, "inner layers fibers", true);
	BuildAndOutputConnectedData(LeftFibers, "left side fiber", true);
	BuildAndOutputConnectedData(RightFibers, "right side fiber", true);

	// And smoothening all fibres, we are done
	InterpolateFibers(InterpolationSubdivision, InterpolationMethod, OutputFibers, false);
	times.push_back(clock() - start);
	Log("Muscle built.\n");
	Log("Muscle fibers: ", FiberCount);
	LogTimes(times);
}

//----------------------------------------------------------------------------
// Extract data from internal structure and put them to output in vtkPolyData
// format
void vtkMuscleDecomposer::DoneMesh()
//----------------------------------------------------------------------------
{
	Log("Creating data on output...");
	vtkPoints * PointIDs = vtkPoints::New(); //delete in here
	vtkCellArray * lines = vtkCellArray::New(); //delete in here
	vtkUnsignedCharArray* mtInfo = vtkUnsignedCharArray::New();
	mtInfo->SetName("tendon_part");

	OutputFibersData(OutputFibers, PointIDs, lines, mtInfo);

	OutputMesh->SetPoints(PointIDs);
	OutputMesh->SetLines(lines);

	if (mtInfo->GetNumberOfTuples() != 0)
		OutputMesh->GetPointData()->AddArray(mtInfo);

	mtInfo->Delete();

	PointIDs->Delete();
	PointIDs = NULL;
	lines->Delete();
	lines = NULL;
#ifdef debug
	OutputNormals();
#endif

	Log("Mesh done.\n");
}

//----------------------------------------------------------------------------
// Clear data from memory
//----------------------------------------------------------------------------
void vtkMuscleDecomposer::ClearMesh()
{
	Log("Clearing data...");
	for (int i = 0; i < (int)Fibers.size(); i++)
	{
		Fibers[i]->PointIDs.clear();
		delete Fibers[i];
		Fibers[i] = NULL;
	}
	Fibers.clear();

	for (int i = 0; i < (int)OutputFibers.size(); i++)
	{
		OutputFibers[i]->PointIDs.clear();
		delete OutputFibers[i];
		OutputFibers[i] = NULL;
	}
	OutputFibers.clear();

	for (int i = 0; i < (int)BottomFibers.size(); i++)
	{
		BottomFibers[i]->PointIDs.clear();
		delete BottomFibers[i];
		BottomFibers[i] = NULL;
	}
	BottomFibers.clear();

	for (int i = 0; i < (int)LeftFibers.size(); i++)
	{
		LeftFibers[i]->PointIDs.clear();
		delete LeftFibers[i];
		LeftFibers[i] = NULL;
	}

	LeftFibers.clear();

	for (int i = 0; i < (int)RightFibers.size(); i++)
	{
		RightFibers[i]->PointIDs.clear();
		delete RightFibers[i];
		RightFibers[i] = NULL;
	}
	RightFibers.clear();

	for (int i = 0; i < (int)PointCloud.size(); i++)
	{
		if (PointCloud[i] != NULL)
		{
			delete PointCloud[i]; // everything inside PointCloud[i] is static...
			PointCloud[i] = NULL;
		}
	}
	PointCloud.clear();
	AllFibers.clear();
	for (int side = 0; side < 2; side++)
	{
		OpenFibers[side].clear();
		OpenFibersSeconds[side].clear();
		AttachmentPoints[side].clear();
	}

	Log("Data cleared.\n");
}

//// =========================================================================
//// Decomposition methods (called in BuildMesh())
//// =========================================================================

//----------------------------------------------------------------------------
// Converts fibers to internal represenation
void vtkMuscleDecomposer::InitFibers()
//----------------------------------------------------------------------------
{
	//make sure we have valid fibres
	std::vector<vtkPolyData*> fixedFibres;
	FixFibers(fixedFibres);
	
	vtkPolyData* Data;		// Data to be processed
	double pCoord[3];		// Point coords
	Point* pt;				// Point in internal structure

	int PointCount;			// Number of added points
	int PointID;			// ID of currently processed point
	int WhereIndex;			// Where will the new fiber go?

	FiberCount = 0;

	for (int i = 0; i < (int)fixedFibres.size(); i++)
	{
		Data = fixedFibres[i];
		PointCount = (int)Data->GetNumberOfPoints();
		PointCloud.reserve(PointCloud.size() + PointCount);

		//if cells are specified, the data is supposed to contain one fibre per one cell
		vtkIdType nLines = Data->GetNumberOfLines();
		Data->GetCellType(0);	//ensure the cells information is built

		vtkUnsignedCharArray* mtInfo =
			vtkUnsignedCharArray::SafeDownCast(Data->GetPointData()->GetArray("tendon_part"));

		int iFib = 0;
		do
		{
			WhereIndex = (int)Fibers.size();
			Fibers.push_back(new Fiber()); // delete in ClearMesh()

			vtkIdType nPts, *pPts;			
			Data->GetCellPoints(iFib, nPts, pPts);
			
			// For all points in Fiber, add them to cloud and to fiber
			for (int i = 0; i < nPts; i++)
			{
				vtkIdType origId = pPts == NULL ? i : pPts[i];

				Data->GetPoint(origId, pCoord); // get point from data
				pt = new Point(pCoord[0], pCoord[1], pCoord[2]); // delete in ClearMesh()

				PointCloud.push_back(pt);
				PointID = (int)PointCloud.size() - 1; // get its global id (in cloud)

				pt->Id = PointID;
				pt->MuscleTendonInfo = mtInfo != NULL ? mtInfo->GetValue(origId) : -1;

				Fibers[WhereIndex]->PointIDs.push_back(PointID); // add point id to specified fiber vector
			}
			FiberCount++;
		} while (++iFib < nLines);

		Data->Delete();
	}
}

//Fix the control fibres so that they contain the same number of segments.
//The fixed fibres are stored in fixedFibres (to be deleted by the caller)
void vtkMuscleDecomposer::FixFibers(std::vector<vtkPolyData*>& fixedFibres)
{
	vtkIdType nMaxSegments = 0;
	bool needFixing = false;

	for (int i = 0; i < (int)ControlFibers.size(); i++)
	{
		vtkPolyData* Data = ControlFibers[i];
		vtkIdType nLines = Data->GetNumberOfLines();

		//first make sure we have polylines
		if (nLines != 0){
			Data->Register(this);
		}
		else
		{
			vtkSmartPointer<vtkLHPPolyLinesBuilder> polyBuilder = vtkSmartPointer<vtkLHPPolyLinesBuilder>::New();
			polyBuilder->SetInputData(Data);
			polyBuilder->Update();

			Data = polyBuilder->GetOutput();
			Data->Register(this);

			nLines = Data->GetNumberOfLines();
		}

		fixedFibres.push_back(Data);

		//and now check all fibres in the Data
		Data->GetCellType(0);	//ensure the cells information is built

		int iFib = 0;
		do
		{
			vtkIdType nPts, *pPts;			
			Data->GetCellPoints(iFib, nPts, pPts);
			
			if (nPts != nMaxSegments)
			{
				needFixing = nMaxSegments != 0;

				if (nPts > nMaxSegments)
					nMaxSegments = nPts;
			}

		} while (++iFib < nLines);
	}

	if (needFixing)
	{
		for (int i = 0; i < (int)ControlFibers.size(); i++)
		{
			vtkPolyData* Data = fixedFibres[i];

			vtkSmartPointer<vtkLHPMuscleFibresResample> resampler = vtkSmartPointer<vtkLHPMuscleFibresResample>::New();
			resampler->SetInputData(Data);
			resampler->SetResolution((int)nMaxSegments);
			resampler->PrecisionModeOff();
			resampler->Update();

			Data->UnRegister(this);			
			Data = resampler->GetOutput();
			Data->Register(this);
			fixedFibres[i] = Data;
		}
	}
}

//UNCOMMENT THIS to get a good cover of large attachment sites
//perfect results for surface fibres but poor (unacceptable) for all other fibres

//#define HOTFIX_SURFACE

//----------------------------------------------------------------------------
// Connects fibers to tendons based on proximity
void vtkMuscleDecomposer::PrepareTendons(int targetFiberCount)
//----------------------------------------------------------------------------
{
	double *totalAreas, *minTriangleArea; // total and minimal triangle areas, for map building
	double **triangleAreas;			      // areas of all triangles
	double totalSideArea;                 // total area of the whole side
	int surfaces, triangles, areas, triangleID; // iterators
	double triangleArea;				  // temp variable

	int mapSize;                          // size of the map created
	int proportionalFiberCount;           // = targetFiberCount now
	int *map, *mapIterator;               // pointers to the map
	int areaIndicesCount;                 // number of indices in the map for single triangle

	vtkPolyData* data;                    // surface of the area
	int overflowCounter;				  // since map is integer numbers and sizes are double, we
	bool overflowDetected;                // need to check and correct if the map overflows due to
										  // rounding of numbers

	int mapSkip, sample;                  // how much we skip in map during sampling; sample ID

	double p0[3], p1[3], p2[3];		      // points of the triangle
	double baryCoords[3], result[3];      // barycentric coords; resulting attachment point coords

	Point* newPoint;					  // new attachment point as Point*
	int newID;							  // ID of our new poin

	vtkPolyData* test = vtkPolyData::New();  // for testing
	vtkPoints* testP = vtkPoints::New();   // for testing

	unsigned int rndseed = (unsigned int)time(NULL);
	rndseed = 1424763019;

	_RPT1(_CRT_WARN, "Random seed = %d\n", rndseed);

	srand(rndseed);					  // init the random generator
	Log("  Processing attachment areas...");

	//TODO: we need to generate a fine grid of attachment points, otherwise,
	//there will be always abrupt bends of fibres (the fibre should connect
	//optimally to point Q, however, argmin(i, ||Pi - Q||) > eps, where Pi
	//are the generated attachement points)
	//
	//TODO: BEWARE with more points, however, the angle part
	//of the cost function will be dominant and there will be places
	//without any fibre (it can be seen for iliacus with parameters
	//SS = 18, VS = 9, IS = 10) because a fibre could slightly bend
	//instead of going farer. Penalizing the angle would easily
	//lead to other problems and, therefore, it would be better
	//to measure curvature of a smooth curve at Pi instead of angle.
	//An ultimate solution (for one head) would be:
	//1) generate as many attachment points as reasonable
	//2) once final solution is reached, triangulate attachment points
	//3) perform relaxation to change the shape of triangles to
	//obtain the most equiangular triangles, i.e., this is in fact
	//smoothing of inner vertices with restriction that they lie on
	//the surface -> this could do the work

#ifdef HOTFIX_SURFACE
	const double SAMPLING_RES = 1.0;	//1 mm sampling resolution
#endif

	// we have 2 sides, muscle start and muscle end
	for (int side = 0; side < 2; side++)
	{
		surfaces = (int)AttachmentsSurfaces[side].size();
		totalAreas = new double[surfaces]; // delete[] here
		minTriangleArea = new double[surfaces]; // delete[] here
		triangleAreas = new double*[surfaces]; // delete[] here
		totalSideArea = 0;

		// each side may have multiple attachment areas
		for (int surface = 0; surface < surfaces; surface++)
		{
			data = AttachmentsSurfaces[side][surface];

			data->BuildCells();
			triangles = (int)data->GetNumberOfCells();
			minTriangleArea[surface] = DBL_MAX;
			totalAreas[surface] = 0;
			triangleAreas[surface] = new double[triangles]; // delete[] here

			// compute area of the surface and find minimum area (smallest triangle)
			for (int triangle = 0; triangle < triangles; triangle++)
			{
				vtkCell* cell = data->GetCell(triangle);
				vtkTriangle* t = dynamic_cast<vtkTriangle*>(cell);
				t->GetPoints()->GetPoint(0, p0);
				t->GetPoints()->GetPoint(1, p1);
				t->GetPoints()->GetPoint(2, p2);
				// triangle area
				triangleArea = vtkTriangle::TriangleArea(p0, p1, p2);
				triangleAreas[surface][triangle] = triangleArea;
				totalAreas[surface] += triangleArea;
				totalSideArea += triangleArea;
				minTriangleArea[surface] = std::min(triangleArea, minTriangleArea[surface]);
			}
		}

		// Now we have an area computed for all the surfaces on this side
		// we can build the map for all surfaces on this side
		for (int surface = 0; surface < surfaces; surface++)
		{
			areas = (int)AttachmentsSurfaces[side][surface]->GetNumberOfCells();

			proportionalFiberCount = targetFiberCount; // well, yeah.. it works better

			// Map size is determined by the smallest triangle, better more than less
			mapSize = (int)ceil(totalAreas[surface] / minTriangleArea[surface]);

#ifdef HOTFIX_SURFACE
			//Check, if the smallest triangle is not too big
			int resmapSize = (int)ceil(totalAreas[surface] / SAMPLING_RES);
			if (mapSize < resmapSize)
			{
				mapSize = resmapSize;
				minTriangleArea[surface] = totalAreas[surface] / resmapSize;
			}
#endif

			// If there are not enough triangles, we cannot really skip trough the map
			// so we adjust the minTriangleArea to be a smaller area than the actual smallest triangle
			if (mapSize < proportionalFiberCount)
			{
				minTriangleArea[surface] = totalAreas[surface] / targetFiberCount;
				mapSize = (int)ceil(totalAreas[surface] / minTriangleArea[surface]);
			}
			// Now we proceed by creating the map. Each map cell is an ID of a triangle
			// The bigger the triangle is, the more cells in map it has, the more probability
			// it will be picked for sampling a new attachment point
			map = mapIterator = new int[mapSize];
			overflowCounter = 0;
			overflowDetected = false;
			// map building
			for (int area = 0; area < areas && !overflowDetected; area++)
			{
				triangleArea = triangleAreas[surface][area];
				// this is the number of cells (indices) in the map for this triangle
				areaIndicesCount = vtkMath::Round(triangleArea / minTriangleArea[surface]);
				for (int counter = 0; counter < areaIndicesCount && !overflowDetected; counter++)
				{
					*mapIterator = area;
					mapIterator++;
					overflowCounter++;
					// overflow may happen, we rather exclude some triangles from the map
					// than face a SEGFAULT later
					overflowDetected = overflowCounter >= mapSize;
				}
			}

			// revert the pointer
			mapSize = (int)(mapIterator - map);
			delete[] triangleAreas[surface]; // cleanup, we dont need this anymore

#ifdef HOTFIX_SURFACE
			proportionalFiberCount = (int)ceil(totalSideArea / SAMPLING_RES);
#endif

			// how many cells in the map will we skip?
			mapSkip = (int)(mapSize / proportionalFiberCount);
			mapSkip = std::max(mapSkip, 1);

			data = AttachmentsSurfaces[side][surface];
			sample = 0;
			// Lines marked with ### are implementation of totally random sampling
			//int samples = 0; // ###
			while (/*samples < proportionalFiberCount ###*/sample < mapSize)
			{
				//sample = (rand() * rand()) % mapSize; // ###
				//samples++; // ###
				triangleID = map[sample];

				// we have a triangle picked
				vtkCell* cell = data->GetCell(triangleID);
				vtkTriangle* t = dynamic_cast<vtkTriangle*>(cell);
				t->GetPoints()->GetPoint(0, p0);
				t->GetPoints()->GetPoint(1, p1);
				t->GetPoints()->GetPoint(2, p2);
				// generate barycoords
				baryCoords[0] = (rand() % 10000 / (double)10000);
				baryCoords[1] = (rand() % 10000 / (double)10000);
				// Make sure we stay inside of the triangle
				if (baryCoords[0] + baryCoords[1] > 1)
				{
					baryCoords[0] = 1 - baryCoords[0];
					baryCoords[1] = 1 - baryCoords[1];
				}
				// barycoords must sum up to 1
				baryCoords[2] = 1 - baryCoords[0] - baryCoords[1];

				vtkMuscleDecomposerUtils::ApplyBarycentric(p0, p1, p2, baryCoords, result);
				// and we have ourselves a new attachment point
				newPoint = new Point(result[0], result[1], result[2]); //delete in ClearMesh()
				_CHECK_POINT(newPoint->DCoord);

				// insert it to cloud
				PointCloud.push_back(newPoint);
				newID = (int)PointCloud.size() - 1;
				newPoint->Id = newID;

				// TODO: For "points.vtk" debug output, keeping it, may be usefull
				testP->InsertNextPoint(result[0], result[1], result[2]);

				AttachmentPoints[side].push_back(newID);
				sample += mapSkip;
			}
			delete[] map; // we will recreate it for the other surfaces / sides
		}
		// there may be a different count of surface on the other side, we will have to reallocate these
		delete[] totalAreas;
		delete[] minTriangleArea;
		delete[] triangleAreas;

		cout << "Attachment Points (cummulative) = " << testP->GetNumberOfPoints() << endl;
	}

	// Debug. But i find it usefull at times.
	vtkPolyDataWriter* writer = vtkPolyDataWriter::New();
	char filename[100];
#if defined(_MSC_VER)
	sprintf_s
#else
	sprintf
#endif
		(filename, "points.vtk");

	writer->SetFileName(filename);
	test->SetPoints(testP);
#if VTK_MAJOR_VERSION < 5
	writer->SetInput(const_cast<vtkPolyData*>(test));
#else
	writer->SetInputData(const_cast<vtkPolyData*>(test));
#endif
	writer->Update();
	writer->Delete();

	testP->Delete();
	test->Delete();
}

//----------------------------------------------------------------------------
// Moves the muscle mesh so its centered in the fibre-defined muscle area
// This method is currently UNUSED as we no longer use the broken cadaver fibers
void vtkMuscleDecomposer::AlignMesh()
//----------------------------------------------------------------------------
{
	//return;
	Log("  Aligning muscle mesh to precise data...");
	double pCoord[3];		// Triangle vertices coords
	double Sub[3];
	double result[3];

	double MaxDistance = DBL_MIN;
	double MinDistance = DBL_MAX;
	int MaxDistanceIndex = 0;
	int MinDistanceIndex = 0;

	if (MuscleMesh == NULL)
		return;
	Log("    Computing direction...");
	int fiberPointCount = 0;

	for (int i = 0; i < (int)Fibers.size(); i++)
	{
		fiberPointCount += (int)Fibers[i]->PointIDs.size();
	}

	int* closestPoints = new int[fiberPointCount]; //delete in here
	double* distances = new double[fiberPointCount]; //delete in here
	std::vector<int> fiberPointIDS; //delete in here

	for (int i = 0; i < fiberPointCount; i++)
	{
		distances[i] = DBL_MAX;
	}

	fiberPointIDS.reserve(fiberPointCount);
	for (int j = 0; j < (int)Fibers.size(); j++)
	{
		for (int k = 0; k < (int)Fibers[j]->PointIDs.size(); k++)
		{
			fiberPointIDS.push_back(Fibers[j]->PointIDs[k]);
		}
	}

	int NumOfPoints = (int)MuscleMesh->GetNumberOfPoints();
	for (int i = 0; i < NumOfPoints; i++)
	{
		MuscleMesh->GetPoint(i, pCoord);
		for (int j = 0; j < fiberPointCount; j++)
		{
			double distance = vtkMuscleDecomposerUtils::PointPointDistance(pCoord, PointCloud[fiberPointIDS[j]]->DCoord);
			if (distance < distances[j])
			{
				distances[j] = distance;
				closestPoints[j] = i;
			}
		}
	}

	result[0] = result[1] = result[2] = 0;
	for (int i = 0; i < fiberPointCount; i++)
	{
		vtkMuscleDecomposerUtils::Subtract(PointCloud[fiberPointIDS[i]]->DCoord, MuscleMesh->GetPoint(closestPoints[i]), Sub);
		result[0] += Sub[0];
		result[1] += Sub[1];
		result[2] += Sub[2];
		if (distances[i] > MaxDistance)
		{
			MaxDistance = distances[i];
			MaxDistanceIndex = i;
		}
		if (distances[i] < MinDistance)
		{
			MinDistance = distances[i];
			MinDistanceIndex = i;
		}
	}
	result[0] /= fiberPointCount;
	result[1] /= fiberPointCount;
	result[2] /= fiberPointCount;

	vtkMuscleDecomposerUtils::Subtract(PointCloud[fiberPointIDS[MaxDistanceIndex]]->DCoord, MuscleMesh->GetPoint(closestPoints[MaxDistanceIndex]), result);

	for (int j = 0; j < fiberPointCount; j++)
	{
		PointCloud[fiberPointIDS[j]]->DCoord[0] += result[0];
		PointCloud[fiberPointIDS[j]]->DCoord[1] += result[1];
		PointCloud[fiberPointIDS[j]]->DCoord[2] += result[2];
	}

	delete[] distances;
	distances = NULL;
	fiberPointIDS.clear();
	delete[] closestPoints;
	closestPoints = NULL;
}

//----------------------------------------------------------------------------
// Computes normals for all fiber points from this neighborhood
//		|\2|6/| k
//		|1\|/5| k
//		|--*--| k
//		|3/|\7| k
//		|/4|8\|  jjjjjjj

void vtkMuscleDecomposer::ComputeNormals()
//----------------------------------------------------------------------------
{
	Log("  Computing normals...");
	std::vector<double*> normals;
	normals.reserve(8);
	double* normal;
	double* point1;
	double* point2;
	double* point3;

	// iterate all surface fibers
	for (int k = 0; k < (int)Fibers[0]->PointIDs.size(); k++)
	{
		for (int j = 0; j < (int)Fibers.size(); j++)
		{
			// this is the point where we are computing the normal
			point1 = PointCloud[Fibers[j]->PointIDs[k]]->DCoord;

			// go arount the neighbors and compute the normals
			if (j > 0) // do the left normals;
			{
				if (k < (int)Fibers[j - 1]->PointIDs.size() - 1)
				{
					normal = new double[3]; //delete in here
					point2 = PointCloud[Fibers[j - 1]->PointIDs[k]]->DCoord;
					point3 = PointCloud[Fibers[j - 1]->PointIDs[k + 1]]->DCoord;
					vtkMuscleDecomposerUtils::ComputeNormal(point1, point2, point3, normal);
					normals.push_back(normal); //1

					normal = new double[3]; //delete in here
					point2 = PointCloud[Fibers[j - 1]->PointIDs[k + 1]]->DCoord;
					point3 = PointCloud[Fibers[j]->PointIDs[k + 1]]->DCoord;
					vtkMuscleDecomposerUtils::ComputeNormal(point1, point2, point3, normal);
					normals.push_back(normal); //2
				}

				if (k > 0)
				{
					normal = new double[3];//delete in here
					point2 = PointCloud[Fibers[j - 1]->PointIDs[k - 1]]->DCoord;
					point3 = PointCloud[Fibers[j - 1]->PointIDs[k]]->DCoord;
					vtkMuscleDecomposerUtils::ComputeNormal(point1, point2, point3, normal);
					normals.push_back(normal); //3

					normal = new double[3];//delete in here
					point2 = PointCloud[Fibers[j]->PointIDs[k - 1]]->DCoord;
					point3 = PointCloud[Fibers[j - 1]->PointIDs[k - 1]]->DCoord;
					vtkMuscleDecomposerUtils::ComputeNormal(point1, point2, point3, normal);
					normals.push_back(normal); //4
				}
			}
			if (j < (int)Fibers.size() - 1) // do the right normals;
			{
				if (k < (int)Fibers[j + 1]->PointIDs.size() - 1)
				{
					normal = new double[3];//delete in here
					point2 = PointCloud[Fibers[j + 1]->PointIDs[k + 1]]->DCoord;
					point3 = PointCloud[Fibers[j + 1]->PointIDs[k]]->DCoord;
					vtkMuscleDecomposerUtils::ComputeNormal(point1, point2, point3, normal);
					normals.push_back(normal); //5

					normal = new double[3];//delete in here
					point2 = PointCloud[Fibers[j]->PointIDs[k + 1]]->DCoord;
					point3 = PointCloud[Fibers[j + 1]->PointIDs[k + 1]]->DCoord;
					vtkMuscleDecomposerUtils::ComputeNormal(point1, point2, point3, normal);
					normals.push_back(normal); //6
				}

				if (k > 0)
				{
					normal = new double[3];//delete in here
					point2 = PointCloud[Fibers[j + 1]->PointIDs[k]]->DCoord;
					point3 = PointCloud[Fibers[j + 1]->PointIDs[k - 1]]->DCoord;
					vtkMuscleDecomposerUtils::ComputeNormal(point1, point2, point3, normal);
					normals.push_back(normal); //7

					normal = new double[3];//delete in here
					point2 = PointCloud[Fibers[j + 1]->PointIDs[k - 1]]->DCoord;
					point3 = PointCloud[Fibers[j]->PointIDs[k - 1]]->DCoord;
					vtkMuscleDecomposerUtils::ComputeNormal(point1, point2, point3, normal);
					normals.push_back(normal); //8
				}
			}
		} // j
		// average normals and assign as normal
		// the normal is same for the whole tier of points (to avoid normal intersections)
		normal = new double[3]; //delete in here
		normal[0] = normal[1] = normal[2] = 0;
		for (int i = 0; i < (int)normals.size(); i++)
		{
			normal[0] += normals[i][0];
			normal[1] += normals[i][1];
			normal[2] += normals[i][2];
			delete[] normals[i];
			normals[i] = NULL;
		}
		vtkMath::Normalize(normal);
		for (int l = 0; l < (int)Fibers.size(); l++)
		{
			PointCloud[Fibers[l]->PointIDs[k]]->SetNormal(normal);
		}
		normals.clear();
		delete normal;
		normal = NULL;
	} // k
}

//----------------------------------------------------------------------------
// Multiplies all normals by -1
void vtkMuscleDecomposer::FlipNormals()
//----------------------------------------------------------------------------
{
	if (DoFlipNormals)
	{
		for (int j = 0; j < (int)Fibers[0]->PointIDs.size(); j++) // NEEDS SAME FIBER POINT COUNTS
		{
			for (int i = 0; i < (int)Fibers.size(); i++) // row by row
			{
				PointCloud[Fibers[i]->PointIDs[j]]->DNormal[0] *= -1;
				PointCloud[Fibers[i]->PointIDs[j]]->DNormal[1] *= -1;
				PointCloud[Fibers[i]->PointIDs[j]]->DNormal[2] *= -1;
			}
		}
	}
}

//----------------------------------------------------------------------------
// Outputs normals - debug purposes only.
//----------------------------------------------------------------------------
void vtkMuscleDecomposer::OutputNormals()
{
	vtkPoints * PointIDs = vtkPoints::New(); //delete in here
	vtkCellArray * lines = vtkCellArray::New(); //delete in here
	vtkIdType * lineIDs = new vtkIdType[2]; //delete in here
	double moved[3];
	double* original;
	double* normal;

	for (int i = 0; i < (int)Fibers.size(); i++)
	{
		for (int j = 0; j < (int)Fibers[i]->PointIDs.size(); j++)
		{
			// insert point to output
			original = PointCloud[Fibers[i]->PointIDs[j]]->DCoord;
			lineIDs[0] = PointIDs->InsertNextPoint(original);
			normal = PointCloud[Fibers[i]->PointIDs[j]]->DNormal;
			vtkMuscleDecomposerUtils::Add(original, normal, moved, -PointCloud[Fibers[i]->PointIDs[j]]->Thickness);
			lineIDs[1] = PointIDs->InsertNextPoint(moved);
			lines->InsertNextCell(2, lineIDs);
		}
	}

	NormalMesh = vtkPolyData::New(); //not deleted (used for debugging only)
	NormalMesh->SetPoints(PointIDs);
	NormalMesh->SetLines(lines);

	delete[] lineIDs;
	lineIDs = NULL;
	PointIDs->Delete();
	lines->Delete();
}

//----------------------------------------------------------------------------
// Computes muscle thinckness for fiber points, sampling the muscle volume
void vtkMuscleDecomposer::ComputeThickness(void)
//----------------------------------------------------------------------------
{
	Log("  Computing muscle thickness...");
	if (RebuildLocator)
	{
		tree->BuildLocator();
		RebuildLocator = false;
	}

	// Now we test for normal orientation
	Log("    Testing normal orientation...");

	int NumOfTests = 16;
	int successfull = 0;
	int fibre, cut;
	int ID;
	double testPoint[3];
	int test;
	for (int i = 0; i < NumOfTests; i++)
	{
		// we pick some random fibers
		fibre = rand() % Fibers.size();
		cut = rand() % Fibers[0]->PointIDs.size();
		ID = Fibers[fibre]->PointIDs[cut];
		// create a test point in the middle of the sample
		vtkMuscleDecomposerUtils::Add(PointCloud[ID]->DCoord, PointCloud[ID]->DNormal, testPoint, 
			-SampleThickness(PointCloud[ID]->DCoord, PointCloud[ID]->DNormal) / 2);
		// and we test if it is inside of the model
		test = tree->InsideOrOutside(testPoint);
		if (test == -1)
		{
			successfull++;
		}
	}

	if (successfull >= NumOfTests / 8 * 7 || DoForceFlipNormals)
	{
		Log("      Normal direction seems OK.");
	}
	else
	{
		Log("      Incorrect orientation detected, normals will be flipped.");
		FlipNormals();
	}		

	times.push_back(clock() - start);
	Log("    Sampling muscle thickness from mesh data...");
	// First we compute muscle thicknes from mesh data where possible

	int nFibers = (int)Fibers.size();

#pragma omp parallel for
	for (int i = 0; i < nFibers; i++)
	{
		for (int j = 0; j < (int)Fibers[i]->PointIDs.size(); j++)
		{
			//BES: 29.5.2017 - if the points do not lie on the surface, the error is too big to ignore => sample in one side only
			double thickness = SampleThickness(PointCloud[Fibers[i]->PointIDs[j]]->DCoord, PointCloud[Fibers[i]->PointIDs[j]]->DNormal, true);
			PointCloud[Fibers[i]->PointIDs[j]]->Thickness = thickness;

			//_RPT2(_CRT_WARN, "#%d.%d = %f\n", i, j, thickness);
		}
	}

	times.push_back(clock() - start);

	//TODO: if the muscle is non-convex, e.g., U shape, there is always some thickness =>
	//hole cannot be identified just by whole muscle thickness but we need to use intervals:
	//space - muscle - space - muscle - ... space and determine, if the distance of the point
	//towards the nearest muscle interval is < THRESHOLD, if not => HOLE

	//TODO: Search For Valid muscle thickness must be changed accordingly

	//TODO: spreading of fibres must be done for each layer, not only superficial one

	if (this->GetDoSpreadFibers())
	{
		Log("    Spreading fibres...");
		// Now we check if the muscle is split somewhere
		bool hasHole;
		double holeT;
		double validTs[2];
		int LeftAnchorID, RightAnchorID;
		Point* dummy;

		std::vector<int> segmentVertices;
		double step;
		double minT;
		double maxT;
		int count;
		double newCoord[3];

		//surface is already interpolated so we have in Fibers data for fibres of two types, A and B,
		//in this fashion: A B B ... B A B ...B A B ... B A
		//A is a fibre defined by the input points of a fibre (P)
		//B is a fibre constructed by interpolating A fibres
		//all the fibres have the same number of points
		//we assume that A fibres lie correctly above the surface, i.e., muscle thickness is > 0
		//B.LeftAnchorId denotes the nearest left A fibre, B.RightAnchorId denotes the nearest right A fibre
		//A.LeftAnchorId = B.LeftAnchodId = -1
		for (int j = 0; j < (int)Fibers[0]->PointIDs.size(); j++) // NEEDS SAME FIBER POINT COUNTS
		{
			hasHole = false;		//BES: 24.2.2015 - fixed bug
			LeftAnchorID = -1;	//we always start from A fibre
			RightAnchorID = -1;

			for (int i = 0; i < (int)Fibers.size(); i++) // row by row
			{
				// now we iterate every row and find hits and misses
				// we test every segment for holes
				dummy = PointCloud[Fibers[i]->PointIDs[j]];

				if (dummy->LeftAnchorID != -1) // it's an intermediate point
				{
					LeftAnchorID = dummy->LeftAnchorID;
					RightAnchorID = dummy->RightAnchorID;
					segmentVertices.push_back(dummy->Id);

					if (dummy->Thickness < MIN_THICKNESS)
					{
						hasHole = true;
						holeT = dummy->t;
					}
				}

				if (dummy->Id == RightAnchorID)
				{
					// we tested the segment, now its time to processed
					if (hasHole)
					{
						// now we search for valid muscle on the lefthand and on the righthand side of the segment
						SearchValidThickness(LeftAnchorID, RightAnchorID, holeT, validTs, BINARY_SEARCH_ITERATIONS);
						// now we resample the T and move the vertices in both groups

						// first, the left side vertices
						minT = PointCloud[LeftAnchorID]->t;
						maxT = validTs[0];
						count = (int)segmentVertices.size() >> 1;
						step = (maxT - minT) / count;
						for (int k = 0; k < count; k++)
						{
							minT += step;
							PointCloud[LeftAnchorID]->ComputeInterpolatedCurveCoord(minT, newCoord, InterpolationMethod != SPLINE_CONSTRAINED && InterpolationMethod != SPLINE_CATMULL_ROM);
							// now we move the point
							PointCloud[segmentVertices[k]]->DCoord[0] = newCoord[0];
							PointCloud[segmentVertices[k]]->DCoord[1] = newCoord[1];
							PointCloud[segmentVertices[k]]->DCoord[2] = newCoord[2];
							PointCloud[segmentVertices[k]]->t = minT;
							double thickness = SampleThickness(newCoord, PointCloud[segmentVertices[k]]->DNormal, true);
							PointCloud[segmentVertices[k]]->Thickness = thickness;
							if (thickness < MIN_THICKNESS)
							{
								Log("      Muscle spreading error at fibre: ", i);
								Log("      Segment: ", k);
								Log("      Control fibres outside of the surface?");
							}
						}
						PointCloud[segmentVertices[count - 1]]->Splits = true;

						// now, the right side vertices
						minT = validTs[1];
						maxT = PointCloud[RightAnchorID]->t;
						count = (int)segmentVertices.size() - count;
						step = (maxT - minT) / count;
						for (int k = count; k < (int)segmentVertices.size(); k++)
						{
							PointCloud[LeftAnchorID]->ComputeInterpolatedCurveCoord(minT, newCoord, InterpolationMethod != SPLINE_CONSTRAINED && InterpolationMethod != SPLINE_CATMULL_ROM);
							// now we move the point
							PointCloud[segmentVertices[k]]->DCoord[0] = newCoord[0];
							PointCloud[segmentVertices[k]]->DCoord[1] = newCoord[1];
							PointCloud[segmentVertices[k]]->DCoord[2] = newCoord[2];
							PointCloud[segmentVertices[k]]->t = minT;
							double thickness = SampleThickness(newCoord, PointCloud[segmentVertices[k]]->DNormal, true);
							PointCloud[segmentVertices[k]]->Thickness = thickness;
							if (thickness < MIN_THICKNESS)
							{
								Log("      Muscle spreading error at fibre: ", i);
								Log("      Segment: ", k);
								Log("      Control fibres outside of the surface?");
							}
							minT += step;
						}
					}
					hasHole = false;
					segmentVertices.clear();
				}
			}
		}
	}

	times.push_back(clock() - start);	
}

//----------------------------------------------------------------------------
// Computes thicknes at specified coord along the specified normal
double vtkMuscleDecomposer::SampleThickness(double* coord, double* normal, bool oneSideOnly)
//----------------------------------------------------------------------------
{
	//TODO: there are several scenarios
	//assuming a convex muscle, the specified point may lie above the surface
	//below the surface or exactly on the surface


	double rayPoint1[3];
	double rayPoint2[3];

	double intersectionPoint1[3];
	double intersectionPoint2[3];

	// create the ray
	vtkMuscleDecomposerUtils::Add(coord, normal, rayPoint1, (oneSideOnly ? 0 : RAY_POINT_DISTANCE));
	vtkMuscleDecomposerUtils::Add(coord, normal, rayPoint2, -RAY_POINT_DISTANCE);

	// search for intersections
	vtkPoints* intersections = vtkPoints::New(); //delete in here	
	tree->IntersectWithLine(rayPoint1, rayPoint2, intersections, NULL);

	double thickness = 0;
	if (intersections->GetNumberOfPoints() > (oneSideOnly ? 0 : 1))
	{
		// thickness is the distance between the two intersecting points
		if (oneSideOnly && intersections->GetNumberOfPoints() % 2 != 0)
		{
			//if coord is inside, we need to add the distance between coord and the first intersection
			intersections->GetPoint(0, intersectionPoint1);
			thickness += vtkMuscleDecomposerUtils::PointPointDistance(coord, intersectionPoint1);
		}

		for (int i = 0; i < intersections->GetNumberOfPoints() - 1; i += 2)
		{
			intersections->GetPoint(i, intersectionPoint1);
			intersections->GetPoint(i + 1, intersectionPoint2);
			// thickness is the distance between the two intersecting points
			thickness += vtkMuscleDecomposerUtils::PointPointDistance(intersectionPoint1, intersectionPoint2);
		}
	}

	intersections->Delete();
	intersections = NULL;

	return thickness;
}

//----------------------------------------------------------------------------
// Searches along vertex slice for 0.8 times the anchor thickness. Searches
// both sides
void vtkMuscleDecomposer::SearchValidThickness(int LeftAnchorID, int RightAnchorID, double holeT, double* result, int iterations)
//----------------------------------------------------------------------------
{
	double weightCoef = 0.8;

	int iteration;
	double leftT, rightT, currentT;
	double thickness;
	double currentCoord[3];

	// left part
	leftT = PointCloud[LeftAnchorID]->t;
	rightT = holeT;
	iteration = 0;
	while (iteration++ < iterations)
	{
		currentT = (leftT + rightT) / 2;
		// the values of a, b, c and d are already at the anchor point, from previous linear interpolation!
		PointCloud[LeftAnchorID]->ComputeInterpolatedCurveCoord(currentT, currentCoord, InterpolationMethod != SPLINE_CONSTRAINED && InterpolationMethod != SPLINE_CATMULL_ROM);
		// normals are all equal on the row, so we may use the anchor point's normal safely
		thickness = SampleThickness(currentCoord, PointCloud[LeftAnchorID]->DNormal);

		if (thickness > weightCoef * PointCloud[LeftAnchorID]->Thickness) // binary search
		{
			result[0] = currentT; // we found the thickness is valid
			leftT = currentT;
		}
		else
		{
			rightT = currentT;
		}
	}

	// right part
	leftT = holeT;
	rightT = PointCloud[RightAnchorID]->t;
	iteration = 0;
	while (iteration++ < iterations)
	{
		currentT = (leftT + rightT) / 2;
		// the values of a, b, c and d are already at the anchor point, from previous linear interpolation!
		// segments are interpolated from left to right
		PointCloud[LeftAnchorID]->ComputeInterpolatedCurveCoord(currentT, currentCoord, InterpolationMethod != SPLINE_CONSTRAINED && InterpolationMethod != SPLINE_CATMULL_ROM);
		// normals are all equal on the row, so we may use the anchor point's normal safely
		thickness = SampleThickness(currentCoord, PointCloud[RightAnchorID]->DNormal);

		if (thickness > weightCoef * PointCloud[RightAnchorID]->Thickness) // binary search
		{
			result[1] = currentT; // only if the found thickness is valid
			rightT = currentT;
		}
		else
		{
			leftT = currentT;
		}
	}
}

//----------------------------------------------------------------------------
// Builds bottom layer of fibers from normals and thickness
void vtkMuscleDecomposer::BuildBottomFibers(double t, std::vector<Fiber*>  & target, 
	double step, std::vector<double>& tScalePerFibre)
//----------------------------------------------------------------------------
{
	//TODO: non-convex muscles are not supported! Although thickness sums only intervals inside the muscle
	//here we assume just one interval
	//TODO: furthermore, if the points does not lie precisely on the surface, there is some shift since
	//we do not sample the interval but interval <P, P + thickness>

	if (t < 1)
		Log("  Building a layer of inner fibers...");
	else
		Log("  Building bottom fibers...");

	int ID = 0;
	int newID = 0;
	Fiber* newFiber;
	Point* newPoint;
	// copy all fibres
	target.reserve(Fibers.size());
	double pCoord[3];
	double offset = 0;

	for (int i = 0; i < (int)Fibers.size(); i++)
	{
		double tCorr = t * tScalePerFibre[i];
		if (tCorr > 1.0) {
			continue;	//skip the fibre
		}

		newFiber = new Fiber(); //delete in ClearMesh()
		newFiber->DLayer = t;
		newFiber->PointIDs.reserve(Fibers[i]->PointIDs.size());
		if (DoOffsetFibers) // suppress the top-down
		{
			offset = -step / 2 + (rand() % 1000) / 1000.0f * step;
		}
		for (int j = 0; j < (int)Fibers[i]->PointIDs.size(); j++)
		{
			ID = Fibers[i]->PointIDs[j];
			// compute point position from normal and thickness

			vtkMuscleDecomposerUtils::Add(PointCloud[ID]->DCoord, PointCloud[ID]->DNormal, pCoord, -PointCloud[ID]->Thickness * (tCorr + offset));
			_CHECK_POINT(pCoord);

			// create point
			newPoint = new Point(pCoord[0], pCoord[1], pCoord[2]); //delete in ClearMesh()
			newPoint->Splits = PointCloud[ID]->Splits;

			// insert it to cloud
			PointCloud.push_back(newPoint);
			newID = (int)PointCloud.size() - 1;
			newPoint->Id = newID;
			// adjust position
			newPoint->AttachmentConnectionID = PointCloud[ID]->AttachmentConnectionID;
			newPoint->Interpolated = PointCloud[ID]->Interpolated;
			// and to the list of IDs for this fibre
			newFiber->PointIDs.push_back(newID);
		}
		target.push_back(newFiber);
		RegisterOpenEnds(newFiber);
	}
}

//----------------------------------------------------------------------------
// Offsets the fibers inside the layer of fibers specified to
// suppress the left-right "rake" effect
void vtkMuscleDecomposer::OffsetFibers(std::vector<Fiber*>  & target)
//----------------------------------------------------------------------------
{
	//TODO: This shifts points of one fibre differently. Is it OK?
	//TODO: What is not OK is that it shifts the points iteratively and in one direction only,
	//so it may happen that a straight fibre is no longer straight, example:
	//Fi-1 = [(0,0,0),(2,3,4),(3,4,8)]
	//Fi =	 [(5,2,0),(5,4,5),(5,5,10)]
	//Fi+1 = [(10,2,0),(10,4,5),(10,5,10)]
	//Fi+2 = [(15,2,0),(15,4,5),(15,5,10)]
	//-> with random = 0.5 for all cases, iterative version results in:
	//Fi-1' = Fi
	//Fi' =		[(2.5,1,0),(3.5,3.5,4.5),(4,4.5,9)]
	//Fi+1' = [(6.25,1.5,0),(6.75,3.75,5),(7.5,5,10)]
	//Fi+2' = [(10.625,1.75,0),(10.875, ...]
	//although originally Fi to Fi+2 was flat (x = const) and all fibres parallel (delta x = const), suddenly they are not
	//with non-iterative version
	//Fi' =		[(2.5,1,0),(3.5,3.5,4.5),(4,4.5,9)]
	//Fi+1' = [(7.5,2,0),(7.5,4,5),(7.5,5,10)]
	//Fi+2' = [(12.5,2,0),(12.5,4,5),(12.5,5,10)]
	//=> Fi+1 and Fi+2 are straight and parallel, Fi is influenced by Fi-1

	Fiber *first, *second;
	double offset;
	double vector[3];
	double length;

	for (int i = 0; i < (int)target.size() - 1; i++)
	{
		first = target[i];
		second = target[i + 1];

		offset = (rand() % 1000) / 1000.0f;

		for (int j = 0; j < (int)first->PointIDs.size(); j++)
		{
			vtkMuscleDecomposerUtils::Subtract(PointCloud[first->PointIDs[j]]->DCoord, PointCloud[second->PointIDs[j]]->DCoord, vector);
			length = vtkMath::Norm(vector);
			if (length > 0.1 && !PointCloud[first->PointIDs[j]]->Splits)
			{
				//result[i] = point2[i] * multiplication + point1[i];
				vtkMuscleDecomposerUtils::Add(PointCloud[first->PointIDs[j]]->DCoord, vector, PointCloud[first->PointIDs[j]]->DCoord, offset);

				_CHECK_POINT(PointCloud[first->PointIDs[j]]->DCoord);
			}
		}
	}
}

//----------------------------------------------------------------------------
// Builds side layers (use after bottom layer had been built)
void vtkMuscleDecomposer::BuildSideLayers(int subdivision, int interpolationType, std::vector<Fiber*>  & topFibers, std::vector<Fiber*>  & bottomFibers)
//----------------------------------------------------------------------------
{
	LogFormated("  Interpolating muscle side points using %s method...", GetSplineType(interpolationType));
	if (topFibers.size() < 2 || bottomFibers.size() < 2) {
		Log("	Not enough side fibres for the interpolation.");
		return;
	}

	Fiber* crossPoints;
	Fiber* newFiber;
	int fiberCount = subdivision;
	int pointCount = (int)topFibers[0]->PointIDs.size();

	LeftFibers.reserve(fiberCount);
	RightFibers.reserve(fiberCount);

	SplineInterpolator* interpolator = new SplineInterpolator(interpolationType); //delete in here

	//we build both sides simultaniously
	for (int j = 0; j < pointCount; j++)
	{
		// left side
		crossPoints = new Fiber(); //delete in here
		crossPoints->PointIDs.reserve(fiberCount);
		for (int i = 1; i >= 0; i--) // pick 2 border fibres
		{
			crossPoints->PointIDs.push_back(topFibers[i]->PointIDs[j]);
		}

		for (int i = 0; i < 2; i++) // pick 2 border fibres
		{
			crossPoints->PointIDs.push_back(bottomFibers[i]->PointIDs[j]);
		}

		InterpolateFiber(crossPoints, subdivision, interpolator, 1, true);

		for (int i = 2; i < (int)crossPoints->PointIDs.size() - 2; i++) // we do not want to have the duplicates
		{
			if ((int)LeftFibers.size() <= i - 2) // add fiber if not there
			{
				newFiber = new Fiber();
				newFiber->DLayer = (double)i / crossPoints->PointIDs.size();
				LeftFibers.push_back(newFiber); //delete in ClearMesh()
			}
			LeftFibers[i - 2]->PointIDs.push_back(crossPoints->PointIDs[i]); // add vertex slice
		}

		crossPoints->PointIDs.clear();
		// right side
		for (int i = (int)topFibers.size() - 2; i < (int)topFibers.size(); i++) // pick 2 border fibres
		{
			crossPoints->PointIDs.push_back(topFibers[i]->PointIDs[j]);
		}

		for (int i = (int)bottomFibers.size() - 1; i >= (int)bottomFibers.size() - 2; i--) // pick 2 border fibres
		{
			crossPoints->PointIDs.push_back(bottomFibers[i]->PointIDs[j]);
		}

		InterpolateFiber(crossPoints, subdivision, interpolator, 1, true);

		for (int i = 2; i < (int)crossPoints->PointIDs.size() - 2; i++) // we do not want to have the duplicates
		{
			if ((int)RightFibers.size() <= i - 2)
			{
				newFiber = new Fiber();
				newFiber->DLayer = (double)i / crossPoints->PointIDs.size();
				RightFibers.push_back(newFiber); //delete in ClearMesh()
			}
			RightFibers[i - 2]->PointIDs.push_back(crossPoints->PointIDs[i]); // add vertex slice
		}

		delete crossPoints;
		crossPoints = NULL;
	}

	// these fibers need to be registered too, so they will be connected as well
	for (int i = 0; i < (int)LeftFibers.size(); i++)
	{
		RegisterOpenEnds(LeftFibers[i]);
	}

	for (int i = 0; i < (int)RightFibers.size(); i++)
	{
		RegisterOpenEnds(RightFibers[i]);
	}
	delete interpolator;
	interpolator = NULL;
}

//----------------------------------------------------------------------------
// Connects attachment points to fibers, using "dist(p1,b) * angle(p1, p2, b)"
// metric and Hungarian method for the Aassignment problem.
// (b) ---------- (p1) ---------- (p2) ----------- ... rest of the fibre
void vtkMuscleDecomposer::ConnectAttachments()
//----------------------------------------------------------------------------
{
	double* point1; // end point of the fibre
	double* point2; // first inner point of the fibre

	double distance, angle; // self explanatory

	int** distances; // a matrix of metrics for the Hungarian method
	bool reverseFibers = false; // debug
	int sideRev; // debug

	hungarian_problem_t prob;

	if (AttachmentPoints[0].size() < 1 || AttachmentPoints[1].size() < 1)
	{
		Log("  ERROR: Not enough tendon points. Some fibers may be connected incorrectly!");
		return;
	}

	for (int side = 0; side < 2; side++)
	{
		if (side == 0 == reverseFibers) // should we reverse the fibers?
		{
			sideRev = 0;
		}
		else
		{
			sideRev = 1;
		}
		// allocate, allocate
		distances = new int*[AttachmentPoints[side].size()];
		for (int i = 0; i < (int)AttachmentPoints[side].size(); i++)
		{
			distances[i] = (int*) new int[OpenFibers[sideRev].size()];
		}

		// build the metrics matrix
#pragma omp parallel for
		for (int i = 0; i < (int)OpenFibers[sideRev].size(); i++)
		{
			for (int j = 0; j < (int)AttachmentPoints[side].size(); j++)
			{
				point1 = PointCloud[OpenFibers[sideRev][i]]->DCoord;
				point2 = PointCloud[OpenFibersSeconds[sideRev][i]]->DCoord;

				distance = vtkMuscleDecomposerUtils::PointPointDistance(point1, PointCloud[AttachmentPoints[side][j]]->DCoord);

				angle = vtkMuscleDecomposerUtils::ComputeBendAngle(point2, point1, PointCloud[AttachmentPoints[side][j]]->DCoord);
				angle /= vtkMath::Pi();

				// this clamping is just cosmetic, but improves results for higher amounts of fibers
				//angle = angle < 0.1 ? 0.1 : angle;

				distance *= angle;

				// the library is capable to work only with integers, sadly
				distances[j][i] = (int)(distance * 100000);
			}
		}

		// let the libhungarian do its magic :)
		hungarian_init(&prob, distances, (int)AttachmentPoints[side].size(), (int)OpenFibers[sideRev].size(), HUNGARIAN_MODE_MINIMIZE_COST);
		hungarian_solve(&prob);
		// taking the solution, we connect the fibres
		for (int i = 0; i < (int)AttachmentPoints[side].size(); i++)
		{
			for (int j = 0; j < (int)OpenFibers[sideRev].size(); j++)
			{
				if (prob.assignment[i][j] == HUNGARIAN_ASSIGNED)
				{
					PointCloud[OpenFibers[sideRev][j]]->AttachmentConnectionID = AttachmentPoints[side][i];
					break;
				}
			}
		}

		// cleanup the room
		hungarian_free(&prob);
		for (int j = 0; j < (int)AttachmentPoints[side].size(); j++)
		{
			delete[] distances[j];
		}
		delete[] distances;
	}
}

//----------------------------------------------------------------------------
// Adds IDs of endpoints and next-to-endpoints to appropriate open fibre lists
void vtkMuscleDecomposer::RegisterOpenEnds(Fiber* fiber)
//----------------------------------------------------------------------------
{
	OpenFibers[0].push_back(fiber->PointIDs[0]);
	OpenFibersSeconds[0].push_back(fiber->PointIDs[1]);
	OpenFibers[1].push_back(fiber->PointIDs[fiber->PointIDs.size() - 1]);
	OpenFibersSeconds[1].push_back(fiber->PointIDs[fiber->PointIDs.size() - 2]);
}

//----------------------------------------------------------------------------
// Removes IDs of endpoints and next-to-endpoints from appropriate open fibre lists
void vtkMuscleDecomposer::UnRegisterOpenEnds(Fiber* fiber)
//----------------------------------------------------------------------------
{
	OpenFibers[0].erase(std::remove(OpenFibers[0].begin(), OpenFibers[0].end(), fiber->PointIDs[0]), OpenFibers[0].end());
	OpenFibersSeconds[0].erase(std::remove(OpenFibersSeconds[0].begin(), OpenFibersSeconds[0].end(), fiber->PointIDs[1]), OpenFibersSeconds[0].end());

	OpenFibers[1].erase(std::remove(OpenFibers[1].begin(), OpenFibers[1].end(), fiber->PointIDs[fiber->PointIDs.size() - 1]), OpenFibers[1].end());
	OpenFibersSeconds[1].erase(std::remove(OpenFibersSeconds[1].begin(), OpenFibersSeconds[1].end(), fiber->PointIDs[fiber->PointIDs.size() - 2]), OpenFibersSeconds[1].end());
}

//----------------------------------------------------------------------------
// Gathers all fibers to AllFibers list
void vtkMuscleDecomposer::GatherAllFibers()
//----------------------------------------------------------------------------
{
	AllFibers.insert(AllFibers.end(), Fibers.begin(), Fibers.end());
	AllFibers.insert(AllFibers.end(), BottomFibers.begin(), BottomFibers.end());
	AllFibers.insert(AllFibers.end(), LeftFibers.begin(), LeftFibers.end());
	AllFibers.insert(AllFibers.end(), RightFibers.begin(), RightFibers.end());
}

//----------------------------------------------------------------------------
// Checks the connection of all fibres for sharp connection angle
bool vtkMuscleDecomposer::CheckConnection(int side)
//----------------------------------------------------------------------------
{
	bool correct = true;

	double* point1; // endpoint
	double* point2; // next-to-endpoint
	double* point3; // tendon attachment point

	int ID1, ID2; // IDs of point1 and point2
	int attachmentID; // ID of point3

	double angle;	// the angle to be checked

	Fiber* fiber; // the processed fiber

	for (int i = 0; i < (int)AllFibers.size(); i++)
	{
		fiber = AllFibers[i];
		// which side do we check?
		if (side == 0)
		{
			ID1 = fiber->PointIDs[0];
			ID2 = fiber->PointIDs[1];
		}
		else
		{
			ID1 = fiber->PointIDs[fiber->PointIDs.size() - 1];
			ID2 = fiber->PointIDs[fiber->PointIDs.size() - 2];
		}

		point1 = PointCloud[ID1]->DCoord;
		point2 = PointCloud[ID2]->DCoord;
		attachmentID = PointCloud[ID1]->AttachmentConnectionID;
		if (attachmentID >= 0)
		{
			point3 = PointCloud[attachmentID]->DCoord;
			// here is the angle
			angle = vtkMuscleDecomposerUtils::ComputeBendAngle(point2, point1, point3);
			// and here is the treshold
			if (angle > PennationTreshold)
			{
				// if the angle is too sharp, the test fails
				return false;
			}
		}
	}
	return true;
}

//----------------------------------------------------------------------------
// Shortens fibers based on the layer of fibers they are in, on one side.
// amount is the maximum amount of shorteninfg, in percent. Only bottom
// fibers will be shortened by this, the others by less, based on how
// deep in the muscle they are
void vtkMuscleDecomposer::ShortenFibersByLayers(double amount, int side)
//----------------------------------------------------------------------------
{
	Fiber* fiber;
	for (int i = 0; i < (int)AllFibers.size(); i++)
	{
		fiber = AllFibers[i];
		ComputeFiberLength(fiber, side);
		CutFiber(fiber, side, amount);
	}
}

//----------------------------------------------------------------------------
// Simply goes trough the fiber and measures how long it is, saving the measurements
// in intermediate points while it measures. Implemented for both sides (0->n and n->0)
void vtkMuscleDecomposer::ComputeFiberLength(Fiber* fiber, int side)
//----------------------------------------------------------------------------
{
	double length = 0;

	if (side == 1)
	{
		PointCloud[fiber->PointIDs[0]]->t = 0;
		for (int i = 1; i < (int)fiber->PointIDs.size(); i++)
		{
			length = vtkMuscleDecomposerUtils::PointPointDistance(PointCloud[fiber->PointIDs[i]]->DCoord, PointCloud[fiber->PointIDs[i - 1]]->DCoord);
			PointCloud[fiber->PointIDs[i]]->t = PointCloud[fiber->PointIDs[i - 1]]->t + length;
		}
	}
	else
	{
		PointCloud[fiber->PointIDs[fiber->PointIDs.size() - 1]]->t = 0;
		for (int i = (int)fiber->PointIDs.size() - 2; i >= 0; i--)
		{
			length = vtkMuscleDecomposerUtils::PointPointDistance(PointCloud[fiber->PointIDs[i]]->DCoord, PointCloud[fiber->PointIDs[i + 1]]->DCoord);
			PointCloud[fiber->PointIDs[i]]->t = PointCloud[fiber->PointIDs[i + 1]]->t + length;
		}
	}
}

//----------------------------------------------------------------------------
// Cuts the fiber by amount percent of its length
void vtkMuscleDecomposer::CutFiber(Fiber* fiber, int side, double amount)
//---------------------------------------------------------------------------
{
	double currentT; // length at the current point
	double targetT; // target length
	double localT; // length inside the last segment (more like parameter t there)

	Point* point1; // segment points
	Point* point2;
	int ID1, ID2; // IDs of segment points

	double vector[3]; // dummy vector
	double k; // specifies how much of amount percent we will use based on a layer deepnes

	// first we cut the segments that are completely further then the target length (both their points are
	// further)
	if (side == 0)
	{
		currentT = PointCloud[fiber->PointIDs[0]]->t;
		k = (fiber->DLayer + 1 / VolumeSubdivision);
		k = k > 1 ? 1 : k;
		targetT = currentT * (1 - amount * k);
		for (int i = 0; i < (int)fiber->PointIDs.size() - 2; i++)
		{
			if (PointCloud[fiber->PointIDs[i + 1]]->t > targetT) // we can cut this segment completely
			{
				UnRegisterOpenEnds(fiber);
				fiber->PointIDs.erase(fiber->PointIDs.begin());
				RegisterOpenEnds(fiber);
			}
		}
		ID1 = fiber->PointIDs[0];
		ID2 = fiber->PointIDs[1];
	}
	else
	{
		currentT = PointCloud[fiber->PointIDs[fiber->PointIDs.size() - 1]]->t;
		targetT = currentT * (1 - amount * fiber->DLayer);
		for (int i = (int)fiber->PointIDs.size() - 1; i > 0; i--)
		{
			//currentT = PointCloud[fiber->PointIDs[i]]->t;
			if (PointCloud[fiber->PointIDs[i - 1]]->t > targetT) // we can cut this segment completely
			{
				UnRegisterOpenEnds(fiber);
				fiber->PointIDs.erase(fiber->PointIDs.end() - 1);
				RegisterOpenEnds(fiber);
			}
		}
		ID1 = fiber->PointIDs[fiber->PointIDs.size() - 1];
		ID2 = fiber->PointIDs[fiber->PointIDs.size() - 2];
	}

	// now, the targetT is in somewhere in the end segment
	point1 = PointCloud[ID1];
	point2 = PointCloud[ID2];
	// we compute how much of the segment we need to cut off
	localT = (targetT - point2->t) / (point1->t - point2->t);
	// and cut it
	vtkMuscleDecomposerUtils::Subtract(point2->DCoord, point1->DCoord, vector);
	vtkMuscleDecomposerUtils::Add(point2->DCoord, vector, vector, localT);
	point1->DCoord[0] = vector[0];
	point1->DCoord[1] = vector[1];
	point1->DCoord[2] = vector[2];
}

//----------------------------------------------------------------------------
// Builds OutputFibers. Physically connects the attachment points.
void vtkMuscleDecomposer::BuildAndOutputConnectedData(std::vector<Fiber*> & Data, char* what, bool connect)
//----------------------------------------------------------------------------
{
	LogFormated("  Building connected %s data...", what);
	int firstID;
	int lastID;

	int outputSize = (int)OutputFibers.size();
	for (int i = 0; i < (int)Data.size(); i++)
	{
		// create output fiber
		OutputFibers.push_back(new Fiber); //delete in ClearMesh()
		FiberCount++;
		// get ID of connections
		firstID = PointCloud[Data[i]->PointIDs[0]]->AttachmentConnectionID;
		lastID = PointCloud[Data[i]->PointIDs[Data[i]->PointIDs.size() - 1]]->AttachmentConnectionID;

		if (firstID > -1 && connect)
			// push first connected tendon point
			OutputFibers[outputSize + i]->PointIDs.push_back(PointCloud[firstID]->Id);

		// then Fiber data
		OutputFibers[outputSize + i]->PointIDs.insert(OutputFibers[outputSize + i]->PointIDs.end(), Data[i]->PointIDs.begin(), Data[i]->PointIDs.end());

		if (lastID > -1 && connect)
			// then push second connected tendon point
			OutputFibers[outputSize + i]->PointIDs.push_back(PointCloud[lastID]->Id);
	}
}

//----------------------------------------------------------------------------
// Outputs a Fiber set
void vtkMuscleDecomposer::OutputFibersData(std::vector<Fiber*> & Data, vtkPoints * PointIDs,
	vtkCellArray * lines, vtkUnsignedCharArray* mtInfo)
	//----------------------------------------------------------------------------
{
	//get the total number of points
	int nTotalPoints = 0;
	for (int i = 0; i < (int)Data.size(); i++) {
		nTotalPoints += (int)Data[i]->PointIDs.size();
	}

	PointIDs->SetNumberOfPoints(nTotalPoints);
	mtInfo->SetNumberOfTuples(nTotalPoints);
	unsigned char* pMti = mtInfo->WritePointer(0, nTotalPoints);

	vtkIdType* lineIDs = NULL; // delete in here
	int nBufSize = 0;

	// go trough all data ready for output

	for (int i = 0, index = 0; i < (int)Data.size(); i++)
	{
		int nPoints = (int)Data[i]->PointIDs.size();
		if (nPoints > nBufSize)
		{
			delete[] lineIDs;
			lineIDs = new vtkIdType[nBufSize = nPoints];
		}

		for (int j = 0; j < nPoints; j++)
		{
			// insert point to output
			const double* coords = PointCloud[Data[i]->PointIDs[j]]->DCoord;
			int mti = PointCloud[Data[i]->PointIDs[j]]->MuscleTendonInfo;

			PointIDs->SetPoint(index, coords);
			if (j == 0 || j == nPoints - 1)
				pMti[index] = 1;	//tendon at the ends OVERIDE
			else if (mti >= 0)
				pMti[index] = mti;
			else
				pMti[index] = 0;	//muscle regularly

			lineIDs[j] = index;
			index++;
		}

		lines->InsertNextCell(nPoints, lineIDs);
	}
	delete[] lineIDs;
	lineIDs = NULL;
}

//// =========================================================================
//// Interpolation methods
//// =========================================================================
//Returns the muscle-tendon info for the point with parameter t from  <0, 1> lying between points identified by ID and nexID identifiers.
int vtkMuscleDecomposer::InterpolateMuscleTendonInfo(int ID, int nextID, double t)
{
	int mtiA = PointCloud[ID]->MuscleTendonInfo;
	int mtiB = PointCloud[nextID]->MuscleTendonInfo;

	if (mtiA < 0 || mtiB < 0)
		return -1;
	else
		return (int)(mtiA*(1 - t) + mtiB*t + 0.5);	//round 0 and 1
}

//----------------------------------------------------------------------------
// Interpolates the surface, reconstructs fibers at surface
void vtkMuscleDecomposer::InterpolateSurface(int subdivision, int interpolationType, std::vector<Fiber*> & Data)
//----------------------------------------------------------------------------
{
	// interpolates the top layer of fibers
	LogFormated("  Interpolating surface points using %s method...", GetSplineType(interpolationType));

	Fiber* crossPoints; // we are interpolating surface across all fibers - these are across all fibers
	std::vector<Fiber*> NewFibers;
	int fiberCount = (int)Data.size();
	int pointCount = (int)Data[0]->PointIDs.size();

	NewFibers.reserve(fiberCount * subdivision);

	SplineInterpolator* interpolator = new SplineInterpolator(interpolationType); // delete in here

	// we go slice by slice
	for (int j = 0; j < (int)Data[0]->PointIDs.size(); j++)
	{
		crossPoints = new Fiber(); // delete in here
		crossPoints->PointIDs.reserve(Data.size());
		// we create a cross slice (from all fibers)
		for (int i = 0; i < fiberCount; i++)
		{
			crossPoints->PointIDs.push_back(Data[i]->PointIDs[j]);
		}
		// interpolate that
		InterpolateFiber(crossPoints, subdivision, interpolator, 0, false);
		// and construct new fibers. Cross slice by cross slice
		for (int i = 0; i < (int)crossPoints->PointIDs.size(); i++)
		{
			if ((int)NewFibers.size() <= i)
			{
				NewFibers.push_back(new Fiber()); // delete later in ClearMesh()
			}
			int ID = crossPoints->PointIDs[i];
			NewFibers[i]->PointIDs.push_back(ID);
		}

		crossPoints->PointIDs.clear();
		delete crossPoints;
		crossPoints = NULL;
	}
	// we need to register these fibers as unconnected
	for (int i = 0; i < (int)NewFibers.size(); i++)
	{
		RegisterOpenEnds(NewFibers[i]);
	}

	// Delete old fibers
	for (int i = 0; i < (int)Data.size(); i++)
	{
		delete Data[i];
	}
	Data.clear();
	// And assign new
	Data.assign(NewFibers.begin(), NewFibers.end());
	delete interpolator;
	interpolator = NULL;
}

//----------------------------------------------------------------------------
// Removes artifacts (short segment). This method is currently UNUSED as we no longer
// use the broken cadaver fibers
void vtkMuscleDecomposer::RemoveArtifacts(double epsilon, std::vector<Fiber*> & data)
//----------------------------------------------------------------------------
{
	Log("  Removing possible data artiffacts...");
	int prevID = 0;
	int currID = 0;
	double segmentLength = 0;
	for (int i = 0; i < (int)Fibers.size(); i++)
	{
		int j = 0;
		while (++j < (int)data[i]->PointIDs.size())
		{
			prevID = data[i]->PointIDs[j - 1];
			currID = data[i]->PointIDs[j];

			segmentLength = vtkMuscleDecomposerUtils::PointPointDistance(this->PointCloud[currID]->DCoord, this->PointCloud[prevID]->DCoord);
			if (segmentLength < epsilon)
				data[i]->PointIDs.erase(data[i]->PointIDs.begin() + j);
		}
	}
}

//----------------------------------------------------------------------------
// Interpolates all the fibers in Data list
void vtkMuscleDecomposer::InterpolateFibers(int subdivision, int interpolationType, std::vector<Fiber*> & Data, bool tendons)
//----------------------------------------------------------------------------
{
	if (tendons)
	{
		LogFormated("  Interpolating tendons using %s method...", GetSplineType(interpolationType));
	}
	else
	{
		LogFormated("\n  Interpolating all individual fibres using %s method...", GetSplineType(interpolationType));
	}
	SplineInterpolator* interpolator = new SplineInterpolator(interpolationType); // delete in here
	for (int i = 0; i < (int)Data.size(); i++)
	{
		// interpolation fiber by fiber
		InterpolateFiber(Data[i], subdivision, interpolator, 0, false);
	}
	delete interpolator;
	interpolator = NULL;
}

//----------------------------------------------------------------------------
// Interpolates (lengthwise) single fiber in Data. Initialised SplineInterpolator
// is needed. Ignores skipOnEnds segments on both sides. If insertSingularity,
// inserts even duplicate points
void vtkMuscleDecomposer::InterpolateFiber(Fiber* fiber, int subdivision, SplineInterpolator *interpolator, int skipOnEnds, bool insertSingularity)
//----------------------------------------------------------------------------
{
	double pCoord[3];
	double* pCoord2; // just a pointer
	Point* newPoint; // delete in clearmesh
	int ID, nextID, newID;
	double distance;
	double step;
	double oldT;
	// we call interpolator to calculate the a,b,c,d coeffs
	interpolator->Interpolate(fiber, PointCloud);

	// reserve memory for faster insertion
	int increase = (int)(fiber->PointIDs.size() * subdivision);
	PointCloud.reserve(PointCloud.size() + increase);
	fiber->PointIDs.reserve(fiber->PointIDs.size() + increase);

	// now we can call point->ComputeInterpolatedCurveCoord(x) for desired x
	for (int j = skipOnEnds; j < (int)fiber->PointIDs.size() - 1 - skipOnEnds; j++)
	{
		ID = fiber->PointIDs[j];
		nextID = fiber->PointIDs[j + 1];
		distance = PointCloud[ID]->DistToNext;
		step = distance / subdivision;
		oldT = PointCloud[ID]->t;

		// if we need to maintain point count, so we insert copies of the start point even when distance == 0
		// (if distance == 0, the solution from interpolator is obviously incorrect)
		if (distance == 0 && insertSingularity)
		{
			for (int k = 0; k < subdivision - 1; k++)
			{
				pCoord2 = PointCloud[ID]->DCoord; // just a pointer
				newPoint = new Point(pCoord2[0], pCoord2[1], pCoord2[2]); //delete in ClearMesh()

				_CHECK_POINT(newPoint->DCoord);

				newPoint->LeftAnchorID = ID;
				newPoint->RightAnchorID = nextID;
				// insert it to cloud
				PointCloud.push_back(newPoint);
				newID = (int)PointCloud.size() - 1;
				newPoint->Id = newID;
				newPoint->MuscleTendonInfo = PointCloud[ID]->MuscleTendonInfo;
				// and to the list of IDs for this fibre
				j++;
				fiber->PointIDs.insert(fiber->PointIDs.begin() + j, newID);
			}
		}
		else // interpolate normally
		{
			double newT = oldT;
			for (int l = 0; l < subdivision - 1; l++) // more stable
			//for (double newT = oldT + step; abs(newT - oldT) < abs(distance); newT += step)
			{
				newT += step;
				// compute point position

				PointCloud[ID]->ComputeInterpolatedCurveCoord(newT, pCoord, interpolator->GetInterpolationType() != SPLINE_CONSTRAINED && interpolator->GetInterpolationType() != SPLINE_CATMULL_ROM); //delete here
				newPoint = new Point(pCoord[0], pCoord[1], pCoord[2]); //delete in ClearMesh()
				_CHECK_POINT(newPoint->DCoord);

				newPoint->LeftAnchorID = ID;
				newPoint->RightAnchorID = nextID;
				newPoint->t = newT;

				// insert it to cloud
				PointCloud.push_back(newPoint);
				newID = (int)PointCloud.size() - 1;
				newPoint->Id = newID;

				//interpolate linearly MuscleTendonInfo
				newPoint->MuscleTendonInfo = InterpolateMuscleTendonInfo(ID, nextID, newT / distance);

				// and to the list of IDs for this fibre
				j++;
				fiber->PointIDs.insert(fiber->PointIDs.begin() + j, newID);
			}
		}
	}
}

//----------------------------------------------------------------------------
// Computes distance between two PointIDs
double vtkMuscleDecomposer::ComputePointDistance(int ID1, int ID2)
//----------------------------------------------------------------------------
{
	Point* First = PointCloud[ID1];
	Point* Second = PointCloud[ID2];

	return vtkMath::Distance2BetweenPoints(First->DCoord, Second->DCoord);
}