/*=========================================================================
Program:   Issue 120 Test Application
Module:    $RCSfile: Issue_120.cxx,v $
Language:  C++
Date:      $Date: 2011-08-07 20:51 $
Version:   $Revision: 0.1.0.0 $
Author:    David Cholt
Notes:
=========================================================================*/

//----------------------------------------------------------------------------
// Include:
//----------------------------------------------------------------------------
#ifndef __unix__
#include <Windows.h>
#define SLASH "\\"
#else
#include <sys/types.h>
#include <dirent.h>
#define SLASH "/"
#define sprintf_s snprintf
#define fopen_s(file,name,mode) ((*(file))=fopen((name),(mode)))==NULL
typedef unsigned long DWORD;
#endif

#include "vtkConfigure.h"

#if VTK_MAJOR_VERSION >= 6
//this is required to initialize rendering
#include <vtkAutoInit.h>
#if VTK_MAJOR_VERSION <= 6
VTK_MODULE_INIT(vtkRenderingOpenGL);
#endif
#endif

#include "vtkPolyDataReader.h"
#include "vtkPolyDataWriter.h"

#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkCamera.h"

#include "vtkRenderer.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkOBJReader.h"

#include <string>
#include "vtkCommand.h"
#include "vtkSmartPointer.h"
#include "vtkCallbackCommand.h"
#include "vtkWindowToImageFilter.h"
#include "vtkPNGWriter.h"
#include "vtkMuscleDecomposer.h"
#include "vtkFileOutputWindow.h"
#include "vtkOutputWindow.h"
#include <vector>

#include "vtkMapperActor.h"
#include "vtkReverseSense.h"
#include "vtkSmartPointer.h"
#include "vtkTimerLog.h"
#include "vtkPointData.h"
#include "vtkCellData.h"

#include "vtkLHPMuscleFibresAnalysisFilter.h"
#include "vtkLHPMuscleFibresResample.h"
#include "vtkLHPMuscleFibresMath.h"
#include "vtkDoubleArray.h"
#include "vtkMassProperties.h"
#include "vtkSphereSource.h"
#include "vtkGlyph3D.h"
#include "vtkLookupTable.h"
#include "vtkLHPPolyLinesBuilder.h"
#include "vtkAppendPolyData.h"
#include "vtkCleanPolyData.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkTransform.h"

#include <sstream>

const char* TENDON_PART_STRING = "tendon_part";

vtkRenderWindow *renWin;
vtkRenderer* ren1;
bool showBones = true;
bool showNormals = false;
bool showFibres = true;
bool showMuscle = true;
bool showInputFibres = true;
bool showInputAttachments = true;
bool showFibreSize = true;
#ifdef _DEBUG
bool showLinesOfAction = true;
#endif

#define VISUALGROUP_BONE					0
#define VISUALGROUP_MUSCLE					1
#define VISUALGROUP_INFIBRES				2
#define VISUALGROUP_INATTACH				3
#define VISUALGROUP_NORMALS					4
#define VISUALGROUP_RESULT					5
#define VISUALGROUP_LINEOFACTION			6

//#define CMPB_JOURNAL	//measurements for CMPB journal (2017)

#ifdef _DEBUG
vtkPolyData* GetLineOfActionObj(const double* dir, const double* point, double size, bool onesideonly = false)
{
	vtkPoints* points = vtkPoints::New();
	if (onesideonly)
		points->InsertNextPoint(point[0], point[1], point[2]);
	else
		points->InsertNextPoint(point[0] - size*dir[0], point[1] - size*dir[1], point[2] - size*dir[2]);
	points->InsertNextPoint(point[0] + size*dir[0], point[1] + size*dir[1], point[2] + size*dir[2]);


	vtkIdType edge[2] = { 0, 1 };

	vtkCellArray* cells = vtkCellArray::New();
	cells->InsertNextCell(2, edge);

	vtkPolyData* poly = vtkPolyData::New();
	poly->SetPoints(points);
	poly->SetLines(cells);
	points->Delete();
	cells->Delete();

	vtkTubeFilter* tube = vtkTubeFilter::New();
#if VTK_MAJOR_VERSION < 5
	tube->SetInput(poly);
#else
	tube->SetInputData(poly);
#endif
	tube->SetNumberOfSides(6);
	tube->SetRadius(2*size / 100);
	poly->Delete();

	vtkSmartPointer< vtkAppendPolyData > append = vtkSmartPointer< vtkAppendPolyData >::New();
	append->AddInputConnection(tube->GetOutputPort());

	if (onesideonly)
	{
		vtkSmartPointer< vtkSphereSource > sphere = vtkSmartPointer< vtkSphereSource >::New();
		sphere->SetRadius(size / 25);

		vtkSmartPointer< vtkTransform > tr = vtkSmartPointer< vtkTransform >::New();
		tr->Translate(point);

		vtkSmartPointer< vtkTransformPolyDataFilter > trf = vtkSmartPointer< vtkTransformPolyDataFilter >::New();
		trf->SetInputConnection(sphere->GetOutputPort());
		trf->SetTransform(tr);

		append->AddInputConnection(trf->GetOutputPort());
	}

	poly = vtkPolyData::New();
	append->SetOutput(poly);
	append->Update();
	append->SetOutput(NULL);

	tube->Delete();
	return poly;
}
#endif

//Detect all the tendon points in retData and construct TENDON_PART_STRING unsigned char scalar 
//field with -1 = UNKNOWN, 0 = MUSCLE and 1 = TENDON point
void MarkTendonPoints(const std::vector<vtkPolyData*>& tendons, vtkPolyData* retData)
{
	vtkSmartPointer<vtkUnsignedCharArray> mtInfo = vtkSmartPointer<vtkUnsignedCharArray>::New();
	vtkDataArray* tpd = retData->GetPointData()->GetArray(TENDON_PART_STRING);
	if (tpd != NULL)
		mtInfo = vtkUnsignedCharArray::SafeDownCast(tpd);
	else
	{
		mtInfo->SetName(TENDON_PART_STRING);
		retData->GetPointData()->AddArray(mtInfo);
	}

	int nTotalPoints = retData->GetNumberOfPoints();
	unsigned char* pMti = mtInfo->WritePointer(0, nTotalPoints);
	memset(pMti, 0, sizeof(unsigned char)*nTotalPoints);	//0 = MUSCLE

	//check every point if it comes from tendon part or not
	size_t count = tendons.size();
	for (int i = 0; i < nTotalPoints; i++)
	{
		double point[3], point2[3];
		retData->GetPoint(i, point);

		bool isTendon = false;
		for (size_t j = 0; j < count && !isTendon; j++)
		{
			int nPoints = tendons[j]->GetNumberOfPoints();
			for (int k = 0; k < nPoints; k++)
			{
				tendons[j]->GetPoint(k, point2);
				if (point[0] == point2[0] && point[1] == point2[1] && point[2] == point2[2]) {
					isTendon = true; break;
				}
			}
		}

		if (isTendon) {
			pMti[i] = 1;
		}
	}
}

void ConnectTendonParts(vtkPolyData* retData)
{
	vtkSmartPointer<vtkUnsignedCharArray> mtInfo =
		vtkUnsignedCharArray::SafeDownCast(retData->GetPointData()->GetArray(TENDON_PART_STRING));

	int nTotalPoints = retData->GetNumberOfPoints();
	unsigned char* pMti = mtInfo->GetPointer(0);

	//merge tendon cells with muscle cells
	//TODO: 

	//connect isolated tendon points
	retData->SetVerts(NULL);	//remove unwanted cells
	retData->BuildCells();
	retData->BuildLinks();

	//for each cell we have state: 
	//0 = NOTHING TO CHANGE
	//1 = PREPEND POINT (at 0)
	//2 = APPEND POINT (at 1)
	//3 = PREPEND POINT (at 0) and APPEND POINT (at 1)
	int nCells = retData->GetNumberOfCells();
	std::pair<int, std::pair<int, int>>* cellState = new std::pair<int, std::pair<int, int>>[nCells];
	memset(cellState, 0, sizeof(cellState[0])*nCells);

	int nPointsAdded = 0;
	for (int i = 0; i < nTotalPoints; i++)
	{
		if (pMti[i] != 1)	//1 = tendon
			continue;

		vtkIdType* pCls;
		unsigned short nCls;

		retData->GetPointCells(i, nCls, pCls);
		if (nCls == 0)
		{
			//the point is in no cell, so it is isolated, thus find its cell
			double x[3], mindist2 = DBL_MAX;
			int closestCell, appendPoint;	//0 = prepend, 1 = append

			retData->GetPoint(i, x);
			for (int j = 0; j < nCells; j++)
			{
				vtkIdType nPts, *pPts;
				retData->GetCellPoints(j, nPts, pPts);
				if ((cellState[j].first & 1) == 0)
				{
					//test if we could prepend the point to the cell
					double dist2 = vtkMath::Distance2BetweenPoints(x, retData->GetPoint(pPts[0]));
					if (dist2 < mindist2)
					{
						mindist2 = dist2;
						closestCell = j;
						appendPoint = 0;
					}
				}

				if ((cellState[j].first & 2) == 0)
				{
					//test if we could append the point to the cell
					double dist2 = vtkMath::Distance2BetweenPoints(x, retData->GetPoint(pPts[nPts - 1]));
					if (dist2 < mindist2)
					{
						mindist2 = dist2;
						closestCell = j;
						appendPoint = 1;
					}
				}
			} //end for cells

			if (appendPoint == 0) {
				cellState[closestCell].first |= 1; cellState[closestCell].second.first = i;
			}
			else {
				cellState[closestCell].first |= 2; cellState[closestCell].second.second = i;
			}

			nPointsAdded++;
		}
	}

	//now create new cells
	if (nPointsAdded != 0)
	{
		const vtkIdType* pInPtr = retData->GetLines()->GetPointer();

		vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
		int sz = retData->GetLines()->GetSize() + nPointsAdded;

		vtkIdType* pOutPtr = cells->WritePointer(nCells, sz);
		for (int i = 0; i < nCells; i++)
		{
			vtkIdType* pNumPtr = pOutPtr;
			*pOutPtr = *pInPtr;	//copy number of points in the cell
			pOutPtr++;

			if ((cellState[i].first & 1) != 0)
			{
				(*pNumPtr)++;	//increase the number of points
				*pOutPtr = cellState[i].second.first;
				pOutPtr++;
			}

			memcpy(pOutPtr, pInPtr + 1, sizeof(*pInPtr)*(*pInPtr));
			pOutPtr += *pInPtr;
			pInPtr += (*pInPtr) + 1;

			if ((cellState[i].first & 2) != 0)
			{
				(*pNumPtr)++;	//increase the number of points
				*pOutPtr = cellState[i].second.second;
				pOutPtr++;
			}
		}

		retData->SetLines(cells);

		//internal Cells and Links structures are no longer valid
		retData->DeleteCells();
		retData->DeleteLinks();
	}

	delete[] cellState;
}

//----------------------------------------------------------------------------
// Gets the file extension from filename
std::string GetFileExtension(const std::string& FileName)
//----------------------------------------------------------------------------
{
	std::string tmp;
	unsigned int i;
	if (FileName.find_last_of(".") != std::string::npos) {
		tmp = FileName.substr(FileName.find_last_of(".") + 1);
		for (i = 0; i < strlen(tmp.c_str()); i++)
			if (tmp[i] >= 0x41 && tmp[i] <= 0x5A)
				tmp[i] = tmp[i] + 0x20;
		return tmp;
	}
	return "";
}

//----------------------------------------------------------------------------
// Checks if the file exists
bool FileExists(const char* FileName)
//----------------------------------------------------------------------------
{
	FILE* fp = NULL;
	if (fopen_s(&fp, FileName, "rb") == 0)
	{
		fclose(fp);
		return true;
	}

	return false;
}

//----------------------------------------------------------------------------
// Checks, if specified directory exists
bool DirectoryExists(const char* path)
//----------------------------------------------------------------------------
{
#ifndef __unix__
	DWORD ftyp = GetFileAttributesA(path);
	if (ftyp == INVALID_FILE_ATTRIBUTES)
		return false;  //something is wrong with your path!

	if (ftyp & FILE_ATTRIBUTE_DIRECTORY)
		return true;   // this is a directory!

	return false;    // this is not a directory!
#else
	DIR* dir = opendir(path);
	if (dir) {
		closedir(dir);
		return true;
	}
	else {
		return false;
	}
#endif
}

std::string CombinePathName(const char* path, const char* identifier, int index, const char* ext = ".vtk")
{
	std::string ret = path;
	ret.append(SLASH);
	ret.append(identifier);

	char filename[10];
	std::snprintf(filename, sizeof(filename) / sizeof(filename[0]), "%d", index);
	ret.append(filename);
	ret.append(ext);
	return ret;
}

//----------------------------------------------------------------------------
// Loads requested data, identified by the identifier and index
vtkPolyData* LoadData(const char* path, const char* identifier, int index)
//----------------------------------------------------------------------------
{
	std::string filename = CombinePathName(path, identifier, index);
	if (!FileExists(filename.c_str())) {	//invalid file
		cout << "ERROR - file " << filename << " does not exists" << endl;
		return NULL;
	}

	vtkSmartPointer< vtkPolyDataReader > reader = vtkPolyDataReader::New();
	reader->Delete();	//decrease the number of references to 1

	reader->SetFileName(filename.c_str());
	if (reader->OpenVTKFile() != 0)
	{
		reader->CloseVTKFile();
		reader->Update();	//this calls OpenVTKFile and CloseVTKFile => must call CloseVTKFile first

		//detach the data from the reader
		vtkPolyData* data = reader->GetOutput();
		data->Register(NULL);
		reader->SetOutput(NULL);

		return data;
	}
	else
	{
		cout << "ERROR - file " << filename << " does not exists" << endl;
		return NULL;
	}
}

//----------------------------------------------------------------------------
// Loads and fix the requested fibres data, identified by the path and index
vtkPolyData* LoadFibresData(const char* path, int index)
{
	vtkPolyData* data = LoadData(path, "f", index);
	if (data == NULL)
		return NULL;

	//we have a valid data, check it	
	vtkSmartPointer< vtkLHPPolyLinesBuilder > mt = vtkLHPPolyLinesBuilder::New();
	mt->Delete();	//decrease the reference count to 1

	mt->SetInputData(data);
	mt->Update();

	//detach output	
	vtkPolyData* newdata = mt->GetOutput();
	newdata->Register(NULL);
	mt->SetOutput(NULL);

	if (data->GetLines() != newdata->GetLines())
	{
		//SAVE the new DATA
		vtkSmartPointer< vtkPolyDataWriter > writer = vtkPolyDataWriter::New();
		writer->Delete();	//decrease the reference count to 1

		writer->SetFileTypeToASCII();
		writer->SetFileName(CombinePathName(path, "f", index, ".new.vtk").c_str());
		writer->SetInputData(newdata);
		writer->Update();
	}

	data->Delete();	//no longer needed
	return newdata;
}

vtkPolyData* LoadFibresData(const char* path)
{
	vtkSmartPointer< vtkAppendPolyData > append = vtkSmartPointer< vtkAppendPolyData >::New();
	for (int i = 0; true; i++)
	{
		vtkPolyData* data = LoadFibresData(path, i);
		if (data == NULL) break;

		append->AddInputData(data);
	}

	//try to load tendons	
	std::vector< vtkPolyData* > tendons;
	for (int i = 0; true; i++)
	{
		vtkPolyData* t0 = LoadData(path, "at", i);
		if (t0 == NULL) break;

		append->AddInputData(t0);
		tendons.push_back(t0);
	}

	for (int i = 0; true; i++)
	{
		vtkPolyData* t0 = LoadData(path, "bt", i);
		if (t0 == NULL) break;

		append->AddInputData(t0);
		tendons.push_back(t0);
	}

	vtkSmartPointer< vtkCleanPolyData > clean = vtkSmartPointer< vtkCleanPolyData >::New();
	clean->SetInputConnection(append->GetOutputPort());
	clean->PointMergingOn();
	clean->Update();

	vtkPolyData* retData = clean->GetOutput();
	retData->Register(NULL);
	clean->SetOutput(NULL);

	//now connect the tendons, if available	
	size_t count = tendons.size();
	if (count != 0)
	{
		MarkTendonPoints(tendons, retData);

		ConnectTendonParts(retData);

		//clean up - not needed because LoadData place everything to deletables
		//for (size_t i = 0; i < count; i++) {
		//	tendons[i]->Delete();
		//}

		//SAVE the new DATA
		vtkSmartPointer< vtkPolyDataWriter > writer = vtkPolyDataWriter::New();
		writer->Delete();	//decrease the reference count to 1

		writer->SetFileTypeToASCII();
		writer->SetFileName(CombinePathName(path, "ft", 0, ".new.vtk").c_str());
		writer->SetInputData(retData);
		writer->Update();
	}

	return retData;
}

//----------------------------------------------------------------------------
// Saves single data object to VTK file
void SaveData(vtkPolyData* data, const char* path, const char* name, int index)
//----------------------------------------------------------------------------
{
	vtkSmartPointer< vtkPolyDataWriter > writer = vtkSmartPointer< vtkPolyDataWriter >::New();

	int len = (int)(10 + strlen(path) + strlen(name));
	char* filename = new char[len];

	sprintf_s(filename, len, "%s%s%s%d%s", path, SLASH, name, index, ".vtk");
	writer->SetFileName(filename);
#if VTK_MAJOR_VERSION < 5
	writer->SetInput(data);
#else
	writer->SetInputData(data);
#endif
	writer->Write();

	delete[] filename;
}

#include "vtkParametricCatmullRomSurface.h"
#include "vtkParametricFunctionSource.h"
#include "vtkArrowSource.h"
#include "vtkPolyDataNormals.h"
#include "vtkOBBTree.h"
#include "vtkCellLocator.h"
#include "demoApp.h"
#include "vtkGenericCell.h"

template < typename T >
void GetAvgNormalT(T* normals, int count, double avgNormal[3])
{
	avgNormal[0] = avgNormal[1] = avgNormal[2] = 0.0;
	if (count == 0)
		return;	//array is empty

	const T* pPtr = normals;
	for (int i = 0; i < count; i++)
	{
		for (int k = 0; k < 3; k++)
		{
			avgNormal[k] += (double)(*pPtr);
			pPtr++;	//move to next
		}
	}

	for (int k = 0; k < 3; k++) {
		avgNormal[k] /= count;
	}
}

//Computes the average vector (normal) of the given vector (normal) field.
void GetAvgNormal(vtkDataArray* normals, double avgNormal[3])
{
	switch (normals->GetDataType())
	{
		vtkTemplateMacro(
			GetAvgNormalT(static_cast<VTK_TT*>(normals->GetVoidPointer(0)),
				normals->GetNumberOfTuples(), avgNormal)
		);
	}
}

//Constructs the Catmull-Rom surface patch interpolating the input fibres.
//ensuring that the orientation of its cells and its normals is consistent
//with the orientation of the cells(and normals) of the muscle surface.
void ConstructFibresPatch(vtkPolyData* fibresData, vtkPolyData* muscleData, vtkPolyData* out)
{
	vtkSmartPointer< vtkParametricCatmullRomSurface > cr = vtkSmartPointer< vtkParametricCatmullRomSurface >::New();
	cr->SetPointsGrid(fibresData);

	vtkSmartPointer<vtkParametricFunctionSource> parametricFunctionSource =
		vtkSmartPointer<vtkParametricFunctionSource>::New();
	parametricFunctionSource->SetUResolution(16);
	parametricFunctionSource->SetVResolution(16);
	parametricFunctionSource->SetParametricFunction(cr);

	vtkSmartPointer<vtkReverseSense> nrminv = vtkSmartPointer<vtkReverseSense>::New();
	nrminv->SetInputConnection(parametricFunctionSource->GetOutputPort());
	nrminv->SetReverseNormals(0);	//just pass the data to output
	nrminv->SetReverseCells(0);

	nrminv->SetOutput(out);
	nrminv->Update();

	//check the orientation of the normals
	//as the input fibres are supposed to be superficial, they should lie close to the 
	//muscle surface, and, furthermore, only to one side of the muscle surface
	//the normals of the fibres patch must go outside of the muscle surface
	//normals of the muscle surface are supposed to be correct
	if (!FibresPatchNormalsOK(out, muscleData))
	{
		//invert normals
		cout << "Wrong orientation of fibres surface normals detected. Normals will be flipped.\n";

		nrminv->SetReverseNormals(1);
		nrminv->SetReverseCells(1);
		nrminv->Update();
	}

	//detach the output polydata from the pipeline
	nrminv->SetOutput(NULL);
}

//Computes the normals at points of the given polydata.
//N.B. the caller is responsible for deleting the returned object.
vtkDataArray* ComputePointNormals(vtkPolyData* polyData)
{
	vtkSmartPointer<vtkPolyDataNormals> nf = vtkSmartPointer<vtkPolyDataNormals>::New();
	nf->SetInputData(polyData);
	nf->SetComputePointNormals(1);
	nf->Update();

	vtkDataArray* ret = nf->GetOutput()->GetPointData()->GetNormals();
	ret->Register(NULL);

	return ret;
}

//Returns true if the orientation of the normals of the fibres patch surface
//is consistent with the orientation of the normals of the muscle surface.
bool FibresPatchNormalsOK(vtkPolyData* fibresPoly, vtkPolyData* musclePoly)
{
	//muscle surface normals are oriented to outwards the muscle
	//as the fibres patch interpolates superficial fibres, it should lie 
	//close to the muscle surface and its normals must be consistently
	//oriented with the normals of the muscle surface in its proximity
	vtkSmartPointer<vtkDataArray> fbNormals = fibresPoly->GetPointData()->GetNormals();
	vtkSmartPointer<vtkDataArray> mscNormals = musclePoly->GetPointData()->GetNormals();
	if (mscNormals == NULL)
	{
		mscNormals = ComputePointNormals(musclePoly);
		mscNormals->Delete();	//decrease reference count to 1
	}

	//as some points of the fibres patch may lie completely "outside" the muscle
	//we need to get some distance tolerance to remove such points
	double obbsz[3];
	vtkOBBNode nodeMuscle;
	memset(&nodeMuscle, sizeof(nodeMuscle), 0);

	vtkOBBTree::ComputeOBB(musclePoly->GetPoints(), nodeMuscle.Corner,
		nodeMuscle.Axes[0], nodeMuscle.Axes[1], nodeMuscle.Axes[2], obbsz);

	double tol2 = obbsz[2] * 0.125;
	tol2 *= tol2;	//square

	//compute an OBB tree for the muscle surface to speed up the detection
	vtkSmartPointer< vtkCellLocator > locator = vtkSmartPointer< vtkCellLocator >::New();
	locator->SetDataSet(musclePoly);
	locator->Update();

	int nPosOK = 0;

	vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();
	double* weights = NULL;
	int nw = 0;

	int nPoints = fibresPoly->GetNumberOfPoints();
	for (int i = 0; i < nPoints; i++)
	{
		int subId;
		vtkIdType cellId;
		double xF[3], xM[3], dist2;
		fibresPoly->GetPoint(i, xF);

		//for the given fibres patch point find the closest point on the muscle surface
		locator->FindClosestPoint(xF, xM, cell, cellId, subId, dist2);
		if (dist2 > tol2)
			continue;	//the point is too far

		//find the relative location of the closest point within the cell
		int nCellPts = cell->GetNumberOfPoints();
		if (nCellPts > nw)
		{
			delete[] weights;
			weights = new double[nw = nCellPts];
		}

		double xMC[3], pcoords[3];
		cell->EvaluatePosition(xM, xMC, subId, pcoords, dist2, weights);

		//interpolate the normals at cell points to get the normal at
		//the point xMC, which is closest to the fibres patch
		double msn[3]{ 0, 0, 0 };
		for (int j = 0; j < nCellPts; j++)
		{
			double ptn[3];
			mscNormals->GetTuple(cell->GetPointId(j), ptn);

			for (int k = 0; k < 3; k++) {
				msn[k] += ptn[k] * weights[j];
			}
		}

		//get the normal of the fibre point
		double fbn[3];
		fbNormals->GetTuple(i, fbn);
		double dot = vtkMath::Dot(fbn, msn);
		if (dot > 0) {
			nPosOK++;
		}
		else {
			nPosOK--;
		}
	} //end for i

	delete[] weights;

	return nPosOK > 0;
}

//Moves the fibres patch in the given direction to make sure that
//the transformed patch does not intersect the muscle.
void TransformFibresPatch(vtkPolyData* fibresPoly, vtkPolyData* muscleData,
	const double dir[3], vtkPolyData* out)
{
	//get OBB of the fibres patch and of the muscle
	double obbsz[3];
	vtkOBBNode nodeFib, nodeMuscle;
	memset(&nodeFib, sizeof(nodeFib), 0);
	memset(&nodeMuscle, sizeof(nodeMuscle), 0);

	vtkOBBTree::ComputeOBB(fibresPoly->GetPoints(), nodeFib.Corner,
		nodeFib.Axes[0], nodeFib.Axes[1], nodeFib.Axes[2], obbsz);

	//obbsz contains the sizes of the OBB bounding box ordered from max to min
	double shiftDistDelta = obbsz[2] / 64;	//the smallest distance we will shift the fibres patch

	vtkOBBTree::ComputeOBB(muscleData->GetPoints(), nodeMuscle.Corner,
		nodeMuscle.Axes[0], nodeMuscle.Axes[1], nodeMuscle.Axes[2], obbsz);

	if (shiftDistDelta == 0.0)
		shiftDistDelta = obbsz[2] / 128;		//this happens if the fibres lie on a plane

//we now need to make sure that the OBBs nodes do not intersect
//and then transform rigidly the current surface patch to be located fully above the muscle
	vtkSmartPointer< vtkTransform > tr = vtkSmartPointer< vtkTransform >::New();
	vtkSmartPointer< vtkTransformPolyDataFilter > trf = vtkSmartPointer< vtkTransformPolyDataFilter >::New();
	trf->SetInputData(fibresPoly);
	trf->SetTransform(tr);
	trf->SetOutput(out);

	double uTr[3];
	for (int k = 0; k < 3; k++) {
		uTr[k] = dir[k] * shiftDistDelta;
	}

	vtkSmartPointer<vtkOBBTree> obbMuscle = vtkSmartPointer<vtkOBBTree>::New();
	while (obbMuscle->DisjointOBBNodes(&nodeFib, &nodeMuscle, NULL) == 0)
	{
		for (int k = 0; k < 3; k++) {
			nodeFib.Corner[k] += uTr[k];
		}

		tr->Translate(uTr);
	}

	trf->Update();
	trf->SetOutput(NULL);		//detach output from the pipeline
}

void RegisterFibresAndMesh(vtkPolyData* fibresData, vtkPolyData* muscleData
#ifdef _DEBUG
	, const char* path
#endif
)
{
	//Algorithm:
	//1) Construct fibres surface patch
	//2) Move the patch along its principal (average) normal so that the patch
	//	 no longer intersects the muscle
	//3) Use position based dynamics to move the vertices of the patch 
	//	 along the oposite direction of the principal normal toward the muscle
	//	 while handing collisions and self-collisions and minimizing
	//	 internal energy (to avoid strange deformations of the patch)
	//4) When the method converges, the patch should lie on the muscle

	//construct Catmull-Rom surface patch interpolating the input superficial fibres
	vtkSmartPointer<vtkPolyData> fibresPoly = vtkSmartPointer<vtkPolyData>::New();
	ConstructFibresPatch(fibresData, muscleData, fibresPoly);

	//determine the principal (average) normal of the constructed patch
	//this normal defines the opposite direction to the direction, in 
	//which the fibres surface should be aligned with the muscle surface.	
	double avgNormal[3];
	GetAvgNormal(fibresPoly->GetPointData()->GetNormals(), avgNormal);
	vtkMath::Normalize(avgNormal);

	vtkSmartPointer<vtkPolyData> fibresPolyTr = vtkSmartPointer<vtkPolyData>::New();
	TransformFibresPatch(fibresPoly, muscleData, avgNormal, fibresPolyTr);

#ifdef _DEBUG	
	//for the debugging purpose, save the constructed fibres patch
	SaveData(fibresPoly, path, "ffPatch", 0);

	//and save its normals so they could be displayed
	vtkSmartPointer< vtkArrowSource > glyphArrow = vtkSmartPointer< vtkArrowSource >::New();
	vtkSmartPointer< vtkGlyph3D > glyphs = vtkSmartPointer< vtkGlyph3D >::New();
	glyphs->SetInputData(fibresPoly);

	glyphs->SetScaleModeToScaleByVector();
	glyphs->SetScaleFactor(3.0);
	glyphs->SetVectorModeToUseNormal();	//must be placed after SetScaleModeToScaleByVector
	glyphs->SetSourceConnection(glyphArrow->GetOutputPort());
	glyphs->Update();

	SaveData(glyphs->GetOutput(), path, "ffPatchNormals", 0);

	//and save OBB of the patch and the muscle
	vtkSmartPointer< vtkOBBTree > obbFibres = vtkSmartPointer< vtkOBBTree >::New();
	obbFibres->SetDataSet(fibresPoly);
	obbFibres->Update();

	vtkSmartPointer< vtkPolyData > obbRepr = vtkSmartPointer< vtkPolyData >::New();
	obbFibres->GenerateRepresentation(0, obbRepr);
	SaveData(obbRepr, path, "ffPatchObb", 0);

	vtkSmartPointer< vtkOBBTree > obbMuscle = vtkSmartPointer< vtkOBBTree >::New();
	obbMuscle->SetDataSet(muscleData);
	obbMuscle->Update();

	obbMuscle->GenerateRepresentation(0, obbRepr);
	SaveData(obbRepr, path, "mMscObb", 0);

	//save the principal direction in which we move the fibres patch
	vtkPolyData* dirPrinc = GetLineOfActionObj(avgNormal, fibresPoly->GetPoint(0), 50, true);
	SaveData(dirPrinc, path, "ffPatchDir", 0);
	dirPrinc->Delete();

	//for the debugging purpose, save the transformed fibres patch
	SaveData(fibresPolyTr, path, "ffPatchTR", 0);
#endif
}

static vtkSphereSource* sphere = NULL;
static vtkPolyData* g_MuscleMesh = NULL;

//----------------------------------------------------------------------------
// Loads and adds all necessary data to the decomposer
void LoadAllData(vtkMuscleDecomposer* musdec, const char* path)
//----------------------------------------------------------------------------
{
	vtkMapperActor* boneMA;

	vtkPolyData* data;
	vtkPolyData* bone;

	char identifier[3];
	vtkPolyData* fibresData = LoadFibresData(path);
	musdec->AddFiberData(fibresData);

	vtkTubeFilter* tuber = vtkTubeFilter::New();
#if VTK_MAJOR_VERSION < 5
	tuber->SetInput(fibresData);
#else
	tuber->SetInputData(fibresData);
#endif

	tuber->SetRadius(0.8);
	tuber->SetNumberOfSides(6);

	//		vtkSphereSource* sphere = vtkSphereSource::New();
	if (sphere == NULL)
		sphere = vtkSphereSource::New();

	sphere->SetRadius(1.2);
	sphere->SetPhiResolution(6);
	sphere->SetThetaResolution(6);

	vtkGlyph3D* glyphs = vtkGlyph3D::New();
#if VTK_MAJOR_VERSION < 5
	glyphs->SetInput(fibresData);
	glyphs->SetSource(sphere->GetOutput());
#else
	glyphs->SetInputData(fibresData);
	glyphs->SetSourceConnection(sphere->GetOutputPort());
#endif

	tuber->SetVaryRadius(0);
	glyphs->SetScaleModeToDataScalingOff();

	//sphere->Delete();

	tuber->Update();	//make Output
	new vtkMapperActor(tuber->GetOutput(), 255, 192, 0, 255, false, ren1,
		VISUALGROUP_INFIBRES, showInputFibres);

	glyphs->Update();	//make Output
	new vtkMapperActor(glyphs->GetOutput(), 0, 0, 255, 255, false, ren1,
		VISUALGROUP_INFIBRES, showInputFibres);

	glyphs->Delete();
	tuber->Delete();
	printf("  Fibers loaded\n");
	cout << "n = " << fibresData->GetNumberOfCells() << ", m = " <<
		fibresData->GetNumberOfPoints() / fibresData->GetNumberOfCells() - 1 << "\n";

	identifier[2] = 0;
	for (int side = 0; side < 2; side++)
	{
		for (int i = 0; true; i++)
		{
			identifier[0] = 'a' + side;
			identifier[1] = 'a';
			data = LoadData(path, identifier, i);
			if (data != NULL)
			{
				musdec->AddAttachmentData(data, NULL, side);

				new vtkMapperActor(data, 72, 0, 255, 255, false, ren1,
					VISUALGROUP_INATTACH, showInputAttachments);
			}
			else
			{
				identifier[0] = 'a' + side;
				identifier[1] = 'o';
				data = LoadData(path, identifier, i);
				if (data == NULL) break;
				identifier[1] = 'b';
				bone = LoadData(path, identifier, i);
				if (bone != NULL)
				{
					printf("    Some pre-generated attachment areas not found, regenerating...\n");
					musdec->AddAttachmentData(data, bone, side);
				}
			}

		}

		for (int i = 0; true; i++)
		{
			identifier[0] = 'a' + side;
			identifier[1] = 'b';
			bone = LoadData(path, identifier, i);
			if (bone == NULL) break;

			boneMA = new vtkMapperActor(bone, 242, 239, 206, 255, false, ren1,
				VISUALGROUP_BONE, showBones);
		}
	}
	printf("  Attachments loaded\n");
	data = LoadData(path, "m", 0);
	if (data != NULL)
	{
		musdec->AddMeshData(data);
		printf("  Mesh loaded\n");

		g_MuscleMesh = data;
	}

	/*	RegisterFibresAndMesh(fibresData, data
	#ifdef _DEBUG
			, path
	#endif
		);
		printf("  Fibres registered to the mesh\n");
		*/
}

//----------------------------------------------------------------------------
// Saves all attachment data and output polylines to VTK files
void SaveAllData(vtkMuscleDecomposer* musdec, const char* path)
//----------------------------------------------------------------------------
{
	char identifier[3];
	std::vector<vtkPolyData*>* attachments;

	SaveData(musdec->GetOutput(), path, "output", 0);

	identifier[2] = 0;
	attachments = musdec->GetAttachmentsSurfaces();
	for (int side = 0; side < 2; side++)
	{
		for (int i = 0; i < (int)attachments[side].size(); i++)
		{
			identifier[0] = 'a' + side;
			identifier[1] = 'a';
			SaveData(attachments[side][i], path, identifier, i);
		}
	}
}

//----------------------------------------------------------------------------
// Saves the current camera view
void SaveCamera(vtkCamera* cam, const char* path)
{
	std::string fname = CombinePathName(path, "camera", 0, ".bin");

	FILE* f = NULL;
	if (0 != fopen_s(&f, fname.c_str(), "wb"))
		return;

	const double* pos = cam->GetPosition();
	fwrite(pos, sizeof(double), 3, f);

	pos = cam->GetFocalPoint();
	fwrite(pos, sizeof(double), 3, f);

	double val = cam->GetViewAngle();
	fwrite(&val, sizeof(double), 1, f);

	pos = cam->GetViewUp();
	fwrite(pos, sizeof(double), 3, f);

	val = cam->GetDistance();
	fwrite(&val, sizeof(double), 1, f);
	fclose(f);
}

//----------------------------------------------------------------------------
// Saves the current camera view
void LoadCamera(vtkCamera* cam, const char* path)
{
	std::string fname = CombinePathName(path, "camera", 0, ".bin");

	FILE* f = NULL;
	if (0 != fopen_s(&f, fname.c_str(), "rb"))
		return;

	double pos[3];
	fread(pos, sizeof(double), 3, f);
	cam->SetPosition(pos);

	fread(pos, sizeof(double), 3, f);
	cam->SetFocalPoint(pos);

	double val;
	fread(&val, sizeof(double), 1, f);
	cam->SetViewAngle(val);

	fread(pos, sizeof(double), 3, f);
	cam->SetViewUp(pos);

	fread(&val, sizeof(double), 1, f);
	cam->SetDistance(val);
	fclose(f);
}


//#include <vld.h>

//------------------------------------------------------------------------
/*static*/ void KeypressCallback(vtkObject* caller,
	long unsigned int vtkNotUsed(eventId), void* clientData, void* vtkNotUsed(callData))
	//------------------------------------------------------------------------
{
	vtkRenderWindowInteractor *iren = static_cast<vtkRenderWindowInteractor*>(caller);
	switch (iren->GetKeyCode())
	{
	case '1':
		vtkMapperActor::ShowHideMAs(VISUALGROUP_BONE, showBones = !showBones);
		break;
	case '2':
		vtkMapperActor::ShowHideMAs(VISUALGROUP_MUSCLE, showMuscle = !showMuscle);
		break;
	case '3':
		vtkMapperActor::ShowHideMAs(VISUALGROUP_RESULT, showFibres = !showFibres);
		break;
	case '4':
		vtkMapperActor::ShowHideMAs(VISUALGROUP_INFIBRES, showInputFibres = !showInputFibres);
		break;
	case '5':
		vtkMapperActor::ShowHideMAs(VISUALGROUP_INATTACH, showInputAttachments = !showInputAttachments);
		break;
	case '6':
		vtkMapperActor::ShowHideMAs(VISUALGROUP_NORMALS, showNormals = !showNormals);
		break;
#ifdef _DEBUG
	case '7':
		vtkMapperActor::ShowHideMAs(VISUALGROUP_LINEOFACTION, showLinesOfAction = !showLinesOfAction);
		break;
		//	case '8':
		//		vtkMapperActor::ShowHideMAs(VISUALGROUP_LINEOFACTION, showLinesOfAction = !showLinesOfAction);
		break;
#endif

	}

	iren->Render();
}

struct STATS
{
	double mean, stdev;
	double min, Q1, median, Q3, max;
	int outliers;
};

double CalcStatsDiff(const STATS& s1, const STATS& s2)
{
	return
		abs(s1.mean - s2.mean) / max(s1.mean, s2.mean) +
		abs(s1.stdev - s2.stdev) / max(s1.stdev, s2.stdev) +
		abs(s1.min - s2.min) / max(s1.min, s2.min) +
		abs(s1.median - s2.median) / max(s1.median, s2.median) +
		abs(s1.max - s2.max) / max(s1.max, s2.max) +
		abs(s1.Q1 - s2.Q1) / max(s1.Q1, s2.Q1) +
		abs(s1.Q3 - s2.Q3) / max(s1.Q3, s2.Q3);
}

#if VTK_MAJOR_VERSION < 5
#define RAD2DEG(x) (x)*vtkMath::RadiansToDgrees
#else
#define RAD2DEG(x) vtkMath::DegreesFromRadians(x)
#endif

double Percentile(const std::vector<double>& data, double p)
{
	//see https://en.wikipedia.org/wiki/Percentile
	const double C = 1.0;	//to be compatible with Excel
	double x = (data.size() + 1 - 2 * C)*p + C;

	int idxA = ((int)x);
	x = x - idxA;	//x = x % 1

	int idxB = idxA;
	--idxA;	//by default it is defined for indices starting from 1 => we need to do change	

	return x == 0 ? data[idxA] : data[idxA] + x*(data[idxB] - data[idxA]);
}

//data must be ordered, no outliers are detected
void CalcStats(const std::vector<double>& data, STATS& out)
{
	size_t n = data.size();

	out.min = data[0];
	out.max = data[n-1];
		
	out.mean = vtkLHPMuscleFibresMath::Avg(&data[0], (int)n);
	out.stdev = vtkLHPMuscleFibresMath::Dev(&data[0], (int)n);
	
	out.Q1 = Percentile(data, 0.25);
	out.median = Percentile(data, 0.5);
	out.Q3 = Percentile(data, 0.75);
}

void Rad2Degrees(STATS& stats)
{
	stats.mean = RAD2DEG(stats.mean);
	stats.stdev = RAD2DEG(stats.stdev);
	stats.min = RAD2DEG(stats.min);
	stats.Q1 = RAD2DEG(stats.Q1);
	stats.median = RAD2DEG(stats.median);
	stats.Q3 = RAD2DEG(stats.Q3);
	stats.max = RAD2DEG(stats.max);
}

//calculates the statistics without outliers
void CalcStats(vtkDoubleArray* data, STATS& out)
{
	//first convert the data into std::vector
	std::vector<double> arr(data->GetPointer(0), data->GetPointer(data->GetNumberOfValues()));
	std::sort(arr.begin(), arr.end());

	size_t lastN, newN = arr.size();	
	do
	{
		lastN = newN;

		double Q1 = Percentile(arr, 0.25);
		double Q3 = Percentile(arr, 0.75);
		double IQR = Q3 - Q1;
		double k = 1.5;	//Tukey's recommended constant to detect outliers

		double lb = Q1 - k * IQR, ub = Q3 + k*IQR;
		auto newEnd = std::remove_if(arr.begin(), arr.end(),
			[lb,ub](double v) { 
					return v < lb || v > ub; 
				}
			);

		arr.erase(newEnd, arr.end());
		newN = arr.size();
	} while (newN != lastN);

	
	CalcStats(arr, out);
	out.outliers = (int)(data->GetNumberOfValues() - newN);
}

double CalcAllStatsDiff(STATS oldstats[3], STATS newstats[3])
{
	return CalcStatsDiff(oldstats[0], newstats[0]) +
		CalcStatsDiff(oldstats[1], newstats[1]) +
		CalcStatsDiff(oldstats[2], newstats[2]);
}

//Saves the statistics into a common file
void SaveStats(const char* fileName, int iter, double diff, int SS, int VS, int IS, int nFibs, const STATS newStats[3])
{
#pragma warning(suppress: 4996)
	FILE* fOutAnalysis = fopen(fileName, "rt");
	if (fOutAnalysis != NULL) {
#pragma warning(suppress: 4996)
		fOutAnalysis = fopen(fileName, "at");
	}
	else
	{
#pragma warning(suppress: 4996)
		fOutAnalysis = fopen(fileName, "wt");
		fprintf(fOutAnalysis, "#iter\tDiff\tSS\tVS\tIS\tFibs\t\tL\tPA prox\tPA dist\n");
	}	

	fprintf(fOutAnalysis, "%d\t%.03f\t%d\t%d\t%d\t%d\toutliers\t%d\t%d\t%d\n", iter, diff, SS, VS, IS,nFibs,newStats[0].outliers, newStats[1].outliers, newStats[2].outliers);
	
	fprintf(fOutAnalysis, "\t\t\t\t\t\tmean\t%.02f\t%.02f\t%.02f\n", newStats[0].mean, newStats[1].mean, newStats[2].mean);
	fprintf(fOutAnalysis, "\t\t\t\t\t\tstdev\t%.02f\t%.02f\t%.02f\n", newStats[0].stdev, newStats[1].stdev, newStats[2].stdev);
	fprintf(fOutAnalysis, "\t\t\t\t\t\tmin\t%.02f\t%.02f\t%.02f\n", newStats[0].min, newStats[1].min, newStats[2].min);
	fprintf(fOutAnalysis, "\t\t\t\t\t\tQ1\t%.02f\t%.02f\t%.02f\n", newStats[0].Q1, newStats[1].Q1, newStats[2].Q1);
	fprintf(fOutAnalysis, "\t\t\t\t\t\tQ2\t%.02f\t%.02f\t%.02f\n", newStats[0].median, newStats[1].median, newStats[2].median);
	fprintf(fOutAnalysis, "\t\t\t\t\t\tQ3\t%.02f\t%.02f\t%.02f\n", newStats[0].Q3, newStats[1].Q3, newStats[2].Q3);
	fprintf(fOutAnalysis, "\t\t\t\t\t\tmax\t%.02f\t%.02f\t%.02f\n", newStats[0].max, newStats[1].max, newStats[2].max);

	fclose(fOutAnalysis);

}

//----------------------------------------------------------------------------// 
// Main entry point of the application
int main(int argc, const char* argv[])
//----------------------------------------------------------------------------
{
	printf(
		"Cholt Muscle Decomposition\n"
		"Copyright (c) 2009-2017 University of West Bohemia\n"
		"---------------------------------------------------\n"
		"Use:\n"
		"vtkMuscleDecomposition.exe input_data_dir [[IS SS VS [PT]][RD]]\n"
		"IS = smoothing factor (default is 10)\n"
		"SS = surface reconstruction factor (default is 5)\n"
		"VS = volume reconstruction factor (default is 5)\n"
		"PT = pennation angle threshold (default is 4/9*PI = 1.3962)\n"
		"RD = required resolution in mm used to determine IS, SS, VS\n"
		"===================================================\n\n"
	);

	//vtkDataObject::SetGlobalUseDataCacheManagerFlag(1);	
	sphere = NULL;
	std::vector<bool> s;
	s.push_back(false);

	// Create usable Renderer
	ren1 = vtkRenderer::New();
	ren1->SetBackground(1, 1, 1);

	/* This is a camera that can be used so the angle and camera settings for different data stays the same
	vtkCamera* cam = vtkCamera::New();
	cam->SetViewUp(0,0,1);
	cam->SetPosition(20, -40, 60);
	cam->SetFocalPoint(0,+10,0);
	ren1->SetActiveCamera(cam);
	ren1->ResetCamera();*/

	// lightly process input params
	if (argc < 2) {
		printf("Source folder not specified!");
		return 1;
	}

	if (!DirectoryExists(argv[1]))
	{
		printf("Source directory does not exist!");
		return 1;
	}


	vtkSmartPointer<vtkFileOutputWindow> fileOutputWindow = vtkFileOutputWindow::New();
	fileOutputWindow->Delete();
	fileOutputWindow->SetFileName("output.txt");

	vtkSmartPointer<vtkOutputWindow> outputWindow = vtkOutputWindow::GetInstance();
	if (outputWindow.GetPointer() != NULL)
	{
		outputWindow->SetInstance(fileOutputWindow.GetPointer());
	}

	// Create muscle decompositor
	vtkMuscleDecomposer *musdec = vtkMuscleDecomposer::New();

	// Load all needed data
	printf("Loading data from '%s'...\n", argv[1]);
	LoadAllData(musdec, argv[1]);

	// Decomposition Parameters
	int SS = 5, VS = 5, IS = 10;
	double PT = vtkMath::Pi() * 4 / 9;

	double obbsz[3];
	vtkOBBNode nodeMuscle;
	memset(&nodeMuscle, sizeof(nodeMuscle), 0);

	vtkOBBTree::ComputeOBB(g_MuscleMesh->GetPoints(),
		nodeMuscle.Corner, nodeMuscle.Axes[0], nodeMuscle.Axes[1], nodeMuscle.Axes[2], obbsz);

	//assume that the longest dimension is along the fibres
	for (int i = 0; i < 3; i++) {
		obbsz[i] = vtkMath::Norm(nodeMuscle.Axes[i]);
	}
	std::sort(&obbsz[0], &obbsz[3]);
	double szCoeff = obbsz[0] / 16;

	sphere->SetRadius(sphere->GetRadius() * szCoeff);

	if (argc > 3)
	{
		IS = atoi(argv[2]);
		SS = atoi(argv[3]);
		VS = atoi(argv[4]);

		if (argc > 5)
		{
			PT = atof(argv[5]);
		}
	}
	else
	{
		double resD = 5; //5 mm resolution
		if (argc == 3)
			resD = atof(argv[2]);

		//automatical determination of SS, VS, and IS parameters		
		VS = (int)(obbsz[0] / resD);
		SS = (int)(obbsz[1] / resD);
		IS = (int)(obbsz[2] / (4 * resD));
	}




	musdec->SetInterpolationSubdivision(IS);
	musdec->SetSurfaceSubdivision(SS);
	musdec->SetVolumeSubdivision(VS);
	musdec->SetPennationTreshold(PT);

	// folowing are default values of the constructor, just to sum them up, no need to use them
	musdec->SetInterpolationMethod(SPLINE_CATMULL_ROM);
	musdec->SetDoOffsetFibers(true);
	musdec->SetDoFlipNormals(true);
	
	//musdec->SetDoOffsetFibers(false);
	//musdec->SetDoForceFlipNormals(true);
	musdec->SetDoSpreadFibers(false);
	//musdec->SetDoConnectFibersToAA(false);
	//musdec->SetVolumeSubdivision(0);
	//musdec->SetInterpolationSubdivision(0);

#ifdef CMPB_JOURNAL
	if (VS == 0 && SS == 0) {
		musdec->SetDoConnectFibersToAA(false);
	}
#endif
	// GO!
	cout << "s = " << SS << ", v = " << VS << ", r = " << IS << "\n";
	musdec->Update();

	//	exit(1);

		// We tube the results so the polylines look more fiber-ish
	vtkTubeFilter* tuber = vtkTubeFilter::New();
#if VTK_MAJOR_VERSION < 5
	tuber->SetInput(musdec->GetOutput(0));
#else
	tuber->SetInputConnection(musdec->GetOutputPort(0));
#endif

	tuber->SetRadius(0.7 * szCoeff);
	tuber->SetNumberOfSides(6);


#if defined (ANALYSE_FIBRES) || defined(CMPB_JOURNAL)

	//BES: 5.9.2014 - Calculate Lines of Axis	
	vtkLHPMuscleFibresResample* resampler = vtkLHPMuscleFibresResample::New();
	vtkLHPMuscleFibresAnalysisFilter* analyser = vtkLHPMuscleFibresAnalysisFilter::New();

	analyser->SetComputeProximalLineOfAction(1);
	analyser->SetComputeDistalLineOfAction(1);
//	analyser->SetProximalLineOfAction(-1, 0, 0);
//	analyser->SetDistalLineOfAction(0, 1, 0);


	analyser->SetComputeProximalPennationAngle(1);
	analyser->SetComputeDistalPennationAngle(1);
	analyser->SetComputeFibresLengths(1);

#ifdef CMPB_JOURNAL
	analyser->SetComputeFibresRadii(0);
	analyser->SetComputeTotalFibresVolume(0);
	analyser->SetComputePCSA(0);
#else
	analyser->SetComputeFibresRadii(1);
	analyser->SetComputeTotalFibresVolume(1);
	analyser->SetComputePCSA(1);
	analyser->SetFibresVolumeCalculationMethod(vtkLHPMuscleFibresAnalysisFilter::VOLUME_CALCULATION_RAVICHANDIRAN);
#endif

	for (int exp = 0; exp < 1; exp++)
	{
		if (exp == 0)
		{
#if VTK_MAJOR_VERSION < 5
			analyser->SetInput(musdec->GetOutput());
#else
			analyser->SetInputConnection(musdec->GetOutputPort());

			//analyser->SetInputData(LoadFibresData(argv[1]));

			/*std::string fileName = argv[1];
			fileName.append("\\LHDLdata.vtk");

			vtkSmartPointer<vtkPolyDataReader> readerPS = vtkSmartPointer<vtkPolyDataReader>::New();
			readerPS->SetFileName(fileName.c_str());
			analyser->SetInputConnection(readerPS->GetOutputPort());
			*/
#endif
	}
		else
		{
#if VTK_MAJOR_VERSION < 5
			resampler->SetInput(musdec->GetOutput());
			analyser->SetInput(resampler->GetOutput());
#else
			resampler->SetInputConnection(musdec->GetOutputPort());
			analyser->SetInputConnection(resampler->GetOutputPort());
#endif
}


#ifdef CMPB_JOURNAL
		double statsDiff = DBL_MAX;
		STATS statsB1[3], statsB2[3];		
		memset(&statsB1[0], 0, sizeof(statsB1));

		STATS* oldStats = statsB1;
		STATS* newStats = statsB2;

		int iter = 0;	//idx of iteration
		do
		{
#endif

			analyser->Update();

#ifdef _DEBUG
			//add line of action to visualization	
			double size = fabs(musdec->GetOutput()->GetBounds()[5] - musdec->GetOutput()->GetBounds()[4]) / 4;

			new vtkMapperActor(
				GetLineOfActionObj(analyser->GetProximalLineOfAction(),
					analyser->GetProximalLineOfActionPos(), 10 * size),
				240, 240, 0, 255, false, ren1,
				VISUALGROUP_LINEOFACTION, showLinesOfAction
			);

			/*		new vtkMapperActor(
						GetLineOfActionObj(analyser->GetProximalLineOfActionLR(), analyser->GetProximalLineOfActionPos(), size),
						92, 92, 0, 255, false, ren1, VISUALGROUP_LINEOFACTION, showLinesOfAction);
						*/
			new vtkMapperActor(
				GetLineOfActionObj(analyser->GetDistalLineOfAction(), analyser->GetDistalLineOfActionPos(), 10 * size),
				0, 240, 240, 255, false, ren1, VISUALGROUP_LINEOFACTION, showLinesOfAction);
			/*
			new vtkMapperActor(
				GetLineOfActionObj(analyser->GetDistalLineOfActionLR(), analyser->GetDistalLineOfActionPos(), size),
				0, 92, 92, 255, false, ren1, VISUALGROUP_LINEOFACTION, showLinesOfAction);
				*/
#endif

			std::string name = argv[1];
			size_t slashPos = name.find_last_of(SLASH);
			if (slashPos != std::string::npos)
				name = name.substr(slashPos + 1);

			vtkPolyData* anlOut = analyser->GetOutput();

			vtkDoubleArray* fibLengths = vtkDoubleArray::SafeDownCast(
				anlOut->GetCellData()->GetArray(vtkLHPMuscleFibresAnalysisFilter::FibresLengthsFieldName));

			vtkDoubleArray* fibPrxPA = vtkDoubleArray::SafeDownCast(
				anlOut->GetCellData()->GetArray(vtkLHPMuscleFibresAnalysisFilter::FibresProximalPAFieldName));

			vtkDoubleArray* fibDistPA = vtkDoubleArray::SafeDownCast(
				anlOut->GetCellData()->GetArray(vtkLHPMuscleFibresAnalysisFilter::FibresDistalPAFieldName));

			if (fibLengths == NULL || fibPrxPA == NULL || fibDistPA == NULL)
				printf("ERROR: analysing of output fibres failed. The output fibres are probably corrupted. ");
			else
			{
#ifdef CMPB_JOURNAL
				std::ostringstream strstream;
				strstream << "FibLengths_" << name << "_" << SS << "_" << VS << "_" << IS << ".txt";

#pragma warning(suppress: 4996)
				FILE* fOutFibLengths = fopen(strstream.str().c_str(), "wt");
				int N = fibLengths->GetNumberOfTuples();
				for (int i = 0; i < N; i++)
				{
					fprintf(fOutFibLengths, "%f\n", fibLengths->GetValue(i));
				}
				fclose(fOutFibLengths);

				strstream.clear();
				strstream << "PrxPA_" << name << "_" << SS << "_" << VS << "_" << IS << ".txt";

#pragma warning(suppress: 4996)
				fOutFibLengths = fopen(strstream.str().c_str(), "wt");
				N = fibPrxPA->GetNumberOfTuples();
				for (int i = 0; i < N; i++)
				{
					fprintf(fOutFibLengths, "%f\n", fibPrxPA->GetValue(i));
				}
				fclose(fOutFibLengths);

				strstream.clear();
				strstream << "DstPA_" << name << "_" << SS << "_" << VS << "_" << IS << ".txt";

#pragma warning(suppress: 4996)
				fOutFibLengths = fopen(strstream.str().c_str(), "wt");
				N = fibDistPA->GetNumberOfTuples();
				for (int i = 0; i < N; i++)
				{
					fprintf(fOutFibLengths, "%f\n", fibDistPA->GetValue(i));
				}
				fclose(fOutFibLengths);
#endif


				const char* fileName = "FibresAnalysis.txt";
#pragma warning(suppress: 4996)
				FILE* fOutAnalysis = fopen(fileName, "rt");
				if (fOutAnalysis != NULL) {
#pragma warning(suppress: 4996)
					fOutAnalysis = fopen(fileName, "at");
				}
				else
				{
#pragma warning(suppress: 4996)
					fOutAnalysis = fopen(fileName, "wt");
					fprintf(fOutAnalysis, "NAME\tN\tResolution\tFBL\tProximal PA\tDistal PA\tPA\tVolumeMuscle"
#ifndef CMPB_JOURNAL
						"\tVolumeLee14\tVolumeRav09\tVolumeKoh15"
#endif
						"\tPCSA VolumeMuscle"
#ifndef CMPB_JOURNAL
						"\tPCSA VolumeLee14\tPCSA VolumeRav09\tPCSA VolumeKoh15"
#endif
						"\n");
				}

				double vols[4], PCSAs[4];

				vtkMassProperties* mass = vtkMassProperties::New();
#if VTK_MAJOR_VERSION < 5
				mass->SetInput(musdec->MuscleMesh);
#else
				mass->SetInputData(musdec->MuscleMesh);
#endif
				vols[0] = mass->GetVolume();

				int nTotalPA = fibPrxPA->GetNumberOfTuples() + fibDistPA->GetNumberOfTuples();
				double* PAtotl = new double[nTotalPA];
				memcpy(PAtotl, fibPrxPA->GetPointer(0), fibPrxPA->GetNumberOfTuples() * sizeof(double));
				memcpy(PAtotl + fibPrxPA->GetNumberOfTuples(), fibDistPA->GetPointer(0), fibDistPA->GetNumberOfTuples() * sizeof(double));

				double averagePA = /*(vtkLHPMuscleFibresMath::Sum(analyser->GetFibresProximalPA())
					+ vtkLHPMuscleFibresMath::Sum(analyser->GetFibresDistalPA()))
					/ (2 * analyser->GetInput()->GetNumberOfCells());*/

					vtkLHPMuscleFibresMath::Avg(PAtotl, nTotalPA);

				PCSAs[0] = vols[0] * cos(averagePA) / vtkLHPMuscleFibresMath::Avg(analyser->GetFibresLengths());

				mass->Delete();

				vols[1] = PCSAs[1] = -1.0;
#ifndef CMPB_JOURNAL	
				vols[2] = analyser->GetTotalFibresVolume();
				PCSAs[2] = analyser->GetPCSA();

				analyser->SetFibresVolumeCalculationMethod(vtkLHPMuscleFibresAnalysisFilter::VOLUME_CALCULATION_CYLINDRICAL);
				analyser->Update();

				vols[3] = analyser->GetTotalFibresVolume();
				PCSAs[3] = analyser->GetPCSA();
#endif

				fprintf(fOutAnalysis, "%s\t%d\t%cR\t%.2f ± %.2f (%.2f - %.2f) cm\t"
					"%.2f ± %.2f (%.2f - %.2f)°\t%.2f ± %.2f (%.2f - %.2f)°\t%.2f ± %.2f (%.2f - %.2f)°\t"
					"%.2f cm3"
#ifndef CMPB_JOURNAL
					"\t%.2f cm3\t%.2f cm3\t%.2f cm3"
#endif
					"\t%.2f cm2"
#ifndef CMPB_JOURNAL
					"\t%.2f cm2\t%.2f cm2\t%.2f cm2"
#endif
					"\n",
					name.c_str(), (int)fibLengths->GetNumberOfTuples(), exp == 0 ? 'L' : 'H',
					vtkLHPMuscleFibresMath::Avg(fibLengths) / 10,
					vtkLHPMuscleFibresMath::Dev(fibLengths) / 10,
					vtkLHPMuscleFibresMath::Min(fibLengths) / 10,
					vtkLHPMuscleFibresMath::Max(fibLengths) / 10,

					RAD2DEG(vtkLHPMuscleFibresMath::Avg(fibPrxPA)), RAD2DEG(vtkLHPMuscleFibresMath::Dev(fibPrxPA)), RAD2DEG(vtkLHPMuscleFibresMath::Min(fibPrxPA)), RAD2DEG(vtkLHPMuscleFibresMath::Max(fibPrxPA)),
					RAD2DEG(vtkLHPMuscleFibresMath::Avg(fibDistPA)), RAD2DEG(vtkLHPMuscleFibresMath::Dev(fibDistPA)), RAD2DEG(vtkLHPMuscleFibresMath::Min(fibDistPA)), RAD2DEG(vtkLHPMuscleFibresMath::Max(fibDistPA)),

					RAD2DEG(vtkLHPMuscleFibresMath::Avg(PAtotl, nTotalPA)), RAD2DEG(vtkLHPMuscleFibresMath::Dev(PAtotl, nTotalPA)),
					RAD2DEG(vtkLHPMuscleFibresMath::Min(PAtotl, nTotalPA)), RAD2DEG(vtkLHPMuscleFibresMath::Max(PAtotl, nTotalPA)),

					vols[0] / 1000,
#ifndef CMPB_JOURNAL
					vols[1] / 1000, vols[2] / 1000, vols[3] / 1000,
#endif
					PCSAs[0] / 100
#ifndef CMPB_JOURNAL
					, PCSAs[1] / 100, PCSAs[2] / 100, PCSAs[3] / 100
#endif
				);


				fclose(fOutAnalysis);

				delete[] PAtotl;
			}

#ifdef CMPB_JOURNAL
			CalcStats(fibLengths, newStats[0]);
			CalcStats(fibPrxPA, newStats[1]); Rad2Degrees(newStats[1]);
			CalcStats(fibDistPA, newStats[2]); Rad2Degrees(newStats[2]);
			statsDiff = CalcAllStatsDiff(oldStats, newStats);
	
			SaveStats((name + "_STATS.txt").c_str(), iter, statsDiff, 				
				SS, VS, IS, fibLengths->GetNumberOfValues(), newStats);

			std::swap(oldStats, newStats);

			iter++;

			musdec->SetSurfaceSubdivision(++SS);
			musdec->SetVolumeSubdivision(++VS); 
			musdec->SetDoConnectFibersToAA(true);
			analyser->Update();						
		} while (statsDiff > 1);	//end
#endif

	} //end for


	vtkPolyData* anlOut = analyser->GetOutput();
	anlOut->GetPointData()->SetActiveScalars(
		vtkLHPMuscleFibresAnalysisFilter::FibresRadiiFieldName
	);

#if VTK_MAJOR_VERSION < 5
	tuber->SetInput(anlOut);
#else
	tuber->SetInputData(anlOut);
#endif

	//tuber->SetVaryRadiusToVaryRadiusByAbsoluteScalar();


	analyser->Delete();
	resampler->Delete();
#endif

	tuber->Update();	//produce Output

	tuber->GetOutput()->GetPointData()->SetActiveScalars(TENDON_PART_STRING);
	vtkMapperActor* resultMA = new vtkMapperActor(tuber->GetOutput(),
		255, 30, 30, 255, false, ren1,
		VISUALGROUP_RESULT, showFibres);

	double scalarsLUTintensity = 1.0;

	vtkActor* resultActor = resultMA->GetActor();
	vtkPolyDataMapper* mapper = vtkPolyDataMapper::SafeDownCast(resultActor->GetMapper());
	mapper->SetScalarVisibility(1);
	mapper->SetColorModeToMapScalars();

	vtkLookupTable* lut = vtkLookupTable::New();
	lut->SetNumberOfTableValues(2);
	lut->SetTableRange(0.0, 1.0);

	lut->SetTableValue(0, 1.0*scalarsLUTintensity, 0.118*scalarsLUTintensity, 0.118*scalarsLUTintensity);
	lut->SetTableValue(1, 0.95*scalarsLUTintensity, 0.95*scalarsLUTintensity, 0.95*scalarsLUTintensity);

	mapper->SetLookupTable(lut);
	lut->Delete();


	// this was for the images of input data - highlights the edges of the muscle mesh
	//vtkMapperActor* muscleMAw = new vtkMapperActor(musdec->MuscleMesh, 0, 0, 0, 128, true, ren1); 
	vtkMapperActor* muscleMA = new vtkMapperActor(musdec->MuscleMesh,
		(int)(255 * 0.9599999785423279), (int)(255 * 0.2899999916553497), (int)(255 * 0.2899999916553497), 128, false, ren1,
		VISUALGROUP_MUSCLE, showMuscle);

/*
	vtkMapperActor* muscleMAW = new vtkMapperActor(musdec->MuscleMesh, 128, 128, 128, 192, true, ren1,
		VISUALGROUP_MUSCLE, showMuscle);
*/

	vtkMapperActor* normalsMA = new vtkMapperActor(musdec->NormalMesh, 56, 84, 0, 255, false, ren1,
		VISUALGROUP_NORMALS, showNormals);



	// Create RenderWindow and add the renderer
	// Create trackball style interactor, assign to renderwindow
	renWin = vtkRenderWindow::New();
	vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
	vtkInteractorStyleTrackballCamera *style = vtkInteractorStyleTrackballCamera::New();
	iren->SetRenderWindow(renWin);	
	renWin->SetSize(1300, 800);
	renWin->AddRenderer(ren1);
	iren->SetInteractorStyle(style);

	// save results
	SaveAllData(musdec, argv[1]);

	vtkCallbackCommand *keypressCallback = vtkCallbackCommand::New();
	keypressCallback->SetClientData(NULL);
	keypressCallback->SetCallback(&KeypressCallback);
	iren->AddObserver(vtkCommand::KeyPressEvent, keypressCallback);

	// Loop on the interactor
	iren->Initialize();

	vtkCamera* cam = ren1->GetActiveCamera();
	ren1->ResetCamera();

	LoadCamera(cam, argv[1]);
	ren1->ResetCameraClippingRange();

	iren->Start();
	
	SaveCamera(ren1->GetActiveCamera(), argv[1]);


	// Free up all objects we created.
	vtkMapperActor::DeleteAllMAs();

	musdec->Delete();
	tuber->Delete();
	ren1->Delete();
	renWin->Delete();
	iren->Delete();
	style->Delete();


	//fix VTK bug
	vtkTimerLog::CleanupLog();

	// Dump memleaks to output window. Only leaks with path specified in the dump are mine (others are VTK)
	//_CrtDumpMemoryLeaks();
	return 0;
}
