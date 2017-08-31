/*=========================================================================
Module:    $RCSfile: vtkMuscleDecomposer.h $
Language:  C++
Date:      $Date: 2013-05-01 17:00 $
Version:   $Revision: 0.1.0.0 $
Author:    David Cholt
Notes:	   Muscle decomposition method using cadaver fibres. The method uses
		   input surface muscle mesh (AddMesh()), input attachment data
		   (AddAttachmentData()) and cadaver fibers (AddFibers()) to decompose
		   the muscle into a set of polylines. Note that the fibers are 
		   crucial for the decomposition process, plase see my master's thesis
		   (in czech) or contact me at cholt@students.zcu.cz for more info.
		   Also please see demoApp.cpp for class usage.
=========================================================================*/

#pragma once
//----------------------------------------------------------------------------
// Include:
//----------------------------------------------------------------------------
#include "vtkObject.h"
#if VTK_MAJOR_VERSION < 5
#include "vtkPolyDataToPolyDataFilter.h"
#define vtkPolyDataAlgorithm vtkPolyDataToPolyDataFilter
#else
#include "vtkPolyDataAlgorithm.h"
#endif

#include "vtkMath.h"
#include "vtkPolyData.h"
#include "vtkCellArray.h"
#include "vtkTriangle.h"
#include "vtkOBBTree.h"
#include "vtkTubeFilter.h"
#include "vtkMAFPolyDataCutOutFilterEx.h"
#include "vtkSmartPointer.h"
#include "vtkSTLWriter.h"

#include <float.h>
#include <stdlib.h>     /* srand, rand */
#include <time.h>		/* testing */
#include <fstream>
#include <iomanip>
#include <string>
#include <map>
#include <vector>
#include <algorithm>

#include "hungarian.h"
#include "vtkMuscleDecomposerUtils.h"


//----------------------------------------------------------------------------
// Define:
//----------------------------------------------------------------------------

#define ATTACHMENT_START 0
#define ATTACHMENT_END 1

#define CONNECT_CLOSEST_TENDON_POINT 0
#define CONNECT_INTENDED_TENDON_POINT 1
#define CONNECT_TENDONS_BY_DIRECTION 2
#define CONNECT_ATTACHMENTS_BY_DIRECTION 3

#define RAY_POINT_DISTANCE 10000
#define MIN_THICKNESS 1.0

#define BINARY_SEARCH_ITERATIONS 20
#define SEARCH_THICKNESS_FACTOR 0.8;

#define SPLINE_NATURAL 0
#define SPLINE_PARABOLIC_RUNOUT 1
#define SPLINE_CUBIC_RUNOUT 2
#define SPLINE_CONSTRAINED 3
#define SPLINE_CATMULL_ROM 4

class VTK_EXPORT vtkMuscleDecomposer : public vtkPolyDataAlgorithm
{
public:
#pragma region Nested data structure classes
	
	/* Represents one fiber in dataset. */
	class Fiber
	{
		friend class vtkMuscleDecomposerInterpolator;
	public:
		std::vector<int>  PointIDs; // Points which make the fiber
		double DLayer; // Layer level
		Fiber(void);
		~Fiber(void);
	};

	
	/* Represents one point in dataset. */
	class Point
	{
	public:
		friend class vtkMuscleDecomposerInterpolator;
		double	DCoord[3];				// Coord
		double  DNormal[3];				// Normal
		bool	Interpolated;				// Is this point interpolated?

		int		Id;						// ID
		int		AttachmentConnectionID; // ID of attachment point

		// for interpolation purposes:
		double	a[3], b[3], c[3], d[3]; // (3 axes)
		double	DistToNext;				// distance to the next point
		double  t;						// parameter value in this point

		// for catmull-rom
		int prevID;
		int nextID;
		int nextNextID;

		// for sampling
		double Thickness;

		// for spreading
		int LeftAnchorID; 
		int RightAnchorID;
		bool Splits;

		//muscle and tendon fibre info
		int MuscleTendonInfo;	//-1 = NOT AVAILABLE, 0 = MUSCLE, 1 = TENDON


		Point(double x, double y, double z);
		~Point();

		// Computes coordinates of point on interpolated curve, using given parameter. 
		// Substracting the segment start t from point's t is neccessary for some 
		// interpolation methods.
		void ComputeInterpolatedCurveCoord(double t, double* result, bool Subtract);
		// Assigns normal to this Point
		void SetNormal(double* normal);
	};

#pragma endregion

#pragma region Nested interpolator
	class SplineInterpolator; // forward declaration
#pragma endregion

#pragma region MuscleDecompositor
		///
		///  Public structures and methods
		///
public: 
		vtkPolyData *MuscleMesh;
		vtkPolyData *NormalMesh;

		static vtkMuscleDecomposer* New();
		vtkTypeMacro(vtkMuscleDecomposer, vtkPolyDataAlgorithm);

		/* Adds single control fiber polyline data */
		void AddFiberData(vtkPolyData *Fibers);
		/* Adds attachment data. If Bone == null, then Data is expected to be prevoiusly cut attachment area */
		void AddAttachmentData(vtkPolyData *Data, vtkPolyData *Bone, int Position);
		/* Adds mesh data for the muscle to be decomposed */
		void AddMeshData(vtkPolyData *Mesh);
		
		/** Get and sets the InterpolationSubdivision parameter **/
		vtkGetMacro(InterpolationSubdivision, int);
		vtkSetMacro(InterpolationSubdivision, int);

		/** Get and sets the SurfaceSubdivision parameter **/
		vtkGetMacro(SurfaceSubdivision, int);
		vtkSetMacro(SurfaceSubdivision, int);

		/** Get and sets the VolumeSubdivision parameter **/
		vtkGetMacro(VolumeSubdivision, int);
		vtkSetMacro(VolumeSubdivision, int);

		/** Get and sets the ArtifactEpsilon parameter **/
		vtkGetMacro(ArtifactEpsilon, double);
		vtkSetMacro(ArtifactEpsilon, double);

		/** Get and sets the InterpolationMethod parameter **/
		vtkGetMacro(InterpolationMethod, int);
		vtkSetMacro(InterpolationMethod, int);

		/** Get and sets the DoFlipNormals parameter **/
		vtkGetMacro(DoFlipNormals, bool);
		vtkSetMacro(DoFlipNormals, bool);
		
		/** Get and sets the DoForceFlipNormals parameter **/
		vtkGetMacro(DoForceFlipNormals, bool);
		vtkSetMacro(DoForceFlipNormals, bool);

		/** Get and sets the DoOffsetFibers parameter **/
		vtkGetMacro(DoOffsetFibers, bool);
		vtkSetMacro(DoOffsetFibers, bool);

		/** Get and sets the DoStreadFibers parameter allowing the method to process muscles with multiple heads **/
		vtkGetMacro(DoSpreadFibers, bool);
		vtkSetMacro(DoSpreadFibers, bool);

		/** Get and sets the DoConnectFibersToAA parameter determining if the fibres 
		should be automatically connected to the attachment areas of the muscle **/
		vtkGetMacro(DoConnectFibersToAA, bool);
		vtkSetMacro(DoConnectFibersToAA, bool);
		
		/** Get and sets the PennationTreshold parameter **/
		vtkGetMacro(PennationTreshold, double);
		vtkSetMacro(PennationTreshold, double);

		/** Get the resulting fiber count **/
		vtkGetMacro(FiberCount, int);

		/** Gets the pointer to the 2D (START and END) array of attachment surfaces (surfaces are int the vectors) */
		vtkGetMacro(AttachmentsSurfaces,  std::vector<vtkPolyData*>*);
				
		///
		///  Used variables and structures, support methods
		///
	protected:
		/** constructor */
		vtkMuscleDecomposer();
		/** destructor */
		~vtkMuscleDecomposer();

		std::vector<vtkPolyData*>	 ControlFibers;  // list of fibers of the muscle

		std::vector<Fiber*>			 Fibers;         // list of fibers of the muscle
		std::vector<Fiber*>			 BottomFibers;   // list of bottom fibers of the muscle
		std::vector<Fiber*>			 LeftFibers;	 // list of left fibers
		std::vector<Fiber*>			 RightFibers;	 // list of right fibers
		std::vector<Fiber*>			 OutputFibers;   // list of output fibers

		std::vector<Fiber*>			 AllFibers;	     // list of all fibres

		std::vector<Point*>			 PointCloud;	 // cloud of fiber points of the muscle

		std::vector<int>				 AttachmentPoints[2];  // an array of attachment points ready to be connected
		std::vector<int>				 OpenFibers[2];		   // an array of unconnected fibers
		std::vector<int>				 OpenFibersSeconds[2]; // an array of neighbors of unconnected fibers

		std::vector<vtkPolyData*>	 AttachmentsSurfaces[2]; // attachment surface models

		// PROPERTIES
		int InterpolationMethod;
		int InterpolationSubdivision;
		int SurfaceSubdivision;
		int VolumeSubdivision;
		double ArtifactEpsilon;
		bool DoFlipNormals;
		bool DoForceFlipNormals;
		bool DoOffsetFibers;
		bool DoSpreadFibers;
		bool DoConnectFibersToAA;
		double PennationTreshold;
				
		int FiberCount;

		clock_t start; // benchmark
		std::vector<clock_t> times;
		
		vtkOBBTree* tree;			// for sampling
		bool RebuildLocator;
		///
		///  Main Filter execution methods
		///
	protected:
		vtkPolyData *OutputMesh;

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

		/** Build internal data structure **/
		bool InitData();
		/** Build internal mesh relationships **/
		void BuildMesh();
		/** Hand over finished mesh **/
		void DoneMesh();
		/** Clears the mesh from memory **/
		void ClearMesh();

		///
		///  Methods used by decomposition (ordered by use where possible)
		///
	protected:
		/** Initialises internal structures for input control fibres **/
		void InitFibers();

		/** Fix the control fibres so that they contain the same number of segments.
		The fixed fibres are stored in fixedFibres (to be deleted by the caller)*/
		void FixFibers(std::vector<vtkPolyData*>& fixedFibres);
				
		/** Connects fibers to tendons based on proximity **/
		void PrepareTendons(int targetFiberCount);
		/** Cuts off a attachment area from a surface using vtkMAFPolyDataCutOutFilterEx class**/
		void GetAttachmentAreaSurface(const vtkPoints* projpts, const vtkPolyData* surface, vtkPolyData* output);

		/** Removes too short segments from fibers /unused/ **/
		void RemoveArtifacts(double epsilon, std::vector<Fiber*> & Data);

		/** Reconstructs fibers on the muscle surface**/
		void InterpolateSurface(int subdivision, int interpolationType, std::vector<Fiber*> & Data);

		/** Computes normals for all Fibres points **/
		void ComputeNormals(void);
		/** Flips normals for all Fibres points **/
		void FlipNormals(void);
	
		/** Moves the muscle mesh so its centered in the fibre-defined muscle hull /not used/ **/
		void AlignMesh();
		
		// Computes muscle thinckness for fiber points
		void ComputeThickness(void);
		// Samples a single muscle thickness in the negative and positive (unless oneSideOnly is true) direction of normal
		double SampleThickness(double* coord, double* normal, bool oneSideOnly = false);
		// Searches for valid muscle thickness for left and right part of the segment, using binary search
		void SearchValidThickness(int LeftAnchorID, int RightAnchorID, double holeT, double* result, int iterations);

		/** Builds bottom layer of fibers from normals and thickness
		\param t parameter <0..1> denoting the layer to build
		\param target output buffer in which the built fibres are stored
		\param step distance between two layers in parameter t used to determine the maximal random shifting along the normal  
		\param tScalePerFibre non-uniform correction scaling of t, each of them >= 1, if the corrected t > 1, the fibre is not built*/
		void BuildBottomFibers(double t, std::vector<Fiber*>  & target, double step, std::vector<double>& tScalePerFibre);

		// Offsets bottom fibres to prevent the rake effect
		void OffsetFibers(std::vector<Fiber*>  & target);

		// Builds side layers after bottom layer has been built
		void BuildSideLayers(int subdivision, int interpolationType, std::vector<Fiber*>  & topFibers, std::vector<Fiber*>  & bottomFibers);

		/** Connects open ends of fibres to attachment points **/
		void ConnectAttachments();
		/** Registers fiber's open ends in appripriate list **/ 
		void RegisterOpenEnds(Fiber* fiber);
		/** Unregisters fiber's open ends in appripriate list **/ 
		void UnRegisterOpenEnds(Fiber* fiber);

		/** Gathers all fibres to AllFibres list **/
		void GatherAllFibers();
		/** Checks the attachment connection for sharp angles **/
		bool CheckConnection(int side);

		/** Shortens all fibres by how deeply the fibers are in the muscle (layer by layer)**/
		void ShortenFibersByLayers(double percent, int side);
		/** Cuts single fiber by amount percent **/
		void CutFiber(Fiber* fiber, int side, double amount);
		/** Goes trough a fibre and measures its length in midpoints**/
		void ComputeFiberLength(Fiber* fiber, int side);
		
		/** Builds OutputMesh **/
		void BuildAndOutputConnectedData(std::vector<Fiber*> & Data, char* what, bool connect);
		
		/** Interpolates (lengthwise) all fibres in Data. **/
		void InterpolateFibers(int subdivision, int interpolationType, std::vector<Fiber*> & Data, bool tendons);
		/** Interpolates (lengthwise) single fiber in Data. Initialised SplineInterpolator is needed. Ignores skipOnEnds segments on both sides. If insertSingularity, inserts even duplicate points **/
		void InterpolateFiber(Fiber* fiber, int subdivision, SplineInterpolator *interpolator, int skipOnEnds, bool insertSingularity);
				
		/** Outputs a Fiber set. */
		void OutputFibersData( std::vector<Fiber*> & Data, vtkPoints * PointIDs, vtkCellArray * lines, vtkUnsignedCharArray* mtInfo);

		/** Outputs normals for debugging purposes */
		void OutputNormals();

		/** Returns the muscle-tendon info for the point with parameter t from  <0, 1> lying between points identified by ID and nexID identifiers. */
		int InterpolateMuscleTendonInfo(int ID, int nextID, double t);

		///
		/// Console logging methods
		///
	private:
		static void Log(char* what)
		{
			printf(what);
			printf("\n");
		}
		static void Log(char* what, int parameter)
		{
			printf(what);
			printf("%d", parameter);
			printf("\n");
		}

		static void Log(char* what, double parameter)
		{
			printf(what);
			printf("%f", parameter);
			printf("\n");
		}

		static void LogFormated(char* what, int parameter)
		{
			printf(what, parameter);
			printf("\n");
		}

		static void LogFormated(char* what, char* parameter)
		{
			printf(what, parameter);
			printf("\n");
		}
		
		///
		/// File logging methods
		///
	public:
		static void LogTimes(std::vector<clock_t> times)
		{
			ifstream t;
			ofstream out;
			t.open("Times.csv", ios::in);
			out.open("Times.csv", ios::out | ios::app);
			if (t.is_open())
			{
				t.close();
			}
			else
			{
				out << setw(8) << "Triangles;" << setw(8) << "Preparation;" << setw(8) << "Surface;" << setw(8) << "Normals;" << setw(8) 
					//<< "Sampling;"<< setw(8) 
					<< "Sampling:OBBTree;" << setw(8)
					<< "Sampling:CheckNormalsOri;" << setw(8)
					<< "Sampling:RayCasting;" << setw(8)
					<< "Sampling:Spreading;" << setw(8)					
					<< "Volume;"  << setw(8) << "Connection;" <<  setw(8) << "Total;" << endl;
			}
			for (int i = 0; i < (int)times.size(); i++)
			{
				out << setw(8) << times[i] <<";";
			}
			out << endl;
			out.close();
		}
#pragma endregion

#pragma region Utils
		/** Computes distance between two PointIDs **/
		double ComputePointDistance(int ID1, int ID2);

		static char* GetSplineType(int SplineType)
		{
			switch (SplineType)
			{
			case SPLINE_NATURAL:
				return "natural runout";
				break;
			case SPLINE_PARABOLIC_RUNOUT:
				return "parabolic runout";
				break;
			case SPLINE_CUBIC_RUNOUT:
				return "cubic runout";
				break;
			case SPLINE_CONSTRAINED:
				return "constrained spline";
				break;
			case SPLINE_CATMULL_ROM:
				return "Catmull-Rom spline";
				break;
			default:
				return "ERROR";
			}
		}
#pragma endregion
};

#if VTK_MAJOR_VERSION < 5
#undef vtkPolyDataAlgorithm
#endif