#pragma once

#ifndef __vtkLHPMuscleFibresResample_h
#define __vtkLHPMuscleFibresResample_h

#include "vtkObject.h"
#if VTK_MAJOR_VERSION < 5
#include "vtkPolyDataToPolyDataFilter.h"
#define vtkPolyDataAlgorithm vtkPolyDataToPolyDataFilter
#else
#include "vtkPolyDataAlgorithm.h"
#endif


class vtkPoints;
class vtkCellArray;

/** Performs resampling of the input muscle fibres. */
class vtkLHPMuscleFibresResample : public vtkPolyDataAlgorithm
{
public:
	vtkTypeMacro(vtkLHPMuscleFibresResample, vtkPolyDataAlgorithm);

	/** Constructs a new analysis filter that calculates nothing by the default. */
	static vtkLHPMuscleFibresResample *New();

	/** Gets/Sets the precision of the output fibres, i.e., the maximal distance between two points 
	N.B. this parameter is used only if PrecisionMode is ON */
	vtkSetMacro(Precision, double);
	vtkGetMacro(Precision, double);

	/** Gets/Sets the requested number of vertices of the output fibres
	If the value is set to 0, all the fibres are resampled to contain the same number of vertices
	N.B. this parameter is used only if PrecisionMode is OFF */
	vtkSetMacro(Resolution, int);
	vtkGetMacro(Resolution, int);

	/** Turn on/off the using of centripal Catmull-Rom parameterization of fibres
	If off, the default chordal Catmull-Rom parameterization is used. */
	vtkBooleanMacro(UseCentripetalCatmullRom, int);
	vtkSetMacro(UseCentripetalCatmullRom, int);
	vtkGetMacro(UseCentripetalCatmullRom, int);

	/** Turn on/off the mode that resamples the fibres according to the specified Precision (OFF by default)
	If off, the Resolution is used. */
	vtkBooleanMacro(PrecisionMode, int);
	vtkSetMacro(PrecisionMode, int);
	vtkGetMacro(PrecisionMode, int);

protected:
	vtkLHPMuscleFibresResample();
	~vtkLHPMuscleFibresResample();

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

	/** Performs the resampling of a fibre defined by its sequence of points.
	\param coords the coordinates of the points defining the fibre (pPts); coords->GetPoint(pPts[0]) are the coordinates of the first point
	\param params the knot values, i.e., the Catmull-Rom parameterization of every point; params(pPts[0]) is the knot value of the first point
	\param segLenghts the lengths between ever pair of adjacent points of the fibre
	\params points the output buffer for points coordinates
	\param cells the output buffer for the connectivity of the produced resampled fibre 
	\param resolution the preferred resolution (to be ignored if PrecisionMode is ON)*/
	virtual void ResampleFibre(vtkIdType nPts, const vtkIdType* pPts, const vtkPoints* coords,
		const double* params, const double* segLengths, vtkPoints* points, vtkCellArray* cells, int resolution);

	/** Performs the sampling of the curve
	\param P the points P0..P3 defining the curve
	\param t the knots values t0..t3 at P0..P3
	\param iv defines the permutation of the points, e.g., {0, 1, 2, 3} means that the order of points is P0, P1, P2, P3
	\param tp the knot values from the interval <t1, t2> determining where the samples should be created
	\param count the number of samples to create
	\param points the output buffer for points coordinates
	\param idsList the output buffer for indices of the constructed points*/
	virtual void SampleCurve(double P[4][3], double t[4], int iv[4],
		double* tp, int count, vtkPoints* points, vtkIdList* idsList);

	/** Copies the fibre defined by its sequence of points to the output.
	\param coords the coordinates of the points defining the fibre (pPts); coords->GetPoint(pPts[0]) are the coordinates of the first point	
	\params points the output buffer for points coordinates
	\param cells the output buffer for the connectivity of the produced resampled fibre */
	void CopyFibre(vtkIdType nPts, const vtkIdType* pPts, const vtkPoints* coords, vtkPoints* points, vtkCellArray* cells);
protected:
	double Precision;					///<Precision of the output data, i.e., the maximal distance between two points
	int Resolution;						///<Resolution of the output data, i.e., the number of vertices per fibre
	int UseCentripetalCatmullRom;		///<if non-zero, the centripetal Catmull-Rom parameterization will be used, otherwise chordal parameterization will be used	
	int PrecisionMode;					///<if non-zero, Precision is used to determine the output, otherwise Resolution
};

#if VTK_MAJOR_VERSION < 5
#undef vtkPolyDataAlgorithm
#endif

#endif //vtkLHPMuscleFibresResample
