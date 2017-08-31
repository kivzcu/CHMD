#pragma once

#ifndef __vtkLHPMuscleFibresMath_h_
#define __vtkLHPMuscleFibresMath_h_

#include <float.h>

class vtkPolyData;
class vtkIntArray;
class vtkDoubleArray;

/** This is a helper class containing the key routines for the analysis of muscle fibres. 
 
 The requirements are: 1) Fibres are represented by poly-lines,  2) One poly-line is stored as one cell,
 i.e., the cell data structure of the polydata object is:

 M_0, i_{0,0}, i_{0,1}, ... i_{0,M0-1},
 M_1, i_{1,0}, i_{1,1}, ... i_{1,M1-1},
 ...
 M_{N-1}, i_{N-1,0}, i_{N-1,1}, ... i_{N-1,M1-1},

 where N is the number fibres, M_f is the number of vertices forming the poly-line representing the f-th fibre,
 i_{f, v} is the index of v-th vertex of f-th fibre,

 3) No vertex is be used repeatedly, and 4) every vertex belongs to some fibre.

 If any of these requirements are not satisfied, the results are unpredictable and the application may even crash.
 */
class vtkLHPMuscleFibresMath
{
public:
	typedef double VCoord[3];

public:
	/** Counts the segments of every fibre in poly and stores this information into out */
	static void CountFibresSegments(const vtkPolyData* poly, vtkIntArray* out);

	/** Counts the segments of every fibre in poly and stores this information into pnVerts.
	N.B. the buffer must be big enough to hold poly->GetNumberOfCells() items. */
	static void CountFibresSegments(const vtkPolyData* poly, int* pnVerts);

	/** Computes the lengths of each segment of the fibres in poly and stores them into out. */
	static void CalculateFibresSegmentsLengths(const vtkPolyData* poly, vtkDoubleArray* out);

	/** Computes the lengths of each segment of the fibres in poly and stores them into pdblLen.
	N.B. the buffer must be big enough to hold poly->GetNumberOfPoints() - GetNumberOfCells() items. */
	static void CalculateFibresSegmentsLengths(const vtkPolyData* poly, double* pdblLen);

	/** Computes the overall lengths of input fibres and stores them into out.	*/
	static void CalculateFibresLengths(const vtkIntArray* segmentsCount, const vtkDoubleArray* segmentsLengths, vtkDoubleArray* out);

	/** Computes the overall lengths of input fibres and stores them into pdblTotLen.
	N.B. the output buffer must be capable to hold :nFibres: items. */
	static void CalculateFibresLengths(int nFibres, const int* pdblSegCount,
		const double* pdblSegLen, double* pdblTotLen);

	/** Computes either chordal or centripetal parameterization of the input fibres and stores it into out. 
	Note: the parameterization of a fibre is from 0 to t_m-1. */
	static void CalculateFibresParameterization(const vtkPolyData* poly, const vtkDoubleArray* segmentsLengths, 
		bool centripetal, vtkDoubleArray* out);

	/** Computes either chordal or centripetal parameterization of the input fibres and stores it into out. 	
	N.B. the output buffer must be capable holding nPoints (poly->GetNumberOfPoints()) parameters. */
	static void CalculateFibresParameterization(const vtkPolyData* poly, 
		const double* pdblSegLen, bool centripetal, double* pPars);
	
	/** Normalizes the given parameterization of fibres into 0 - 1. 
	It is safe to call this method with out being parms for in-place normalization. */
	static void NormalizeFibresParameterization(const vtkPolyData* poly,
		const vtkDoubleArray* params, vtkDoubleArray* out);

	/** Normalizes the given parameterization of fibres into 0 - 1.
	It is safe to call this method with out being parms for in-place normalization. 
	N.B. he output buffer must be capable holding nPoints (poly->GetNumberOfPoints()) parameters. */
	static void NormalizeFibresParameterization(const vtkPolyData* poly,
		const double* params, double* out);
	
	/** Gets the minimum*/
	static double Min(const vtkDoubleArray* inArray);

	/** Gets the minimum*/
	inline static double Min(const double* values, int count) {
		double ret = DBL_MAX;
		for (int i = 0; i < count; i++) {
			if (values[i] < ret)
				ret = values[i];
		}

		return ret;
	}

	/** Gets the minimum*/
	static double Max(const vtkDoubleArray* inArray);

	/** Gets the minimum*/
	inline static double Max(const double* values, int count) {
		double ret = -DBL_MAX;
		for (int i = 0; i < count; i++) {
			if (values[i] > ret)
				ret = values[i];
		}

		return ret;
	}

	/** Sums all values in the input array*/
	static double Sum(const vtkDoubleArray* inArray);

	/** Sums all values in the input array*/
	inline static double Sum(const double* values, int count) {
		double sum = 0.0;
		for (int i = 0; i < count; i++) {
			sum += values[i];
		}

		return sum;
	}
	
	/** Sums all values in the input array*/
	static double Avg(const vtkDoubleArray* inArray);

	/** Sums all values in the input array*/
	inline static double Avg(const double* values, int count) {
		return Sum(values, count) / count;
	}

	/** Get the (corrected sample) standard deviation*/
	static double Dev(const vtkDoubleArray* inArray);

	/** Get the (corrected sample) standard deviation*/
	static double Dev(const double* values, int count);
};

#endif //__vtkLHPMuscleFibresMath_h_
