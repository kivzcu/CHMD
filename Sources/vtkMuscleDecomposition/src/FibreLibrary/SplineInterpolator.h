#ifndef __SplineInterpolator_h
#define __SplineInterpolator_h

#include "Point.h"
#include "Fiber.h"
#include "vtkstd/vector"
#include "vtkMath.h";

#define SPLINE_NATURAL 0
#define SPLINE_PARABOLIC_RUNOUT 1
#define SPLINE_CUBIC_RUNOUT 2
#define SPLINE_CONSTRAINED 3

class SplineInterpolator
{
public:
	SplineInterpolator(int splineType);
	~SplineInterpolator();

	void Interpolate(Fiber* fiber, vtkstd::vector<Point*> & points);

private:
	int splineType;

	Fiber* fiber;   // Pointer to fiber data
	vtkstd::vector<Point*> pointCloud; // Point data
	int n; // number of points

	vtkstd::vector<double> subdiagonal;			// Diagonal under main diagonal (for all axes)
	vtkstd::vector<double> superdiagonal;		// Diagonal above main diagonal (for all axes)
	vtkstd::vector<double> diagonal[3];			// Main diagonal (for 3 axes - modified by calculations)
	vtkstd::vector<double> result[3];			// Resulting vector, containing S_i..S_n (for 3 axes)
	vtkstd::vector<double> rightHand[3];		// Right hand vector (for 3 axes)

	vtkstd::vector<double> derivations[3];		// Derivations for Constrained Splines

	void BuildParameter();						// Builds independent variable (parameter) t
	void BuildEquation();						// Creates sub, sup and diagonal vectors and righthand vectors
	void SolveEquation();						// Solves equation and associates weights to points
	void ComputeWeights();						// Computes weights for all points using eq. result

	void ComputeFibreDerivations();
	void ComputeConstrainedWeights();
public:
	int GetInterpolationType();
};

#endif