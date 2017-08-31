#pragma once

#include "vtkMath.h"

class vtkMuscleDecomposerUtils
{
public:
	
	// Computes Line-Point distance between line defined by linepoint1 and linepoint2, and point point
	static double LinePointDistance(double* linepoint1, double* linepoint2, double* point);

	// Computes Line-Line distance between lines defined by (p11, p12) and (p21, p22)
	static double LineLineDistance(double* p11, double* p12, double* p21, double* p22);

	// Subtracts point1 from point2 to result
	static void Subtract(double* point1, double* point2, double* result);

	/// Computes result[i] = point2[i] * multiplication + point1[i]; 
	static void Add(double* point1, double* point2, double* result, double multiplication);

	static double PointPointDistance(double* point1, double* point2);

	// Computes normal of from 3 points (point1 being the normal origin)
	static void ComputeNormal(double* point1, double* point2, double* point3, double* normal);

	static void ApplyBarycentric(double* point1, double* point2, double* point3, double* barycentric, double* result);

	static double ComputeBendAngle(double* p1, double* p2, double* p3);
};

