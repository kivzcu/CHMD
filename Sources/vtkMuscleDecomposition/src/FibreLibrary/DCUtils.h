#ifndef DCUtils_h
#define DCUtils_h

#include "vtkMath.h"

class DCUtils
{
public:
	DCUtils(void);
	~DCUtils(void);
	// Computes Line-Line distance between lines defined by (p11, p12) and (p21, p22)
	static double LineLineDistance(double* p11, double* p12, double* p21, double* p22);
	// Substracts point1 from point2 to result
	static void Substract(double* point1, double* point2, double* result);
	static void Add(double* point1, double* point2, double* result, double multiplication);
	static double PointPointDistance(double* point1, double* point2);
	// Computes normal of from 3 points (point1 being the normal origin)
	static void ComputeNormal(double* point1, double* point2, double* point3, double* normal);
	static int PointIsWithinBounds(double point[3], double bounds[6], double delta[3]);
};

#endif