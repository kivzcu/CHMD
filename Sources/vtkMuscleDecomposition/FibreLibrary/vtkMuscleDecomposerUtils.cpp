#include "vtkMuscleDecomposerUtils.h"

//// =========================================================================
//// Utils & Math
//// =========================================================================

// Computes Line-Point distance between line defined by linepoint1 and linepoint2, and point point
double vtkMuscleDecomposerUtils::LinePointDistance(double* linepoint1, double* linepoint2, double* point)
{
	// (|(point-linepoint1)x(point-linepoint2)|)/(|linepoint2-linepoint1|)
	double cross[3];
	double sub1[3], sub2[3], sub3[3];
	Subtract(point, linepoint1, sub1);
	Subtract(point, linepoint2, sub2);
	Subtract(linepoint2, linepoint1, sub3);
	vtkMath::Cross(sub1, sub2, cross);
	
	return vtkMath::Norm(cross) / vtkMath::Norm(sub3);
}

// Computes Line-Line distance between lines defined by (p11, p12) and (p21, p22)
double vtkMuscleDecomposerUtils::LineLineDistance(double* p11, double* p12, double* p21, double* p22)
{
	/*
	The distance between two skew lines with equations
	x	=	x_1+(x_2-x_1)s
	x	=	x_3+(x_4-x_3)t

	D=(|c·(axb)|)/(|axb|)

	by defining
	a	=	x_2-x_1
	b	=	x_4-x_3
	c   =	x_3-x_1
	*/
	double a[3];
	double b[3];
	double c[3];
	double aCrossB[3];

	Subtract(p12, p11, a);
	Subtract(p22, p21, b);
	Subtract(p21, p11, c);

	
	vtkMath::Cross(a, b, aCrossB);

	double result = abs(vtkMath::Dot(c, aCrossB) / sqrt(aCrossB[0]*aCrossB[0] + aCrossB[1]*aCrossB[1] + aCrossB[2]*aCrossB[2]));
	return result;
}

// Subtracts point1 from point2 to result
void vtkMuscleDecomposerUtils::Subtract(double* point1, double* point2, double* result)
{
	for (int i = 0; i < 3; i++)
	{
		result[i] = point2[i] - point1[i];
	}
}

// Computes result[i] = point2[i] * multiplication + point1[i];
void vtkMuscleDecomposerUtils::Add(double* point1, double* point2, double* result, double multiplication)
{
	for (int i = 0; i < 3; i++)
	{
		result[i] = point2[i] * multiplication + point1[i];
	}
}

double vtkMuscleDecomposerUtils::PointPointDistance(double* point1, double* point2)
{
	return sqrt(vtkMath::Distance2BetweenPoints(point1, point2));
}

// Computes normal of from 3 points (point1 being the normal origin)
void vtkMuscleDecomposerUtils::ComputeNormal(double* point1, double* point2, double* point3, double* normal)
{
	double u[3];
	double v[3];
	double result[3];
	Subtract(point2, point1, u);
	Subtract(point3, point1, v);
	vtkMath::Cross(u, v, result);
	vtkMath::Normalize(result);
	normal[0] = result[0];
	normal[1] = result[1];
	normal[2] = result[2];
}

void vtkMuscleDecomposerUtils::ApplyBarycentric(double* point1, double* point2, double* point3, double* barycentric, double* result)
{
	result[0] = point1[0] * barycentric[0] + point2[0] * barycentric[1] + point3[0] * barycentric[2];
	result[1] = point1[1] * barycentric[0] + point2[1] * barycentric[1] + point3[1] * barycentric[2];
	result[2] = point1[2] * barycentric[0] + point2[2] * barycentric[1] + point3[2] * barycentric[2];
}

double vtkMuscleDecomposerUtils::ComputeBendAngle(double* p1, double* p2, double* p3)
{
	double vector1[3], vector2[3];
	double dot;

	Subtract(p1, p2, vector1);
	Subtract(p2, p3, vector2);

	vtkMath::Normalize(vector1);
	vtkMath::Normalize(vector2);
	dot = vtkMath::Dot(vector1, vector2);
	return acos(dot);
}