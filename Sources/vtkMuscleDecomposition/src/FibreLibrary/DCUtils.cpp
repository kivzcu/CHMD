#include "DCUtils.h"

DCUtils::DCUtils(void)
{
}

DCUtils::~DCUtils(void)
{
}

int DCUtils::PointIsWithinBounds(double point[3], double bounds[6], double delta[3])
{
  if(!point || !bounds || !delta)
    {
    return 0;
    }
  for(int i=0;i<3;i++)
    {
    if(point[i]+delta[i] < bounds[2*i] || point[i]-delta[i] > bounds[2*i+1])
      {
      return 0;
      }
    }
  return 1;
}

// Computes Line-Line distance between lines defined by (p11, p12) and (p21, p22)
double DCUtils::LineLineDistance(double* p11, double* p12, double* p21, double* p22)
{
	double* a = new double[3];
	double* b = new double[3];
	double* c = new double[3];

	Substract(p12, p11, a);
	Substract(p22, p21, b);
	Substract(p21, p11, c);

	double* aCrossB = new double[3];
	vtkMath::Cross(a, b, aCrossB);

	return abs(vtkMath::Dot(c, aCrossB) / sqrt(aCrossB[0]*aCrossB[0] + aCrossB[1]*aCrossB[1] + aCrossB[2]*aCrossB[2]));
}
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

// Substracts point1 from point2 to result
void DCUtils::Substract(double* point1, double* point2, double* result)
{
	for (int i = 0; i < 3; i++)
	{
		result[i] = point2[i] - point1[i];
	}
}

void DCUtils::Add(double* point1, double* point2, double* result, double multiplication)
{
	for (int i = 0; i < 3; i++)
	{
		result[i] = point2[i] * multiplication + point1[i];
	}
}

double DCUtils::PointPointDistance(double* point1, double* point2)
{
	return sqrt(vtkMath::Distance2BetweenPoints(point1, point2));
}

// Computes normal of from 3 points (point1 being the normal origin)
void DCUtils::ComputeNormal(double* point1, double* point2, double* point3, double* normal)
{
	double u[3];
	double v[3];
	double result[3];
	Substract(point2, point1, u);
	Substract(point3, point1, v);
	vtkMath::Cross(u, v, result);
	vtkMath::Normalize(result);
	normal[0] = result[0];
	normal[1] = result[1];
	normal[2] = result[2];
}