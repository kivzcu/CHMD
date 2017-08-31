#include "Point.h"

int Point::count = 0;

//----------------------------------------------------------------------------
Point::Point(double *pCoord)
	//----------------------------------------------------------------------------
{
	DCoord[0] = pCoord[0]; //x
	DCoord[1] = pCoord[1]; //y
	DCoord[2] = pCoord[2]; //z
	double t;
	BBoundary = false;
	BMuscle = false;
	TendonConnectionID = -1;
	interpolated = false;
	distToNext = 0;
	count++;
}

//----------------------------------------------------------------------------
Point::~Point()
//----------------------------------------------------------------------------
{
	count--;
}

//----------------------------------------------------------------------------
// Computes coordinates of point on interpolated curve, using given x-coordinate
//----------------------------------------------------------------------------
void Point::ComputeInterpolatedCurveCoord(double t, double * result, bool substract)
{
	if (!this->interpolated)
	{
		result[0] = result[1] = result[2] = 0;
		return;
	}
	double h = t;
	if (substract)
	{
		h = (t - this->t);
	}
	double u = 0;

	for (int axis = 0; axis < 3; axis++)
	{
		u = a[axis] * h * h * h + b[axis] * h * h + c[axis] * h + d[axis];
		result[axis] = u;
	}
}

// Assigns normal to this Point
void Point::SetNormal(double* normal)
{
	this->DNormal[0] = normal[0];
	this->DNormal[1] = normal[1];
	this->DNormal[2] = normal[2];
}