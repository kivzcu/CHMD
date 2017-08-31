#include "vtkMuscleDecomposerSplineInterpolator.h"
#include <assert.h>
#include <float.h>

//linux finite wrapper
#ifdef __unix__
	#define _CHECK_VALUE(x) assert(finite(x));
#else
	#define _CHECK_VALUE(x) assert(_finite(x));
#endif

//----------------------------------------------------------------------------
vtkMuscleDecomposer::SplineInterpolator::SplineInterpolator(int splineType)
//----------------------------------------------------------------------------
{
	this->splineType = splineType;
}

//----------------------------------------------------------------------------
vtkMuscleDecomposer::SplineInterpolator::~SplineInterpolator()
//----------------------------------------------------------------------------
{
	derivations->clear();
	diagonal->clear();
	rightHand->clear();
	result->clear();
	subdiagonal.clear();
	superdiagonal.clear();
}

//----------------------------------------------------------------------------
// Interpolates a single fiber
void vtkMuscleDecomposer::SplineInterpolator::Interpolate(vtkMuscleDecomposer::Fiber* fiber, std::vector<vtkMuscleDecomposer::Point*> & points)
//----------------------------------------------------------------------------
{
	this->fiber = fiber;
	this->pointCloud = points;
	n = (int)fiber->PointIDs.size() - 2;

	
	switch (this->splineType)
	{
		case SPLINE_NATURAL:
		case SPLINE_PARABOLIC_RUNOUT:
		case SPLINE_CUBIC_RUNOUT:
		default:
			BuildParameter(false); // classic methods use chordal parametrisation, not centripetal
			BuildEquation();	// these solve the equations
			SolveEquation();
			ComputeWeights();	// compute and assign a,b,c,d weights to input points
			break;
		case SPLINE_CONSTRAINED:
			BuildParameter(true);
			ComputeFibreDerivations();
			ComputeConstrainedWeights();
			break;
		case SPLINE_CATMULL_ROM:
			BuildParameter(true);
			ComputeCatmullRom();
			break;
	}
}

//----------------------------------------------------------------------------
// Builds Spline interpolation parameter for all the points in the fibre, and
// associates them to the points in cloud using chordal parametrisation.
// The parameter is not normalised to t=<0,1>
//----------------------------------------------------------------------------
void vtkMuscleDecomposer::SplineInterpolator::BuildParameter(bool useCentripetal)
{
	double currDistance = 0;
	double segmentLength = 0;
	double prevDistance = 0;
	
	int currID = this->fiber->PointIDs[0];
	int prevID = 0;

	this->pointCloud[currID]->t = 0;

	for(int i = 1; i < (int)this->fiber->PointIDs.size(); i++)
	{
		prevID = currID;
		currID = this->fiber->PointIDs[i];

		segmentLength = sqrt(vtkMath::Distance2BetweenPoints(this->pointCloud[currID]->DCoord, this->pointCloud[prevID]->DCoord));
		if (useCentripetal)
		{
			segmentLength = sqrt(segmentLength);
		}

		if (segmentLength == 0) // this shouldn't happen, unless there is an error in the data
		{
			segmentLength = DBL_MIN;
		}
		currDistance = prevDistance + segmentLength;
		pointCloud[currID]->t = currDistance;
		prevDistance = currDistance;
	}
}

//----------------------------------------------------------------------------
// Builds Spline interpolation equations described in
// "Cubic Spline Interpolation, MAE 5093, Charles O’Neill, 28 May 2002"
// modified to be parametric as At = b (tridiagonal matrix equation set)
// for all three axes (X, Y, Z) in linear time.
//----------------------------------------------------------------------------
void vtkMuscleDecomposer::SplineInterpolator::BuildEquation()
{
	// prev, this, next refer to values in matrix/vector, not curve point - not anymore
	double t_prev, t_this, t_next = 0;  // independent parameter (t)
	double u_prev, u_this, u_next = 0;  // dependant value (x, y or z value)
	double h_prev = 0;
	double h_this = 0;

	// Matrix A is the same for all 3 axes (calculation depends only on t)
	superdiagonal.resize(n); // the last one is actually unused
	subdiagonal.resize(n); // the first one is actually unused

	for (int axis = 0; axis < 3; axis++) // we have 3 axes, 3 right sides, and 3 diagonals (modified by solver)
	{
		diagonal[axis].resize(n);
		rightHand[axis].resize(n);
	}

	for (int i = 0; i < n; i++)
	{
		t_prev = pointCloud[fiber->PointIDs[i]]->t;
		t_this = pointCloud[fiber->PointIDs[i+1]]->t;
		t_next = pointCloud[fiber->PointIDs[i+2]]->t;
		h_prev = t_this - t_prev;
		h_this = t_next - t_this;

		// matrices are the same for both axes, we can do them in one pass
		// refer to paper above for details
		// because we modify diagonal in solving process, it is copied for every axis
		for (int axis = 0; axis < 3; axis++)
		{
			diagonal[axis][i] = 2 * (h_prev + h_this) ;
			// different spline types have different border conditions
			if ((i == 0 || i == n - 1))
			{
				switch (this->splineType)
				{
					case SPLINE_CUBIC_RUNOUT:

						if (i == 0)
						{
							diagonal[axis][i] = 2 * (h_prev + h_this) + 2 * h_prev;
							superdiagonal[0] = 0;
						}
						else
						{
							diagonal[axis][i] = 2 * (h_prev + h_this) + 2 * h_this;
							subdiagonal[subdiagonal.size() - 1] = 0;
						}
						break;
					case SPLINE_PARABOLIC_RUNOUT:
						if (i == 0)
						{
							diagonal[axis][i] = 2 * (h_prev + h_this)+ h_prev;
						}
						else
						{
							diagonal[axis][i] = 2 * (h_prev + h_this) + h_this;
						}
						break;
					case SPLINE_NATURAL:
					default:
						break;
				}
			}
		}
		subdiagonal[i] = h_prev;
		superdiagonal[i] = h_this;

		// right hand vectors differ for each axis, we need three separate passes
		for (int axis = 0; axis < 3; axis++)
		{
			u_prev = pointCloud[fiber->PointIDs[i]]->DCoord[0 + axis];
			u_this = pointCloud[fiber->PointIDs[i+1]]->DCoord[0 + axis];
			u_next = pointCloud[fiber->PointIDs[i+2]]->DCoord[0 + axis];
			rightHand[axis][i] = ((u_next - u_this) / h_this - (u_this - u_prev) / h_prev) * 6;
		}
		// equation built for 3 axes and with low memory usage
	}
}

//----------------------------------------------------------------------------
// Solves equation for all axes using TDMA algorithm
void vtkMuscleDecomposer::SplineInterpolator::SolveEquation()
//----------------------------------------------------------------------------
{
	for (int axis = 0; axis < 3; axis++) // we have 3 axes
	{
		result[axis].resize(n);

		// super: accessed on indices 0..n-2 (n-1 is max index)
		// sub: accessed on indices 1..n-1 (n-1 is max index)
		double m;

		for (int i = 1; i < n; i++) // modify coefficients, forward pass
		{
			m = subdiagonal[i]/diagonal[axis][i-1];
			diagonal[axis][i] = diagonal[axis][i] - m * superdiagonal[i-1];
			rightHand[axis][i] = rightHand[axis][i] - m * rightHand[axis][i-1];
		}

		result[axis][n-1] = rightHand[axis][n-1]/diagonal[axis][n-1];

		for (int i = n - 2; i >= 0; i--) // indexed from zero, solve equations, backward pass
			result[axis][i]=(rightHand[axis][i] - superdiagonal[i] * result[axis][i+1]) / diagonal[axis][i];

		// result[axis] now holds S_i .. S_n
		diagonal[axis].clear();
		rightHand[axis].clear();
	}
	superdiagonal.clear();
	subdiagonal.clear();
}

//----------------------------------------------------------------------------
// Computes weights for all 3 axes
void vtkMuscleDecomposer::SplineInterpolator::ComputeWeights()
//----------------------------------------------------------------------------
{
	int ID;
	// next rows now refer to position in curve
	double S_this, S_next;	// computed by TDMA solver
	double t_this, t_next;  // computed in BuildParameter
	double u_this, u_next;  // coords of interpolated points
	double h;

	double S2, S3, SNm1, SNm2;	// for boundary conditions

	// First we add end point variables to the result
	for(int axis = 0; axis < 3; axis++)
	{
		switch (this->splineType)
		{
		case SPLINE_NATURAL: // natural splines have S=0 for end points
			result[axis].insert(result[axis].begin(), 0);
			result[axis].push_back(0);
			break;

		case SPLINE_PARABOLIC_RUNOUT: // parabolic runout splines have S equal to next(previous) S on end points
			S2 = result[axis][0];
			result[axis].insert(result[axis].begin(), S2); // first point copy
			SNm1 = result[axis][result[axis].size()-1];
			result[axis].push_back(SNm1); // last point copy
			break;
		case SPLINE_CUBIC_RUNOUT: // cubic runout splines have S1 = 2*S2 - S3 and SN = 2*SN-1 - SN-2
			S2 = result[axis][0];
			S3 = result[axis][1];
			result[axis].insert(result[axis].begin(), 2*S2 - S3); // first point calc
			SNm1 = result[axis][result[axis].size()-1];
			SNm2 = result[axis][result[axis].size()-2];
			result[axis].push_back(2*SNm1 - SNm2); // last point calc
			break;
		}
	}

	// Next, we compute weights for all points in fibre
	for (int i = 0; i < (int)fiber->PointIDs.size() - 1; i++) // last point doesnt have it's segment
	{
		ID = fiber->PointIDs[i];
		for(int axis = 0; axis < 3; axis++)
		{
			// again refer to paper above for details
			S_this = result[axis][i];
			S_next = result[axis][i+1];

			t_this = pointCloud[fiber->PointIDs[i]]->t;
			t_next = pointCloud[fiber->PointIDs[i+1]]->t;
			h = t_next - t_this;

			u_this = pointCloud[fiber->PointIDs[i]]->DCoord[axis];
			u_next = pointCloud[fiber->PointIDs[i+1]]->DCoord[axis];

			pointCloud[ID]->a[axis] = (S_next - S_this) / (6 * h);
			pointCloud[ID]->b[axis] = S_this / 2;
			pointCloud[ID]->c[axis] = ((u_next - u_this) / h) - ((2 * h * S_this + h * S_next) / 6);
			pointCloud[ID]->d[axis] = u_this;
			pointCloud[ID]->DistToNext = h;

			_CHECK_VALUE(h);
		}
		pointCloud[ID]->Interpolated = true;
	}
}

//----------------------------------------------------------------------------
// Calculates derivation aproximation of fibre data derivatives in all axes
// with respect to parameter t
void vtkMuscleDecomposer::SplineInterpolator::ComputeFibreDerivations()
//----------------------------------------------------------------------------
{
	int size = (int)this->fiber->PointIDs.size();
	double t_prev, t_this, t_next = 0;  // independent parameter (t)
	double u_prev, u_this, u_next = 0;  // dependant value (x, y or z value)
	double h_prev, h_this = 0;
	double h_u_prev, h_u_this = 0;
	double derivation;

	for (int axis = 0; axis < 3; axis++)
	{
		derivations[axis].resize(size);
		for (int i = 1; i < size - 1; i++)  // border derivations will be calculated later
		{
			t_prev = pointCloud[fiber->PointIDs[i-1]]->t;
			t_this = pointCloud[fiber->PointIDs[i]]->t;
			t_next = pointCloud[fiber->PointIDs[i+1]]->t;
			h_prev = t_this - t_prev;
			h_this = t_next - t_this;

			u_prev = pointCloud[fiber->PointIDs[i-1]]->DCoord[axis];
			u_this = pointCloud[fiber->PointIDs[i]]->DCoord[axis];
			u_next = pointCloud[fiber->PointIDs[i+1]]->DCoord[axis];
			h_u_prev = u_this - u_prev;
			h_u_this = u_next - u_this;

			double temp = (h_this/h_u_this) * (h_prev/h_u_prev);
			// if sope sign changes at the point, the derivation should be = 0
			if (temp > 0)
			{
				derivation = 2 / ((h_this/h_u_this) + (h_prev/h_u_prev));
			} else
			{
				derivation = 0;
			}

			derivations[axis][i] = derivation;
		}
		// border derivation calc on fibre start
		t_this = pointCloud[fiber->PointIDs[0]]->t;
		t_next = pointCloud[fiber->PointIDs[1]]->t;
		h_this = t_next - t_this;
		u_this = pointCloud[fiber->PointIDs[0]]->DCoord[axis];
		u_next = pointCloud[fiber->PointIDs[1]]->DCoord[axis];
		h_u_this = u_next - u_this;

		derivations[axis][0] = ((3 * h_u_this) / (2 * h_this)) - (derivations[axis][1] / 2);
		// border derivation calc on fibre end
		t_this = pointCloud[fiber->PointIDs[size-1]]->t;
		t_prev = pointCloud[fiber->PointIDs[size-2]]->t;
		h_prev = t_this - t_prev;
		u_this = pointCloud[fiber->PointIDs[size-1]]->DCoord[axis];
		u_prev = pointCloud[fiber->PointIDs[size-2]]->DCoord[axis];
		h_u_prev = u_this - u_prev;

		derivations[axis][size-1] = ((3 * h_u_prev) / (2 * h_prev)) - (derivations[axis][size - 2] / 2);
	}
	std::vector<double> x = derivations[0];;
	//x = //.begin(), derivation[0].end());
	std::vector<double> y = derivations[1];
}

//----------------------------------------------------------------------------
// Calculates parameters a, b, c, d for constrained splines
// as described in "Constrained Cubic Spline Interpolation
// for Chemical Engineering Applications" by CJC Kruger
// As the first derivatives are known, we dont need to solve any equations.
// (second derivatives are computed from first derivatives)
void vtkMuscleDecomposer::SplineInterpolator::ComputeConstrainedWeights()
//----------------------------------------------------------------------------
{
	int ID;
	double t_this, t_next = 0;  // independent parameter (t)
	double u_this, u_next = 0;  // dependant value (x, y or z value)
	double h_this = 0;
	double h_u_this = 0;
	double this_2ndDerivation;
	double next_2ndDerivation;
	double a, b, c, d;

	for (int i = 0; i < (int)fiber->PointIDs.size() - 1; i++) // last point doesnt have it's segment
	{
		ID = fiber->PointIDs[i];

		for(int axis = 0; axis < 3; axis++)
		{
			// refer to paper above for details
			t_next = pointCloud[fiber->PointIDs[i+1]]->t;
			t_this = pointCloud[fiber->PointIDs[i]]->t;

			h_this = t_next - t_this;

			u_next = pointCloud[fiber->PointIDs[i+1]]->DCoord[axis];
			u_this = pointCloud[fiber->PointIDs[i]]->DCoord[axis];

			h_u_this = u_next - u_this;

			this_2ndDerivation = (-2 * (derivations[axis][i+1] + 2*derivations[axis][i]) / (h_this)) + (6*(h_u_this) / (h_this * h_this));
			next_2ndDerivation = (2 * (2 * derivations[axis][i+1] + derivations[axis][i]) / (h_this)) - (6*(h_u_this) / (h_this * h_this));

			d = (next_2ndDerivation - this_2ndDerivation) / (6 * h_this);
			c = (t_next * this_2ndDerivation - t_this * next_2ndDerivation ) / (2 * h_this);
			b = (h_u_this - c*(t_next*t_next - t_this*t_this) - d*(t_next*t_next*t_next - t_this*t_this*t_this))/(h_this);
			a = u_this - (b*t_this) - (c * t_this * t_this) - (d * t_this * t_this * t_this);

			pointCloud[ID]->a[axis] = d;
			pointCloud[ID]->b[axis] = c;
			pointCloud[ID]->c[axis] = b;
			pointCloud[ID]->d[axis] = a;
			pointCloud[ID]->DistToNext = h_this;

			_CHECK_VALUE(h_this);
		}
		pointCloud[ID]->Interpolated = true;
	}
}

//----------------------------------------------------------------------------
// Prepares fibre for non-uniform Catmull-Rom spline interpolation.
// Unfortunately, its computation is very different from traditional spline
// (employs pyramid recursive computation), hence those insane computations later.
void vtkMuscleDecomposer::SplineInterpolator::ComputeCatmullRom()
//----------------------------------------------------------------------------
{
	int ID; // id of processed curve
	double t0, t1, t2, t3;  // independent parameter (t)
	double t10, t20, t30, t21, t31, t32; // helpers

	double* P03A; double* P13A; double* P23A; double* P33A; // 4 knots, 3Axis
	double P0, P1, P2, P3; // 4 knots, one axis (2nd is t)
	double a, b, c, d; // see spline interpolation
	
	double tmp[3]; // static
	double* pt0; // pointer only, helper
	double* pt1; // pointer only, helper

	for (int i = 0; i < (int)fiber->PointIDs.size() - 1; i++) // last point doesnt have it's segment, 
		//(i-1)----(!i!)----(i+1)----(i+2), computing segment i->i+1, storing data to i
	{
		ID = fiber->PointIDs[i];
		if (i > 0)
		{
			P03A = pointCloud[fiber->PointIDs[i-1]]->DCoord;
			t0 = pointCloud[fiber->PointIDs[i-1]]->t;
		} else // i == 0 -> we need a pt_-1 point -> mirror pt_1 around pt_0
		{	// [TechNote TN-06-001: Catmull-Rom Splines, Jim Armstrong, Singularity, January 2006, online]
			t0 = -pointCloud[fiber->PointIDs[i+1]]->t;
			
			pt0 = pointCloud[fiber->PointIDs[i]]->DCoord;
			pt1 = pointCloud[fiber->PointIDs[i+1]]->DCoord;
			tmp[0] = pt1[0] - pt0[0]; //diff
			tmp[1] = pt1[1] - pt0[1];
			tmp[2] = pt1[2] - pt0[2];
			
			tmp[0] = pt0[0] - tmp[0]; //mirror
			tmp[1] = pt0[1] - tmp[1];
			tmp[2] = pt0[2] - tmp[2];

			P03A = tmp;
		}

		P13A = pointCloud[fiber->PointIDs[i]]->DCoord;
		t1 = pointCloud[fiber->PointIDs[i]]->t;

		P23A = pointCloud[fiber->PointIDs[i+1]]->DCoord;
		t2 = pointCloud[fiber->PointIDs[i+1]]->t;
		
		if (i <  (int)fiber->PointIDs.size() - 2) // next to last
		{
			P33A = pointCloud[fiber->PointIDs[i+2]]->DCoord;
			t3 = pointCloud[fiber->PointIDs[i+2]]->t;
		} else // we are on the last point, -> we need a pt_2 point (not in data) -> mirror pt_0 around pt_1
		{
			t3 = t2 + (t2 - t1); // this seems to work (extend the parameter beyond, it not used for actual interpolation)
			
			pt0 = pointCloud[fiber->PointIDs[i]]->DCoord;
			pt1 = pointCloud[fiber->PointIDs[i+1]]->DCoord;
			tmp[0] = pt1[0] - pt0[0]; // diff
			tmp[1] = pt1[1] - pt0[1];
			tmp[2] = pt1[2] - pt0[2];
			
			tmp[0] = pt1[0] + tmp[0]; // mirror (flipped vector here at the end)
			tmp[1] = pt1[1] + tmp[1];
			tmp[2] = pt1[2] + tmp[2];

			P33A = tmp;
		}
		// some helpers (simplifies the later insanity)
		t10 = t1-t0;
		t20 = t2-t0;
		t30 = t3-t0;
		t21 = t2-t1;
		t31 = t3-t1;
		t32 = t3-t2;
		
		for(int axis = 0; axis < 3; axis++)
		{
			P0 = P03A[axis];
			P1 = P13A[axis];
			P2 = P23A[axis];
			P3 = P33A[axis];

			// Insane? Yeah.... :\ Courtesy of Wolfram Mathematica and a lot of manual tweaking. But it works.
			// Should be the simplest de-recursion (so the alg. can fit the finished interpolator) of recursive pyramid algorithm
			// Why? Remember, we do not compute the interpolation here, just the coeficients according to powers of t.
			// [A Recursive Evaluation Algorithm for a Class of Catmull-Rom Splines, J.Barry et. al, Computer Graphics '88]

			a = (P1/(t20*t21*t21) - P2/(t20*t21*t21) - P0/(t10*t20*t21) + P1/(t10*t20*t21) + P1/(t21*t21*t31) - P2/(t21*t21*t31) - P2/(t21*t31*t32) + P3/(t21*t31*t32));

			b = (-((P1*t0)/(t20*t21*t21)) + (P2*t0)/(t20*t21*t21) + (P2*t1)/(t20*t21*t21) - (2*P1*t2)/(t20*t21*t21) + (P2*t2)/(t20*t21*t21) - (P1*t0)/(t10*t20*t21) + (P0*t1)/(t10*t20*t21) + 
				(2*P0*t2)/(t10*t20*t21) - (2*P1*t2)/(t10*t20*t21) - (P1*t1)/(t21*t21*t31) + (2*P2*t1)/(t21*t21*t31) - (P1*t2)/(t21*t21*t31) - (P1*t3)/(t21*t21*t31) + (P2*t3)/(t21*t21*t31) + 
				(2*P2*t1)/(t21*t31*t32) - (2*P3*t1)/(t21*t31*t32) - (P3*t2)/(t21*t31*t32) + (P2*t3)/(t21*t31*t32));

			c = (-((P2*t0*t1)/(t20*t21*t21)) + (2*P1*t0*t2)/(t20*t21*t21) - (P2*t0*t2)/(t20*t21*t21) - (P2*t1*t2)/(t20*t21*t21) + (P1*t2*t2)/(t20*t21*t21) + (2*P1*t0*t2)/(t10*t20*t21) - 
				(2*P0*t1*t2)/(t10*t20*t21) - (P0*t2*t2)/(t10*t20*t21) + (P1*t2*t2)/(t10*t20*t21) - (P2*t1*t1)/(t21*t21*t31) + (P1*t1*t2)/(t21*t21*t31) + (P1*t1*t3)/(t21*t21*t31) - 
				(2*P2*t1*t3)/(t21*t21*t31) + (P1*t2*t3)/(t21*t21*t31) - (P2*t1*t1)/(t21*t31*t32) + (P3*t1*t1)/(t21*t31*t32) + (2*P3*t1*t2)/(t21*t31*t32) - (2*P2*t1*t3)/(t21*t31*t32));

			d = (P2*t0*t1*t2)/(t20*t21*t21) - (P1*t0*t2*t2)/(t20*t21*t21) - (P1*t0*t2*t2)/(t10*t20*t21) + (P0*t1*t2*t2)/(t10*t20*t21) + (P2*t1*t1*t3)/(t21*t21*t31) - (P1*t1*t2*t3)/(t21*t21*t31) -
				(P3*t1*t1*t2)/(t21*t31*t32) + (P2*t1*t1*t3)/(t21*t31*t32);


			pointCloud[ID]->a[axis] = a;
			pointCloud[ID]->b[axis] = b;
			pointCloud[ID]->c[axis] = c;
			pointCloud[ID]->d[axis] = d;
			pointCloud[ID]->DistToNext = t21;

//			_CHECK_VALUE(t21);
		}
		pointCloud[ID]->Interpolated = true;
	}
}

//----------------------------------------------------------------------------
// Just a getter
int vtkMuscleDecomposer::SplineInterpolator::GetInterpolationType()
//----------------------------------------------------------------------------
{
	return this->splineType;
}
