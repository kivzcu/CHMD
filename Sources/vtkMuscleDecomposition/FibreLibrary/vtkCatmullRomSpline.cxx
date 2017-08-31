/*=========================================================================

	Program:   Visualization Toolkit
	Module:    vtkCatmullRomSpline.cxx

	Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
	All rights reserved.
	See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

		 This software is distributed WITHOUT ANY WARRANTY; without even
		 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
		 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkCatmullRomSpline.h"

#include "vtkObjectFactory.h"
#include "vtkPiecewiseFunction.h"
#include <cassert>

vtkStandardNewMacro(vtkCatmullRomSpline);

//----------------------------------------------------------------------------
// Construct a Cardinal Spline.
vtkCatmullRomSpline::vtkCatmullRomSpline()
{
	this->LeftValue = this->RightValue = 0.00390625; //= 1 / 256
}

//Gets the coefficents for the given t.
//If t is outside the supported range, it is clamped.
//If the spline is invalid, the method returns false, otherwise true.
bool vtkCatmullRomSpline::GetCoefficients(double& t, double*& coefficients)
{
	int size = this->PiecewiseFunction->GetSize();
	if (size == 0)	//empty spline
		return false;

	// check to see if we need to recompute the spline
	if (this->ComputeTime < this->GetMTime()) {	
		this->Compute();
	}	

	if (this->Closed != 0 && size > 1)
	{
		size = size + 1;
	}

	// clamp the function at both ends
	if (t < this->Intervals[0])
	{
		t = this->Intervals[0];
	}
	if (t > this->Intervals[size - 1])
	{
		t = this->Intervals[size - 1];
	}

	// find pointer to cubic spline coefficient using bisection method
	int index = this->FindIndex(size, t);
	coefficients = this->Coefficients + 4*index;
	return true;
}

//----------------------------------------------------------------------------
// Evaluate a 1D Spline
double vtkCatmullRomSpline::Evaluate(double t)
{
	double *coefficients;			
	if (!GetCoefficients(t, coefficients))
		return 0.0;

	// evaluate intervals value y
	return (t * (t * (t * *(coefficients)
		+ *(coefficients + 1))
		+ *(coefficients + 2))
		+ *(coefficients + 3));
}

// Description:
// Evaluate the first derivative of 1D Catmull-Rom spline.
/*virtual*/ double vtkCatmullRomSpline::EvaluateDerivative(double t)
{
	double *coefficients;
	if (!GetCoefficients(t, coefficients))
		return 0.0;

	// evaluate intervals value y
	return (t * (t * *(coefficients)*3
		+ *(coefficients + 1)*2)
		+ *(coefficients + 2)
		);
}

// Description:
// Evaluate the second derivative of 1D Catmull-Rom spline.
/*virtual*/ double vtkCatmullRomSpline::EvaluateDerivative2(double t)
{
	double *coefficients;
	if (!GetCoefficients(t, coefficients))
		return 0.0;

	// evaluate intervals value y
	return (t * *(coefficients) * 6
		+ *(coefficients + 1) * 2);
}

//----------------------------------------------------------------------------
// Compute Catmull-Rom Splines for each dependent variable
void vtkCatmullRomSpline::Compute()
{
	//free the previous cofficients
	delete[] this->Intervals;
	delete[] this->Coefficients;
	this->Intervals = NULL;
	this->Coefficients = NULL;

	// get the size of the independent variables
	int size = this->PiecewiseFunction->GetSize();
	if (size == 0)
	{
		vtkErrorMacro("Cannot compute a spline with no points. ");
	}
	else
	{
		double* ts = this->PiecewiseFunction->GetDataPointer();
		if (size == 1)
		{
			//single point specified
			this->Intervals = new double[1]{ *ts };
			this->Coefficients = new double[4]{ 0, 0, 0, *(ts + 1) };
		}
		else
		{
			//we have here at least two points
			int szAlloc = (this->Closed == 0) ? size : size + 1;

			this->Intervals = new double[szAlloc];
			this->Coefficients = new double[4 * szAlloc];

			// copy the independent variables (intervals). 
			for (int i = 0; i < size; i++) {
				this->Intervals[i] = *(ts + 2 * i);
			}

			if (this->Closed != 0)
			{
				//if the spline is closed, add an extra interval
				if (this->ParametricRange[0] != this->ParametricRange[1])
				{
					//ParametricRange is explicitely set, i.e., we assume that the caller MUST 
					//have specified a valid knot value for the end-point
					this->Intervals[size] = this->ParametricRange[1];
				}
				else
				{
					//ParametricRange is not set, thus compute the right knot value as usual
					this->Intervals[size] = ComputeRightKnot(this->Intervals[size - 1], this->Intervals[size - 2]);
				}
			}

			// get the dependent variable values
			this->Fit1D();
		}
	}

	// update compute time
	this->ComputeTime = this->GetMTime();
}

//----------------------------------------------------------------------------
// Compute the coefficients for a 1D spline. 
void vtkCatmullRomSpline::Fit1D()
{
	int nPts = this->PiecewiseFunction->GetSize();
	double* ts = this->PiecewiseFunction->GetDataPointer();	//knots
	double* xs = ts + 1;	//points
	int span = 2;	//span between two values

	//Let the fibre be formed by points Q0, Q1, ... Qn-1 having Catmull-Rom parameterization Q0_t, Q1_t, ...Qn-1_t
	//Cubic Catmul-Rom curve C is defined by 4 points P0..P3 and their knots (parameterization) t0..t3
	//where knots ti+1 are typically defined as ti+1 = ti + ||Pi+1 - Pi||^alpha 
	//whereby alpha = 1/2 for centripetal CR, 1 for chordal
	//see: http://en.wikipedia.org/wiki/Centripetal_Catmull%E2%80%93Rom_spline
	//Note: the curve is valid only from P1 to P2 (points P0 and P3 are used to define it), i.e.,
	//the valid range of the parameter t is <t1...t2>	

	double P[4];	//buffer for the points P0..P3
	double t[4];	//buffer for the knots values t0..t3

	typedef double Vect4[4];
	Vect4* M = new Vect4[4];

	//P1 is the first point, P2 is the second
	P[1] = xs[0]; t[1] = ts[0]; P[2] = xs[span]; t[2] = ts[span];

	//load or construct P0	
	P[0] = (this->Closed == 0) ? P[1] : xs[(nPts - 1) * span];
	t[0] = ComputeLeftKnot(t[1], t[2]);

	int n = (this->Closed == 0) ? nPts - 1 : nPts;
	for (int j = 0; j < n; j++)
	{
		//load (or construct) the fourth point
		int iv0 = (j + 0) % 4;
		int iv1 = (j + 1) % 4;
		int iv2 = (j + 2) % 4;
		int iv3 = (j + 3) % 4;

		//for j == 0, P1 is Q0
		//for j == 1, P1 is Q1
		// ...
		//for j == nPts - 2, P1 is Qn-2 and P2 is Qn-1
		if (j + 2 < nPts)
		{
			P[iv3] = xs[(j + 2)*span];
			t[iv3] = ts[(j + 2)*span];
		}
		else
		{
			//for an open spline, we are currently at segment Qn-2 .. Qn-1 => construct the fourth point
			//for a closed spline, we are either at segment Qn-2 .. Qn-1 => the fourth point is Q0 and its knot
			//		value is either specified in ParametericRange[1] (SHOULD BE) is calculated from tn-2 and tn-1
			//or at segment Qn-1 .. Q0 => the fourth point is Q1 and its knot value is either specified in RightValue
			//or calculated from previous two knot values
			if (this->Closed == 0)
			{
				P[iv3] = P[iv2];
				t[iv3] = ComputeRightKnot(t[iv1], t[iv2]);
			}
			else if (j + 2 == nPts)
			{
				//segment Qn-2 .. Qn-1
				P[iv3] = xs[0];

				if (this->ParametricRange[0] != this->ParametricRange[1]) {
					t[iv3] = this->ParametricRange[1];
				}
				else {
					t[iv3] = ComputeRightKnot(t[iv1], t[iv2]);
				}
			}
			else
			{
				//segment Qn-1 .. Q0
				P[iv3] = xs[span];
				t[iv3] = ComputeRightKnot(t[iv1], t[iv2]);
			}
		}

		//OK, we have here the full definition of the Catmull-Rom cubic spline in P0..P3 with knots in t0..t3
		//and the valid curve segment is from P1 to P2

		//Catmull-Rom cubic curve C(t) formula is:
		//C(t) = (t2-t)/(t2-t1)*B1 + (t-t1)/(t2-t1)*B2
		//B1 = (t2-t)/(t2-t0)*A1 + (t-t0)/(t2-t0)*A2
		//B2 = (t3-t)/(t3-t1)*A2 + (t-t1)/(t3-t1)*A3
		//A1 = (t1-t)/(t1-t0)*P0 + (t-t0)/(t1-t0)*P1
		//A2 = (t2-t)/(t2-t1)*P1 + (t-t1)/(t2-t1)*P2
		//A3 = (t3-t)/(t3-t2)*P2 + (t-t2)/(t3-t2)*P3
		//
		//Using Matlab C(t) can be expressed as V*[t^3 t^2 t 1] where V is a nasty vector
		//[
		//((P0 - P1)/((t0 - t1)*(t0 - t2)) - (P1 - P2)/((t0 - t2)*(t1 - t2)))/(t1 - t2) - ((P1 - P2)/((t1 - t2)*(t1 - t3)) - (P2 - P3)/((t1 - t3)*(t2 - t3)))/(t1 - t2)
		//(P1*t2^2 - 2*P2*t1^2 - P1*t3^2 + 2*P3*t1^2 + P2*t3^2 - P3*t2^2 + P1*t1*t2 - P1*t1*t3 + P2*t1*t3 - P3*t1*t2)/((t1 - t2)^2*(t1 - t3)*(t2 - t3)) - (P0*t1^2 - P1*t0^2 - 2*P0*t2^2 + P2*t0^2 + 2*P1*t2^2 - P2*t1^2 + P0*t1*t2 - P1*t0*t2 + P2*t0*t2 - P2*t1*t2)/((t0 - t1)*(t0 - t2)*(t1 - t2)^2)
		//(P2*t1^3 - P3*t1^3 - P1*t1*t2^2 + P1*t1*t3^2 + P1*t2*t3^2 - P1*t2^2*t3 - 2*P2*t1*t3^2 + P2*t1^2*t3 + 2*P3*t1*t2^2 - P3*t1^2*t2)/((t1 - t2)^2*(t1 - t3)*(t2 - t3)) - (P0*t2^3 - P1*t2^3 + P0*t1*t2^2 - 2*P0*t1^2*t2 - P1*t0*t2^2 + 2*P1*t0^2*t2 + P2*t0*t1^2 - P2*t0^2*t1 - P2*t0^2*t2 + P2*t1^2*t2)/((t0 - t1)*(t0 - t2)*(t1 - t2)^2)
		//- (t2*((t2*(P0*t1 - P1*t0))/((t0 - t1)*(t0 - t2)) - (t0*(P1*t2 - P2*t1))/((t0 - t2)*(t1 - t2))))/(t1 - t2) + (t1*((t3*(P1*t2 - P2*t1))/((t1 - t2)*(t1 - t3)) - (t1*(P2*t3 - P3*t2))/((t1 - t3)*(t2 - t3))))/(t1 - t2)
		//]
		//
		//V can be defined as [P0 P1 P2 P3]*M where M is a 4x4 matrix of coefficients. These coefficients are
		//M11 = 1/((t0 - t1)*(t0 - t2)*(t1 - t2))						M12 = -(t1^2 + t1*t2 - 2*t2^2)/((t0 - t1)*(t0 - t2)*(t1 - t2)^2)													M13 = -(t2*(- 2*t1^2 + t1*t2 + t2^2))/((t0 - t1)*(t0 - t2)*(t1 - t2)^2)																																			M14 = -(t1*t2^2)/((t0 - t1)*(t0 - t2)*(t1 - t2))
		//M21 = -(t0 - t3)/((t0 - t1)*(t1 - t2)^2*(t1 - t3))		M22 = (2*t0*t1 + t0*t2 + t1*t2 - t1*t3 - 2*t2*t3 - t1^2)/((t0 - t1)*(t1 - t2)^2*(t1 - t3))			M23 = (- 2*t0^2*t2 + t0*t2^2 + t2^3)/((t0 - t1)*(t0 - t2)*(t1 - t2)^2) - (t2^2*t3 + t1*t2^2 - t2*t3^2 - t1*t3^2)/((t1 - t2)^2*(t1 - t3)*(t2 - t3))				M24 = -(t2*(t1^2*t3 - t0*t1*t2 - t0*t1*t3 + t0*t2*t3))/((t0 - t1)*(t1 - t2)^2*(t1 - t3))
		//M31 = (t0 - t3)/((t0 - t2)*(t1 - t2)^2*(t2 - t3))		M32 = -(2*t0*t1 + t0*t2 - t1*t2 - t1*t3 - 2*t2*t3 + t2^2)/((t0 - t2)*(t1 - t2)^2*(t2 - t3))		M33 = (t1^3 + t1^2*t3 - 2*t1*t3^2)/((t1 - t2)^2*(t1 - t3)*(t2 - t3)) - (- t0^2*t1 - t2*t0^2 + t0*t1^2 + t2*t1^2)/((t0 - t1)*(t0 - t2)*(t1 - t2)^2)				M34 = -(t1*(t0*t2^2 + t0*t1*t3 - t0*t2*t3 - t1*t2*t3))/((t0 - t2)*(t1 - t2)^2*(t2 - t3))
		//M41 = -1/((t1 - t2)*(t1 - t3)*(t2 - t3))						M42 = -(- 2*t1^2 + t1*t2 + t2^2)/((t1 - t2)^2*(t1 - t3)*(t2 - t3))												M43 = -(t1*(t1^2 + t1*t2 - 2*t2^2))/((t1 - t2)^2*(t1 - t3)*(t2 - t3))																																				M44 = (t1^2*t2)/((t1 - t2)*(t1 - t3)*(t2 - t3))								
		//
		// It may be seen that these coefficients have the form nominator / denominator and 
		// the terms in the nominators/denominators are shared between them.
		
		double dt01 = t[iv0] - t[iv1], dt02 = t[iv0] - t[iv2], dt03 = t[iv0] - t[iv3],
			dt12 = t[iv1] - t[iv2], dt13 = t[iv1] - t[iv3],
			dt23 = t[iv2] - t[iv3];

		double st11 = t[iv1] * t[iv1], st22 = t[iv2] * t[iv2];

		double a5 = -(t[iv1] + 2 * t[iv2]);
		double a6 = 2 * t[iv1] + t[iv2];

		double t01 = t[iv0] * t[iv1], t12 = t[iv1] * t[iv2];
		double denom = dt01*dt02*dt12;
		M[0][0] = 1 / denom; M[0][1] = a5 / denom; M[0][2] = t[iv2] * a6 / denom; M[0][3] = -t[iv1] * st22 / denom;

		denom = dt01*dt12*dt12*dt13;
		M[1][0] = -dt03 / denom;
		M[1][1] = (t[iv0] * a6 + t[iv1] * dt23 - 2 * t[iv2] * t[iv3] - st11) / denom;
		M[1][2] = (t12*(dt12 - dt03 - 2 * t[iv0]) - t[iv3] * (t[iv0] * dt12 - st11 - st22)) / denom;
		M[1][3] = (t12 * t[iv3] * dt01 + t[iv0] * st22 * dt13) / denom;
		
		denom = dt02*dt12*dt12*dt23;
		M[2][0] = dt03 / denom;
		M[2][1] = (-(t[iv2] * (dt01 - t[iv3] + dt23) + t[iv1] * (dt03 + t[iv0]))) / denom;
		M[2][2] = (t01 * (t[iv1] + t[iv2] + t[iv3]) + t[iv0] * t[iv2] * dt23 - t12*(dt12 + 3 * t[iv3])) / denom;
		M[2][3] = -(t01 * t[iv2] * dt23 + st11 * t[iv3] * dt02) / denom;
		
		denom = dt12*dt13*dt23;
		M[3][0] = -1 / denom; M[3][1] = a6 / denom; M[3][2] = t[iv1] * a5 / denom; M[3][3] = st11*t[iv2] / denom;
						
		double* vP = this->Coefficients + (4 * j);
		for (int k = 0; k < 4; k++) {
			vP[k] = P[iv0] * M[0][k] + P[iv1] * M[1][k] + P[iv2] * M[2][k] + P[iv3] * M[3][k];
		}		
	} //end for

	delete[] M;
}


//----------------------------------------------------------------------------
void vtkCatmullRomSpline::DeepCopy(vtkSpline *s)
{
	vtkCatmullRomSpline *spline = vtkCatmullRomSpline::SafeDownCast(s);

	if (spline != NULL)
	{
		//nothing to do
	}

	// Now do superclass
	this->vtkSpline::DeepCopy(s);
}

//----------------------------------------------------------------------------
void vtkCatmullRomSpline::PrintSelf(ostream& os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os, indent);
}

