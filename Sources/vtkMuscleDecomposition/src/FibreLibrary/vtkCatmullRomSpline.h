/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkCatmullRomSpline.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkCatmullRomSpline - computes an interpolating spline using a
// non-uniform cubic Catmull-Rom spline basis.

// .SECTION Description
// vtkCatmullRomSpline is a concrete implementation of vtkSpline using a
// non-uniform cubic Catmull-Rom spline basis. Catmull-Rom spline is
// defined by N points P0 .. PN-1 and N knots values at these points
// t0 .. tn-1. Points and knots (Pi-1, ti-1), (Pi, ti), (Pi+1,ti+1), 
// P(i+2, ti+2) define the i-the cubic polynomial segment of the spline
// with t from ti to ti+1. If the spline is open, the the spline starts
// LeftConstraint and RightConstraint can be specified to determine
// how the spline should interpolate the first and the last point.
//
// The spline formula is defined by:
// P. J. Barry and R. N. Goldman. A recursive evaluation algorithm for 
// a class of Catmull–Rom splines. SIGGRAPH Computer Graphics, 22(4):199-204, 1988.
//
// Warning: it must be guaranteed that t0 < t1 < ... ti < ti+1 < ... < tn-1

// .SECTION See Also
// vtkSpline vtkCardinalSpline vtkKochanekSpline


#ifndef vtkCatmullRomSpline_h
#define vtkCatmullRomSpline_h

#include "vtkSpline.h"

//compatibility for VTK <= 7.0.0
#ifndef VTK_DELETE_FUNCTION
#define VTK_DELETE_FUNCTION =delete
#endif

class vtkCatmullRomSpline : public vtkSpline
{
public:
  static vtkCatmullRomSpline *New();

  vtkTypeMacro(vtkCatmullRomSpline,vtkSpline);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description
  // Compute Cardinal Splines for each dependent variable
  void Compute ();	

  // Description:
  // Evaluate a 1D Catmull-Rom spline.
  /*virtual*/ double Evaluate (double t);
	
	// Description:
	// Evaluate the first derivative of 1D Catmull-Rom spline.
	virtual double EvaluateDerivative(double t);

	// Description:
	// Evaluate the second derivative of 1D Catmull-Rom spline.
	virtual double EvaluateDerivative2(double t);

	// Description:
	// Set the type of constraint of the left(right) end points. These
	// constraints are available:
	//	
	// 0: an extra point is automatically constructed at left(right) 
	// part with its position and knot value derived from the 
	// data of the first(last) two points as P_{-1} = P0, 
	// t_{-1} = t0 + (t0 - t1)*tension, to determine the shape of spline from 
	// P0 to P1 (Pn-2 to Pn-1). Tension is specified in LeftValue(RightValue).
	// If the tension is zero, the spline will interpolate
	// the points P1 .. Pn-2 only whereas the points P0 and Pn-1 
	// will influence its shape. Default values for the tension is 1/256.
	//
	// 1: knot values of two extra points (or in closed spline for the repeated
	// points) are given in LeftValue and RightValue. Zero values mean that
	// the spline will not interpolate to the end points.
	// N.B. contraint 1 is the default.
	vtkSetClampMacro(LeftConstraint, int, 0, 1);
	vtkSetClampMacro(RightConstraint, int, 0, 1);

  // Description:
  // Deep copy of cardinal spline data.
  /*virtual*/ void DeepCopy(vtkSpline *s);

protected:
  vtkCatmullRomSpline();
  ~vtkCatmullRomSpline() {}

	//Gets the coefficents for the given t.
	//If t is outside the supported range, it is clamped.
	//If the spline is invalid, the method returns false, otherwise true.
	bool GetCoefficients(double& t, double*& coefficients);

	//Compute the coefficients for a 1D spline.
  void Fit1D ();
	
	//calculates the knot value for the point Pc given knot values ta and tb
	inline double ComputeRightKnot(double ta, double tb) {
		return this->RightConstraint != 0 ? this->RightValue : tb + (tb - ta)*this->RightValue;
	}
	
	//calculates the knot value for the point Pa given knot values tb and tc
	inline double ComputeLeftKnot(double tb, double tc) {
		return this->LeftConstraint != 0 ? this->LeftValue : tb - (tc - tb)*this->LeftValue;		
	}	

private:
  vtkCatmullRomSpline(const vtkCatmullRomSpline&) VTK_DELETE_FUNCTION;
  void operator=(const vtkCatmullRomSpline&) VTK_DELETE_FUNCTION;
};

#endif

