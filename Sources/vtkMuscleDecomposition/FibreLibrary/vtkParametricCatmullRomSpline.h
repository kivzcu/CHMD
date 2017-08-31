/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkParametricCatmullRomSpline.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkParametricCatmullRomSpline - parametric function for 1D interpolating splines
// .SECTION Description
// vtkParametricCatmullRomSpline is a parametric function for 1D interpolating Catmull-Rom
// splines. Catmul-Rom defines its parameterization (knots) as ti+1 = ti + ||Pi+1 - Pi||^alpha 
// whereby alpha = 1/2 for centripetal CR, 1 for chordal, and 0 for uniform spline
// see: http://en.wikipedia.org/wiki/Centripetal_Catmull%E2%80%93Rom_spline
//
// vtkParametricCatmullRomSpline maps the single parameter u into a 3D point (x,y,z).
// This family of 1D splines is guaranteed to be parameterized in the interval [0,1]. 
// Attempting to evaluate outside this interval will cause the parameter u to be 
// clamped in the range [0,1].
//
//
// .SECTION Caveats
// If you wish to tessellate the spline, use the class vtkParametricFunctionSource.
//
// .SECTION See Also
// vtkCatmullRomSpline

#ifndef vtkParametricCatmullRomSpline_h
#define vtkParametricCatmullRomSpline_h

class vtkSpline;
class vtkPoints;

#include "vtkParametricFunction.h"

//compatibility for VTK <= 7.0.0
#ifndef VTK_DELETE_FUNCTION
#define VTK_DELETE_FUNCTION =delete
#endif

class vtkParametricCatmullRomSpline : public vtkParametricFunction
{
public:
	enum SplineTypes
	{
		Uniform = 0,
		Chordal = 1,
		Centripetal = 2,
		Custom = 3,
	};

public:
  vtkTypeMacro(vtkParametricCatmullRomSpline,vtkParametricFunction);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Construct the spline with the following parameters:
  // MinimumU = 0, MaximumU = 1, JoinU = 0 (unless the spline is
  // closed, then JoinU = 1), TwistU = 0, DerivativesSupplied = 0
  // (the other vtkParametricFunction parameters are ignored).
  static vtkParametricCatmullRomSpline *New();

  // Description
  // Return the parametric dimension of the class.
  virtual int GetDimension() {return 1;}

  // Description:
  // Evaluate the spline at parametric coordinate u[0] returning
  // the point coordinate Pt[3].
  virtual void Evaluate(double u[3], double Pt[3], double Du[9]);

  // Description:
  // Evaluate a scalar value at parametric coordinate u[0] and Pt[3].
  // The scalar value is just the parameter u[0].
  virtual double EvaluateScalar(double u[3], double Pt[3], double Du[9]);

  // Description:
  // Specify the list of points defining the spline. Do this by
  // specifying a vtkPoints array containing the points. Note that
  // the order of the points in vtkPoints is the order that the
  // splines will be fit.
  void SetPoints(vtkPoints*);
  vtkGetObjectMacro(Points,vtkPoints);

  // Description:
  // Another API to set the points. Set the number of points and then set the
  // individual point coordinates.
  void SetNumberOfPoints(vtkIdType numPts);
  void SetPoint(vtkIdType index, double x, double y, double z);

  // Description:
  // Control whether the spline is open or closed. A closed spline forms
  // a continuous loop: the first and last points are the same
  vtkSetMacro(Closed,int);
  vtkGetMacro(Closed,int);
  vtkBooleanMacro(Closed,int);

  // Description:
  // Controls how the spline is parameterised.
  // Default is Centripetal.
  vtkSetMacro(SplineType, SplineTypes);
  vtkGetMacro(SplineType, SplineTypes);  

	// Description:
	// Specifies alpha parameter for custom parameterized Catmull-Rom splines.
	vtkSetMacro(Alpha, double);
	vtkGetMacro(Alpha, double);

protected:
  vtkParametricCatmullRomSpline();
  ~vtkParametricCatmullRomSpline();

  // Points definition
  vtkPoints *Points;

  // The interpolating splines for each of the x-y-z coordinates
  vtkSpline *XSpline;
  vtkSpline *YSpline;
  vtkSpline *ZSpline;

  // Supplemental variables
  int    Closed;
  enum SplineTypes SplineType;
	double Alpha;	//<Alpha parameter to calculate the knot values of points

  // Initializing the spline
  unsigned long InitializeTime;
  int Initialize();

  // Internal variable for managing parametric coordinates
  double Length;
  double ClosedLength;

private:
  vtkParametricCatmullRomSpline(const vtkParametricCatmullRomSpline&) VTK_DELETE_FUNCTION;
  void operator=(const vtkParametricCatmullRomSpline&) VTK_DELETE_FUNCTION;
};

#endif
