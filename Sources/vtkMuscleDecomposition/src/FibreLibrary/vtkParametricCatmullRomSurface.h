/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkParametricCatmullRomSurface.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkParametricCatmullRomSurface - parametric function for 2D interpolating surface
// .SECTION Description
// vtkParametricCatmullRomSurface is a parametric function that maps (u,v) 
// parameters into a 3D point (x, y, z). Both parameters are in the interval [0,1]. 
// Attempting to evaluate outside this interval will cause the parameter to be 
// clamped in the range [0,1].
//
// The Catmull-Rom surface is defined by a set of points structured in 
// an irregular grid:
//
// P0, P1 ... Pa-1
// Pa, Pa+1 .. Pb-1
// Pb, Pb+1 .. Pc-1
//	...
// Pm, Pm+1, .. Pn-1
//
// The points in the column (e.g., P0, Pa, Pb ...) define the Catmull-Rom spline
// parametrized by the parameter v, while the points in the rows (e.g., P0, P1, ...)
// define the Catmull-Rom spline parametrized by the parameter u. The number of 
// points in the rows may differ, the number of points in the columns is constant.
// The grid is specified as a vtkPolyData where i-th row is described by i-th
// polyline (stored in Lines cells).
//
// .SECTION Caveats
// If you wish to tessellate the spline, use the class vtkParametricFunctionSource.
//
// .SECTION See Also
// vtkParametericCatmullRomSpline vtkPolyData

#ifndef vtkParametricCatmullRomSurface_h
#define vtkParametricCatmullRomSurface_h

#include "vtkParametricCatmullRomSpline.h"
#include "vtkPolyData.h"
#include <vector>

//compatibility for VTK <= 7.0.0
#ifndef VTK_DELETE_FUNCTION
#define VTK_DELETE_FUNCTION =delete
#endif

class vtkParametricCatmullRomSurface : public vtkParametricFunction
{
public:
  vtkTypeMacro(vtkParametricCatmullRomSurface,vtkParametricFunction);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Construct the surface with the following parameters:
  // MinimumU = 0, MaximumU = 1, JoinU = 0 (unless the surface is
  // closed in the parameter u, then JoinU = 1), 
	// MinimumV = 0, MaximumV = 1, JoinU = V (unless the surface is
	// closed in the parameter V, then JoinV = 1), 
	// TwistU = 0, TwistV = 0, DerivativesSupplied = 0,
	// Centripetal splines will be used
  static vtkParametricCatmullRomSurface *New();

  // Description
  // Return the parametric dimension of the class.
  virtual int GetDimension() {return 2;}

	// Description:
	// Catmull-Rom surface.
	//
	// This function performs the mapping \f$f(u,v) \rightarrow (x,y,x)\f$, returning it
	// as Pt. It also returns the partial derivatives Du and Dv.
	// \f$Pt = (x, y, z), Du = (dx/du, dy/du, dz/du), Dv = (dx/dv, dy/dv, dz/dv)\f$ .
	// Then the normal is \f$N = Du X Dv\f$ .
	// N.B. When called from u, v cycles, it is more efficient to have 
	// u in the external cycle and v in the internal one.
	virtual void Evaluate(double uvw[3], double Pt[3], double Duvw[9]);

  // Description:
  // Evaluate a scalar value at parametric coordinate u[0] and Pt[3].
  // The scalar value is just the parameter u[0].
  virtual double EvaluateScalar(double u[3], double Pt[3], double Du[9]);

  // Description:
  // Specify the irregular grid of points defining the surface. 
	// Order of lines cells defines the splines in the parameter u,
	// order of the points in the cell defines the splines in the parameter v.
	vtkSetObjectMacro(PointsGrid, vtkPolyData);
  vtkGetObjectMacro(PointsGrid, vtkPolyData);

  // Description:
  // Controls how the splines are parameterised.
  // Default is Centripetal.
  vtkSetMacro(SplineType, vtkParametricCatmullRomSpline::SplineTypes);
  vtkGetMacro(SplineType, vtkParametricCatmullRomSpline::SplineTypes);

	// Description:
	// Specifies alpha parameter for custom parameterized Catmull-Rom splines.
	vtkSetMacro(Alpha, double);
	vtkGetMacro(Alpha, double);

protected:
  vtkParametricCatmullRomSurface();
  ~vtkParametricCatmullRomSurface();

	// Initializing all the splines in the parameter u
	// Returns zero, if something goes wrong.
	int InitializeSplinesU();

	//Releases the splines in the parameter u.
	void DestroySplinesU();

	//Initializing the spline in the parameter v
	// Returns zero, if something goes wrong.
	int InitializeSplineV(double u);


  // Points definition
  vtkPolyData* PointsGrid;
  
  // Supplemental variables  
  enum vtkParametricCatmullRomSpline::SplineTypes SplineType;
	double Alpha;	//<Alpha parameter to calculate the knot values of points

  
  unsigned long InitializeTime;
	std::vector< vtkParametricCatmullRomSpline* > SplinesU;	//<splines for the paremeter u
  

	double LastU;	//<LastU parameter
	vtkParametricCatmullRomSpline* SplineV;	//<spline in the parameter v
	vtkParametricCatmullRomSpline* SplineDerV;	//<spline in the parameter v for interpolating derivatives in u
	

private:
  vtkParametricCatmullRomSurface(const vtkParametricCatmullRomSurface&) VTK_DELETE_FUNCTION;
  void operator=(const vtkParametricCatmullRomSurface&) VTK_DELETE_FUNCTION;
};

#endif
