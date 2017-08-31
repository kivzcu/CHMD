/*=========================================================================

	Program:   Visualization Toolkit
	Module:    vtkParametricCatmullRomSurface.cxx

	Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
	All rights reserved.
	See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

		 This software is distributed WITHOUT ANY WARRANTY; without even
		 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
		 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkParametricCatmullRomSurface.h"
#include "vtkObjectFactory.h"
#include "vtkPolyData.h"
#include "vtkMath.h"

vtkStandardNewMacro(vtkParametricCatmullRomSurface);

//----------------------------------------------------------------------------
vtkParametricCatmullRomSurface::vtkParametricCatmullRomSurface()
{
	this->MinimumU = this->MinimumV = 0;
	this->MaximumU = this->MaximumV = 1.0;
	this->JoinU = this->JoinV = 0;

	this->PointsGrid = NULL;
	
	this->SplineType = vtkParametricCatmullRomSpline::Centripetal;
	this->Alpha = 0.5;	//corresponds to centripetal

	this->SplineV = NULL;
	this->SplineDerV = NULL;
	this->LastU = 0.0;

	this->DerivativesAvailable = 1;

	this->InitializeTime = 0;
}

//----------------------------------------------------------------------------
vtkParametricCatmullRomSurface::~vtkParametricCatmullRomSurface()
{
	if (this->PointsGrid != NULL)
	{
		this->PointsGrid->Delete();
		this->PointsGrid = NULL;
	}	

	if (this->SplineV != NULL)
	{
		this->SplineV->Delete();
		this->SplineV = NULL;
	}

	if (this->SplineDerV != NULL)
	{
		this->SplineDerV->Delete();
		this->SplineDerV = NULL;
	}

	DestroySplinesU();
}

//----------------------------------------------------------------------------
void vtkParametricCatmullRomSurface::Evaluate(double uvw[3], double Pt[3], double Duvw[9])
{
	// make sure everything has been set up
	if (this->InitializeTime < this->GetMTime())
	{
		if (!this->InitializeSplinesU())
		{
			return;
		}
	}
	
	//clamp u, v parameters to the range
	double u = (uvw[0] < this->MinimumU ? this->MinimumU : 
		(uvw[0] > this->MaximumU ? this->MaximumU : uvw[0]));
	double v = (uvw[1] < this->MinimumV ? this->MinimumV : 
		(uvw[1] > this->MaximumV ? this->MaximumV : uvw[1]));

	if (this->SplineV == NULL || this->LastU != u)
	{
		if (!this->InitializeSplineV(u))
		{
			return;
		}
	}

	if (this->DerivativesAvailable != 0)
		this->SplineDerV->Evaluate(&uvw[1], Duvw, Duvw);	//derivative in U
	
	this->SplineV->Evaluate(&uvw[1], Pt, &Duvw[3]);		
}

//----------------------------------------------------------------------------
double vtkParametricCatmullRomSurface::EvaluateScalar(double*, double*, double *)
{	
	return 0;
}

//----------------------------------------------------------------------------
// Configure the splines for evaluation
int vtkParametricCatmullRomSurface::InitializeSplinesU()
{
	// Check to make sure we have valid input
	if (this->PointsGrid == NULL)
	{
		vtkErrorMacro("Please specify the grid of points.");
		return 0;
	}

	DestroySplinesU();
		
	vtkIdType nCells = this->PointsGrid->GetNumberOfCells();
	if (nCells == 0 || this->PointsGrid->GetLines() == NULL)
	{
		vtkErrorMacro("Please specify at least one point.");
		return 0;
	}
	
	//build cells, if not already built
	this->PointsGrid->GetCellType(0);
	this->SplinesU.reserve(nCells);

	for (int i = 0; i < nCells; i++)
	{
		vtkIdType nPts, *pPts;
		this->PointsGrid->GetCellPoints(i, nPts, pPts);

		vtkParametricCatmullRomSpline* spline = vtkParametricCatmullRomSpline::New();
		spline->SetSplineType(this->GetSplineType());
		spline->SetAlpha(this->GetAlpha());
		spline->SetClosed(this->JoinU);	//set wheter it is closed

		//copy the points
		spline->SetNumberOfPoints(nPts);
		for (int j = 0; j < nPts; j++)
		{
			const double* xyz = this->PointsGrid->GetPoint(pPts[j]);
			spline->SetPoint(j, xyz[0], xyz[1], xyz[2]);
		}

		this->SplinesU.push_back(spline);
	}

	this->InitializeTime = this->GetMTime();
	return 1;
}

//----------------------------------------------------------------------------
int vtkParametricCatmullRomSurface::InitializeSplineV(double u)
{
	double uvw[3], xyz[3], Duvw[9];
	uvw[0] = u;

	int sz = (int)this->SplinesU.size();
	if (sz == 0) {
		return 0;	//should not happen
	}

	//initialize spline v
	if (this->SplineV != NULL) {
		this->SplineV->Delete();
	}

	if (this->SplineDerV != NULL) {
		this->SplineDerV->Delete();
	}

	this->SplineV = vtkParametricCatmullRomSpline::New();
	this->SplineV->SetSplineType(this->GetSplineType());
	this->SplineV->SetAlpha(this->GetAlpha());
	this->SplineV->SetClosed(this->JoinV);
	this->SplineV->SetDerivativesAvailable(this->GetDerivativesAvailable());

	if (this->DerivativesAvailable)
	{
		this->SplineDerV = vtkParametricCatmullRomSpline::New();
		this->SplineDerV->SetSplineType(this->GetSplineType());
		this->SplineDerV->SetAlpha(this->GetAlpha());		
		this->SplineDerV->SetDerivativesAvailable(0);
		this->SplineDerV->SetNumberOfPoints(sz);
	}
	
	this->SplineV->SetNumberOfPoints(sz);
	for (int i = 0; i < sz; i++)
	{		
		this->SplinesU[i]->Evaluate(uvw, xyz, Duvw);
		this->SplineV->SetPoint(i, xyz[0], xyz[1], xyz[2]);

		if (this->SplineDerV != NULL)
			this->SplineDerV->SetPoint(i, Duvw[0], Duvw[1], Duvw[2]);
	}

	this->LastU = u;
	return 1;
}


//----------------------------------------------------------------------------
void vtkParametricCatmullRomSurface::DestroySplinesU()
{
	auto size = this->SplinesU.size();
	for (size_t i = 0; i < size; i++)
	{
		this->SplinesU[i]->Delete();
		this->SplinesU[i] = NULL;
	}

	this->SplinesU.clear();
}

//----------------------------------------------------------------------------
void vtkParametricCatmullRomSurface::PrintSelf(ostream& os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os, indent);

	os << indent << "PointsGrid: ";
	if (this->PointsGrid)
	{
		os << this->PointsGrid << "\n";
	}
	else
	{
		os << "(none)\n";
	}

	os << indent << "Spline Type: " << this->SplineType << "\n";
}
