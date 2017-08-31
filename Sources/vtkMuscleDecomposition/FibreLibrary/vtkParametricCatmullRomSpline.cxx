/*=========================================================================

	Program:   Visualization Toolkit
	Module:    vtkParametricCatmullRomSpline.cxx

	Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
	All rights reserved.
	See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

		 This software is distributed WITHOUT ANY WARRANTY; without even
		 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
		 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkParametricCatmullRomSpline.h"
#include "vtkObjectFactory.h"
#include "vtkPoints.h"
#include "vtkCatmullRomSpline.h"
#include "vtkMath.h"

vtkStandardNewMacro(vtkParametricCatmullRomSpline);

//----------------------------------------------------------------------------
vtkParametricCatmullRomSpline::vtkParametricCatmullRomSpline()
{
	this->MinimumU = 0;
	this->MaximumU = 1.0;
	this->JoinU = 0;

	this->Points = NULL;

	this->XSpline = vtkCatmullRomSpline::New();
	this->YSpline = vtkCatmullRomSpline::New();
	this->ZSpline = vtkCatmullRomSpline::New();

	this->Closed = 0;
	this->SplineType = vtkParametricCatmullRomSpline::Centripetal;
	this->Alpha = 0.5;	//corresponds to centripetal

	this->DerivativesAvailable = 1;

	this->InitializeTime = 0;
}

//----------------------------------------------------------------------------
vtkParametricCatmullRomSpline::~vtkParametricCatmullRomSpline()
{
	if (this->Points)
	{
		this->Points->Delete();
		this->Points = NULL;
	}
	if (this->XSpline)
	{
		this->XSpline->Delete();
		this->XSpline = NULL;
	}
	if (this->YSpline)
	{
		this->YSpline->Delete();
		this->YSpline = NULL;
	}
	if (this->ZSpline)
	{
		this->ZSpline->Delete();
		this->ZSpline = NULL;
	}
}

//----------------------------------------------------------------------------
void vtkParametricCatmullRomSpline::SetNumberOfPoints(vtkIdType numPts)
{
	if (!this->Points)
	{
		vtkPoints* pts = vtkPoints::New(VTK_DOUBLE);
		this->SetPoints(pts);
		pts->Delete();
	}
	this->Points->SetNumberOfPoints(numPts);
	this->Modified();
}

//----------------------------------------------------------------------------
void vtkParametricCatmullRomSpline::SetPoint(
	vtkIdType index, double x, double y, double z)
{
	if (this->Points)
	{
		this->Points->SetPoint(index, x, y, z);
		this->Modified();
	}
}

//----------------------------------------------------------------------------
void vtkParametricCatmullRomSpline::SetPoints(vtkPoints *pts)
{
	if (pts != this->Points)
	{
		if (this->Points != NULL)
		{
			this->Points->Delete();
		}
		this->Points = pts;
		if (this->Points != NULL)
		{
			this->Points->Register(this);
		}
		this->Modified();
	}
}

//----------------------------------------------------------------------------
void vtkParametricCatmullRomSpline::Evaluate(double U[3], double Pt[3], double DU[9])
{
	// make sure everything has been set up
	if (this->InitializeTime < this->GetMTime())
	{
		if (!this->Initialize())
		{
			return;
		}
	}

	double t = (U[0] < 0.0 ? 0.0 : (U[0] > 1.0 ? 1.0 : U[0]));
	if (this->Closed)
	{
		t *= this->ClosedLength;
	}
	else
	{
		t *= this->Length;
	}

	if (this->Length == 0 && this->Points && this->Points->GetNumberOfPoints() > 0)
	{
		this->Points->GetPoint(0, Pt);
		return;
	}

	// Evaluate the spline at the parameter t
	Pt[0] = this->XSpline->Evaluate(t);
	Pt[1] = this->YSpline->Evaluate(t);
	Pt[2] = this->ZSpline->Evaluate(t);

	if (this->DerivativesAvailable != 0)
	{
		DU[0] = vtkCatmullRomSpline::SafeDownCast(this->XSpline)->EvaluateDerivative(t);
		DU[1] = vtkCatmullRomSpline::SafeDownCast(this->YSpline)->EvaluateDerivative(t);
		DU[2] = vtkCatmullRomSpline::SafeDownCast(this->ZSpline)->EvaluateDerivative(t);
	}
}

//----------------------------------------------------------------------------
double vtkParametricCatmullRomSpline::EvaluateScalar(double u[3], double*, double *)
{
	// make sure everything has been set up
	if (this->InitializeTime < this->GetMTime())
	{
		if (!this->Initialize())
		{
			return 0.0;
		}
	}

	return u[0]; //simply parametric value
}

//----------------------------------------------------------------------------
// Configure the splines for evaluation
int vtkParametricCatmullRomSpline::Initialize()
{
	// Check to make sure splines are available
	if (!this->XSpline || !this->YSpline || !this->ZSpline)
	{
		vtkErrorMacro("Please specify splines");
		return 0;
	}
	if (!this->Points)
	{
		vtkErrorMacro("Please specify points");
		return 0;
	}

	// Make sure that the splines are consistent with this instance
	this->XSpline->SetClosed(this->GetClosed());	
	this->YSpline->SetClosed(this->GetClosed());
	this->ZSpline->SetClosed(this->GetClosed());
	
	vtkIdType npts = this->Points->GetNumberOfPoints();
	if (npts < 1)
	{
		vtkErrorMacro("Please specify at least one point.");
		return 0;
	}

	this->Length = 0;
	this->ClosedLength = 0;

	if (npts < 2)
	{
		// If number of points == 1, then we simply generate a single point.
		return 1;
	}

	this->XSpline->RemoveAllPoints();
	this->YSpline->RemoveAllPoints();
	this->ZSpline->RemoveAllPoints();

	double AB[2][3];	//buffer containing both end-points of the segment
	this->Points->GetPoint(0, AB[0]);

	this->XSpline->AddPoint(0, AB[0][0]);
	this->YSpline->AddPoint(0, AB[0][1]);
	this->ZSpline->AddPoint(0, AB[0][2]);

	if (this->SplineType == vtkParametricCatmullRomSpline::Uniform)
	{		
		//uniform Catmull-Rom spline
		for (int j = 1; j < npts; j++)
		{
			this->Points->GetPoint(j, AB[0]);	//get the next point		
			this->Length += 1.0;
			
			this->XSpline->AddPoint(this->Length, AB[0][0]);
			this->YSpline->AddPoint(this->Length, AB[0][1]);
			this->ZSpline->AddPoint(this->Length, AB[0][2]);		
		}

		if (this->Closed != 0)
			this->ClosedLength = this->Length + 1;
	}
	else
	{
		//non-uniform Catmull-Rom spline


	double alfa;
	switch (this->SplineType)
	{
	case vtkParametricCatmullRomSpline::Chordal:
		alfa = 1; 
		break;
	case vtkParametricCatmullRomSpline::Centripetal:
		alfa = 0.5; 
		break;
	default:
		alfa = this->Alpha;
		break;
	}

	bool len1stOK = false;
	double lenFirst, lenLast;	//first segment distance

	int iFirstPt = 0;
	for (int j = 1; j < npts; j++)
	{
		this->Points->GetPoint(j, AB[1 - iFirstPt]);	//get the next point		

		double len = pow(sqrt(vtkMath::Distance2BetweenPoints(AB[0], AB[1])), alfa);
		if (len == 0.0)
			continue;	//skip redundant points as they would cause troubles

		if (!len1stOK) {
			lenFirst = len; len1stOK = true;
		}
		lenLast = len;

		this->Length += len;

		this->XSpline->AddPoint(this->Length, AB[1 - iFirstPt][0]);
		this->YSpline->AddPoint(this->Length, AB[1 - iFirstPt][1]);
		this->ZSpline->AddPoint(this->Length, AB[1 - iFirstPt][2]);

		iFirstPt = 1 - iFirstPt;
	}

	if (this->Closed == 0)
	{
		//open spline
		this->XSpline->SetParametricRange(0.0, this->Length);
		this->YSpline->SetParametricRange(0.0, this->Length);
		this->ZSpline->SetParametricRange(0.0, this->Length);

		//automatically calculated constraints
		this->XSpline->SetLeftConstraint(0);
		this->YSpline->SetLeftConstraint(0);
		this->ZSpline->SetLeftConstraint(0);

		this->XSpline->SetRightConstraint(0);
		this->YSpline->SetRightConstraint(0);
		this->ZSpline->SetRightConstraint(0);
		}
	else
		{
			this->Points->GetPoint(0, AB[1 - iFirstPt]); //get the first point

			double len = pow(sqrt(vtkMath::Distance2BetweenPoints(AB[0], AB[1])), alfa);
			if (len == 0.0)
			{
				//last point coincides with the first one
				this->XSpline->RemovePoint(this->Length);
				this->YSpline->RemovePoint(this->Length);
				this->ZSpline->RemovePoint(this->Length);
				
				this->ClosedLength = this->Length;
			}
			else
			{
				this->ClosedLength = this->Length + len;	
				lenLast = len;
			}

			//constraints are specified
			this->XSpline->SetParametricRange(0.0, this->ClosedLength);
			this->YSpline->SetParametricRange(0.0, this->ClosedLength);
			this->ZSpline->SetParametricRange(0.0, this->ClosedLength);

			this->XSpline->SetLeftConstraint(1);
			this->XSpline->SetLeftValue(-lenLast);

			this->YSpline->SetLeftConstraint(1);
			this->YSpline->SetLeftValue(-lenLast);

			this->ZSpline->SetLeftConstraint(1);
			this->ZSpline->SetLeftValue(-lenLast);

			this->XSpline->SetRightConstraint(1);
			this->XSpline->SetRightValue(this->ClosedLength + lenFirst);

			this->YSpline->SetRightConstraint(1);
			this->YSpline->SetRightValue(this->ClosedLength + lenFirst);

			this->ZSpline->SetRightConstraint(1);
			this->ZSpline->SetRightValue(this->ClosedLength + lenFirst);
		}
	}

	this->InitializeTime = this->GetMTime();
	return 1;
}


//----------------------------------------------------------------------------
void vtkParametricCatmullRomSpline::PrintSelf(ostream& os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os, indent);

	os << indent << "Points: ";
	if (this->Points)
	{
		os << this->Points << "\n";
	}
	else
	{
		os << "(none)\n";
	}

	os << indent << "X Spline: ";
	if (this->XSpline)
	{
		os << this->XSpline << "\n";
	}
	else
	{
		os << "(none)\n";
	}
	os << indent << "Y Spline: ";
	if (this->YSpline)
	{
		os << this->YSpline << "\n";
	}
	else
	{
		os << "(none)\n";
	}
	os << indent << "Z Spline: ";
	if (this->ZSpline)
	{
		os << this->ZSpline << "\n";
	}
	else
	{
		os << "(none)\n";
	}

	os << indent << "Closed: " << (this->Closed ? "On\n" : "Off\n");	
	os << indent << "Spline Type: " << this->SplineType << "\n";
}
