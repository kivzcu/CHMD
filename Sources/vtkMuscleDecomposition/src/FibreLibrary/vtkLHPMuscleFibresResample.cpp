/*=========================================================================
Program: Musculoskeletal Modeling (VPHOP WP10)
Module: vtkLHPMuscleFibresResample.cpp

Authors: Josef Kohout
==========================================================================
Copyright (c) 2014 University of West Bohemia (www.zcu.cz)
See the COPYINGS file for license details
=========================================================================
*/
#include "vtkLHPMuscleFibresResample.h"
#include "vtkObjectFactory.h"
#include "vtkPolyData.h"
#include "vtkCellArray.h"
#include "vtkMath.h"
#include "vtkLHPMuscleFibresMath.h"
#include "vtkDoubleArray.h"
#include "vtkCommand.h"
#include <unordered_map>
#include <vtkPointData.h>
#include <vtkUnsignedCharArray.h>
using namespace std;
static unordered_map<vtkIdType,vtkIdType> map;

vtkStandardNewMacro(vtkLHPMuscleFibresResample);


vtkLHPMuscleFibresResample::vtkLHPMuscleFibresResample()
{
	this->UseCentripetalCatmullRom = 1;	//ON
	this->Precision = 0.1;		//0.1 mm by the default
	this->Resolution = 0;		//auto
	this->PrecisionMode = 0;	//use Resolution by default
}


vtkLHPMuscleFibresResample::~vtkLHPMuscleFibresResample()
{
	
}

#if VTK_MAJOR_VERSION < 5
#define VTK6_RESULT_OK
#define VTK6_RESULT_FAIL
/**
Executes the data operation.

@param [in,out]	output	If non-null, the output.
*/
/*virtual*/ void vtkLHPMuscleFibresResample::ExecuteData(vtkDataObject *output)
{
	vtkPolyData* inputPoly = GetInput();
	vtkPolyData* outputPoly = vtkPolyData::SafeDownCast(output);
#else

#define VTK6_RESULT_OK		1
#define VTK6_RESULT_FAIL	0

// This is the method of the algorithm in which the algorithm should fill in the output ports
/*virtual*/ int vtkLHPMuscleFibresResample::RequestData(vtkInformation* request,
	vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
	vtkPolyData* inputPoly = vtkPolyData::GetData(inputVector[0]);
	vtkPolyData* outputPoly = vtkPolyData::GetData(outputVector);
#endif

	//check whether output is valid	
	if (inputPoly == NULL || inputPoly->GetPoints() == NULL ||
		inputPoly->GetPoints()->GetNumberOfPoints() == 0)
	{
		vtkErrorMacro(<< "Invalid input for vtkLHPMuscleFibresAnalysisFilter.");
		return VTK6_RESULT_FAIL;   //we have no valid input
	}
	
	if (outputPoly == NULL)
	{
		vtkErrorMacro(<< "Invalid output for vtkLHPMuscleFibresAnalysisFilter.");
		return VTK6_RESULT_FAIL;   //we have no valid output
	}

	//Force calling inputPoly->BuildCells, if it has not been already built 	
	inputPoly->GetCellType(0);

	//calculate the Catmul-Rom parameterization
	vtkDoubleArray* segLengths = vtkDoubleArray::New();
	vtkLHPMuscleFibresMath::CalculateFibresSegmentsLengths(inputPoly, segLengths);

	vtkDoubleArray* params = vtkDoubleArray::New();
	vtkLHPMuscleFibresMath::CalculateFibresParameterization(inputPoly, segLengths, this->UseCentripetalCatmullRom != 0, params);

	vtkPoints* points = inputPoly->GetPoints()->NewInstance();
	vtkCellArray* cells = vtkCellArray::New();

	//and now process every fibres	
	int nFibres = inputPoly->GetNumberOfCells();

	//determine the resolution
	int res = this->Resolution;
	if (this->PrecisionMode == 0 && res == 0)
	{
		for (int i = 0; i < nFibres; i++)
		{
			vtkIdType nPts, *pPts;
			inputPoly->GetCellPoints(i, nPts, pPts);
			if (nPts > res) {
				res = nPts;
			}
		}
	}

	//process the fibres
	const double* pParams = params->GetPointer(0);
	for (int i = 0, index = 0; i < nFibres; i++)
	{
		vtkIdType nPts, *pPts;
		inputPoly->GetCellPoints(i, nPts, pPts);		
		
		ResampleFibre(nPts, pPts, inputPoly->GetPoints(), pParams,
			segLengths->GetPointer(index - i), points, cells, res);

		index += nPts;
	}

	outputPoly->SetPoints(points);
	outputPoly->SetLines(cells);

	cells->Delete();
	points->Delete();

	params->Delete();
	segLengths->Delete();
	return VTK6_RESULT_OK;
}

//Copies the fibre defined by its sequence of points to the output.
void vtkLHPMuscleFibresResample::CopyFibre(vtkIdType nPts, const vtkIdType* pPts,
	const vtkPoints* coords, vtkPoints* points, vtkCellArray* cells)
{
	vtkIdList* idsList = vtkIdList::New();
	idsList->Allocate(nPts);	//reserve the buffer	

	for (vtkIdType i = 0; i < nPts; i++)
	{
		idsList->InsertNextId(points->InsertNextPoint(
			const_cast<vtkPoints*>(coords)->GetPoint(pPts[i])));
	}

	cells->InsertNextCell(idsList);
	idsList->Delete();
}

//Performs the resampling of a fibre defined by its sequence of points.
/*virtual*/ void vtkLHPMuscleFibresResample::ResampleFibre(vtkIdType nPts, const vtkIdType* pPts, const vtkPoints* coords,
	const double* params, const double* segLengths, vtkPoints* points, vtkCellArray* cells, int resolution)
{
	double* tp = NULL;
	int tpsize = 0;	

	if (this->PrecisionMode == 0)
	{
		//check if the resolution is invalid
		if (resolution < 2)
			resolution = 2;

		//if the fibre has the requested resolution, we are ready
		if (resolution == nPts) {
			CopyFibre(nPts, pPts, coords, points, cells);
			return;	//we are ready
		}
	
		//otherwise we need to create sampling plan
		double dt = (params[pPts[nPts - 1]] - params[pPts[0]]) / (resolution - 1);
		
		tp = new double[tpsize = resolution];
		tp[0] = params[0];

		for (int i = 1; i < resolution; i++) {
			tp[i] = tp[0] + dt*i;
		}

		tp[tpsize - 1] = params[pPts[nPts - 1]];
	}

	double limit2 = this->Precision*this->Precision;

	vtkIdList* idsList = vtkIdList::New();
	idsList->Allocate(nPts);	//reserve the buffer	
	

	//Let the fibre be formed by points Q0, Q1, ... Qn-1 having Catmull-Rom parameterization Q0_t, Q1_t, ...Qn-1_t
	//Cubic Catmul-Rom curve C is defined by 4 points P0..P3 and their knots (parameterization) t0..t3
	//where knots ti+1 are defined as ti+1 = ti + ||Pi+1 - Pi||^alpha 
	//whereby alpha = 1/2 for centripetal CR, 1 for chordal
	//see: http://en.wikipedia.org/wiki/Centripetal_Catmull%E2%80%93Rom_spline
	//Note: the curve is valid only from P1 to P2 (points P0 and P3 are used to define it), i.e.,
	//the valid range of the parameter t is <t1...t2>
	
	double P[4][3];	//buffer for the points P0..P3
	double t[4];		//buffer for the knots values t0..t3

	const_cast<vtkPoints*>(coords)->GetPoint(pPts[0], P[1]); t[1] = params[pPts[0]];	//P1 is the first point of the fibre
	const_cast<vtkPoints*>(coords)->GetPoint(pPts[1], P[2]); t[2] = params[pPts[1]];	//P2 is the second point of the fibre

	//construct an artificial point P0 to define the Catmull-Rom cubic spline in the
	//in the opposite direction from Q1 as Q0A = 2*Q0 - Q1
	for (int k = 0; k < 3; k++) {
		P[0][k] = 2 * P[1][k] - P[2][k];
	}

	t[0] = t[1] - t[2];

	int tp_index = 0;
	for (int j = 0; j < nPts - 1; j++)
	{
		//load the fourth point			
		int iv[4] = { (j + 0) % 4, (j + 1) % 4, (j + 2) % 4, (j + 3) % 4, };

		if (j + 2 < nPts)
		{
			const_cast<vtkPoints*>(coords)->GetPoint(pPts[j + 2], P[iv[3]]);
			t[iv[3]] = params[pPts[j + 2]];
		}
		else
		{
			//construct an artificial point P0 to define the Catmull-Rom cubic spline in the
			//in the opposite direction from Q2 as Q3A = 2*Q2 + Q1
			for (int k = 0; k < 3; k++) {
				P[iv[3]][k] = 2 * P[iv[2]][k] - P[iv[1]][k];
			}

			t[iv[3]] = 2 * t[iv[2]] - t[iv[1]];			
		}

		//OK, we have here the full definition of the Catmull-Rom cubic spline in P0..P3 with knots in t0..t3
		//and the valid curve segment is from P1 to P2
		if (this->PrecisionMode == 0)
		{
			while (tp_index < tpsize && tp[tp_index] < t[iv[1]]) {
				++tp_index;
			}

			int tp_nexindex = tp_index;
			while (tp_nexindex < tpsize && tp[tp_nexindex] < t[iv[2]]) {
				++tp_nexindex;	//move to next
			}
				
			if (tp_nexindex > tp_index)
			{
				//there is at least one sample to get
				SampleCurve(P, t, iv, &tp[tp_index], tp_nexindex - tp_index, points, idsList);
				tp_index = tp_nexindex;
			}
		}
		else
		{
			//insert the point P1 into the output
			vtkIdType id = points->InsertNextPoint(P[iv[1]]);
			idsList->InsertNextId(id);

			double dist = vtkMath::Distance2BetweenPoints(P[iv[1]], P[iv[2]]);
			if (dist > limit2)
			{
				int nSteps = ((int)sqrt(dist / limit2)) + 1;
				if (nSteps + 1 > tpsize)
				{
					delete[] tp;	//resize the buffer
					tp = new double[tpsize = nSteps + 1];
				}

				double delta_t = (t[iv[2]] - t[iv[1]]) / (nSteps + 1);
				for (int i = 1; i <= nSteps; i++){
					tp[i] = t[iv[1]] + i*delta_t;
				}

				//there is at least one sample to get
				SampleCurve(P, t, iv, &tp[1], nSteps, points, idsList);
			}
		}
	} //end for

	//insert the last point	
	vtkIdType id = points->InsertNextPoint(P[nPts % 4]);
	idsList->InsertNextId(id);

//done.
	
	delete[] tp;

	cells->InsertNextCell(idsList);
	idsList->Delete();
}

//Performs the sampling of the curve
/*virtual*/ void vtkLHPMuscleFibresResample::SampleCurve(double P[4][3], double t[4], int iv[4],
	double* tp, int count, vtkPoints* points, vtkIdList* idsList)
{
	//we need to introduce some intermediate points

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
	//M(1
	// It may be seen that these coefficients have the form nominator / denominator and the terms in the nominators/denominators are shared between them.
	//So after some tuning we get the code below. Although not so nice as the recursive formula, it requires "only" 48(+), 21(-), 91 (*) and 18 (/)
	//floating point operations in the pre-processing and afterwards s*9(*) and s*9*(+) where s is the number of steps, whilst the recursive formula
	//requires s*18(+), 4(-) + s*36(-), 18(*)+s*36(*), 18(/) + s*18(/), thus the non-recursive needs less operations with s >= 2. However,			
	//considering the fact that division (FDIV) in double is ten times slower than multiplication, which itself is twice slower than +/- operations (in year 2015)
	//even for a single point generation, the non-recursive formula is 1.39x faster (and indeed, the speed-up grows with s, for s being 10, it is 5.97)

	double dt01 = t[iv[0]] - t[iv[1]], dt02 = t[iv[0]] - t[iv[2]], dt03 = t[iv[0]] - t[iv[3]],
		dt12 = t[iv[1]] - t[iv[2]], dt13 = t[iv[1]] - t[iv[3]],
		dt23 = t[iv[2]] - t[iv[3]];

	double st11 = t[iv[1]] * t[iv[1]], st22 = t[iv[2]] * t[iv[2]];

	double a5 = -(t[iv[1]] + 2 * t[iv[2]]);
	double a6 = 2 * t[iv[1]] + t[iv[2]];

	double t01 = t[iv[0]] * t[iv[1]], t12 = t[iv[1]] * t[iv[2]];
	double denom = dt01*dt02*dt12;

	typedef double Vect4[4];
	Vect4 M[4];
	M[0][0] = 1 / denom; M[0][1] = a5 / denom; M[0][2] = t[iv[2]] * a6 / denom; M[0][3] = -t[iv[1]] * st22 / denom;

	denom = dt01*dt12*dt12*dt13;
	M[1][0] = -dt03 / denom;
	M[1][1] = (t[iv[0]] * a6 + t[iv[1]] * dt23 - 2 * t[iv[2]] * t[iv[3]] - st11) / denom;
	M[1][2] = (t12*(dt12 - dt03 - 2 * t[iv[0]]) - t[iv[3]] * (t[iv[0]] * dt12 - st11 - st22)) / denom;
	M[1][3] = (t12 * t[iv[3]] * dt01 + t[iv[0]] * st22 * dt13) / denom;

	denom = dt02*dt12*dt12*dt23;
	M[2][0] = dt03 / denom;
	M[2][1] = (-(t[iv[2]] * (dt01 - t[iv[3]] + dt23) + t[iv[1]] * (dt03 + t[iv[0]]))) / denom;
	M[2][2] = (t01 * (t[iv[1]] + t[iv[2]] + t[iv[3]]) + t[iv[0]] * t[iv[2]] * dt23 - t12*(dt12 + 3 * t[iv[3]])) / denom;
	M[2][3] = -(t01 * t[iv[2]] * dt23 + st11 * t[iv[3]] * dt02) / denom;

	denom = dt12*dt13*dt23;
	M[3][0] = -1 / denom; M[3][1] = a6 / denom; M[3][2] = t[iv[1]] * a5 / denom; M[3][3] = st11*t[iv[2]] / denom;

	//t0*(2*t1+t2)+t1*(t2-t3) -2*t2*t3 - t1^2
	//t1*t2*((t1-t2)-(t0-t3)-2*t0)-t3*(t0*(t1-t2)-t1^2-t2^2)
	//t2*(t1*t3*(t0-t1)+t0*t2*(t1-t3))

	//-(t2*((t0-t1)-t3+(t2-t3))+t1*((t0-t3)+t0))
	//t0*t1*(t1+t2+t3)+t0*t2*(t2-t3)-t1*t2*((t1-t2)+3*t3)
	//-t1*(t0*t2*(t2-t3)+t1*t3*(t0-t2))

	double vP[4][3]; //vP = [P0 P1 P2 P3]*M
	for (int k = 0; k < 3; k++)
	{
		for (int j = 0; j < 4; j++) {
			vP[j][k] = P[iv[0]][k] * M[0][j] + P[iv[1]][k] * M[1][j] + P[iv[2]][k] * M[2][j] + P[iv[3]][k] * M[3][j];
		}
	}
	
	for (int i = 0; i < count; i++)
	{
		double curr_t = tp[i];

		double newPoint[3];
		for (int k = 0; k < 3; k++) {	//Horner
			newPoint[k] = vP[3][k] + curr_t*(vP[2][k] + curr_t*(vP[1][k] + curr_t*vP[0][k]));
		}

		//insert the new point
		//idsList->InsertNextId(points->InsertNextPoint(newPoint));
		vtkIdType id = points->InsertNextPoint(newPoint);
		idsList->InsertNextId(id);
	}
}