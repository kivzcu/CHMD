/*=========================================================================
Program: Musculoskeletal Modeling (VPHOP WP10)
Module: vtkLHPMuscleFibresAnalysisFilter.cpp

Authors: Josef Kohout
==========================================================================
Copyright (c) 2014 University of West Bohemia (www.zcu.cz)
See the COPYINGS file for license details
=========================================================================
*/

#include "vtkLHPMuscleFibresAnalysisFilter.h"
#include "vtkObjectFactory.h"
#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkCellArray.h"
#include "vtkMath.h"
#include "vtkSmartPointer.h"
#include "vtkDoubleArray.h"
#include "vtkIntArray.h"
#include "vtkLHPMuscleFibresMath.h"
//#include "mafDbg.h"

#include <float.h>
#include <assert.h>

vtkStandardNewMacro(vtkLHPMuscleFibresAnalysisFilter);



/*static*/ const char* vtkLHPMuscleFibresAnalysisFilter::FibresLengthsFieldName = "FibresLenghts";
/*static*/ const char* vtkLHPMuscleFibresAnalysisFilter::FibresProximalPAFieldName = "FibresProximalPA";
/*static*/ const char* vtkLHPMuscleFibresAnalysisFilter::FibresDistalPAFieldName = "FibresDistalPA";

/*static*/ const char* vtkLHPMuscleFibresAnalysisFilter::FibresChordalParameterizationFieldName = "FibresChordalParameterization";
/*static*/ const char* vtkLHPMuscleFibresAnalysisFilter::FibresCentripetalParameterizationFieldName = "FibresCentripetalParameterization";
/*static*/ const char* vtkLHPMuscleFibresAnalysisFilter::FibresTangentVectorsFieldName = "FibresTangentVectors";

/*static*/ const char* vtkLHPMuscleFibresAnalysisFilter::FibresSegmentsCountFieldName = "FibresSegmentsCount";
/*static*/ const char* vtkLHPMuscleFibresAnalysisFilter::FibresSegmentsLengthsFieldName = "FibresSegmentsLengths";
/*static*/ const char* vtkLHPMuscleFibresAnalysisFilter::FibresProximalTangentVectorsFieldName = "FibresProximalTangentVectors";
/*static*/ const char* vtkLHPMuscleFibresAnalysisFilter::FibresDistalTangentVectorsFieldName = "FibresDistalTangentVectors";

/*static*/ const char* vtkLHPMuscleFibresAnalysisFilter::FibresRadiiFieldName = "FibresRadii";
/*static*/ const char* vtkLHPMuscleFibresAnalysisFilter::FibresVolumeFieldName = "FibresVolume";

/**
 Default constructor.
 */
vtkLHPMuscleFibresAnalysisFilter::vtkLHPMuscleFibresAnalysisFilter()
{
	this->ComputeFibresLengths = 0;
	this->ComputeProximalLineOfAction = 0;
	this->ComputeDistalLineOfAction = 0;
	this->ComputeProximalPennationAngle = 0;
	this->ComputeDistalPennationAngle = 0;
	this->ComputeFibresParameterization = 0;
	this->ComputeFibresTangentVector = 0;
	this->ComputeFibresVolumes = 0;
	this->ComputeTotalFibresVolume = 0;
	this->ComputePCSA = 0;

	this->UseCentripetalCatmullRom = 1;
	this->FibresVolumeCalculationMethod = VOLUME_CALCULATION_RAVICHANDIRAN;

	m_LastExecuteInputTimeStamp = 0;	//the input does not changed
	m_SettingsChangedFlags = 0;
	m_CachedDataFlags = 0;

	this->m_FibresSegmentsCount = NULL;
	this->m_FibresSegmentsLengths = NULL;
	this->m_FibresLengths = NULL;
	this->m_FibresParameterization = NULL;
	this->m_FibresCRParameterization = NULL;
	this->m_FibresTangentVectors = NULL;
	this->m_FibresProximalAvgTgVec = NULL;
	this->m_FibresDistalAvgTgVec = NULL;
	this->m_FibresProximalPA = NULL;
	this->m_FibresDistalPA = NULL;
	this->m_FibresRadii = NULL;
	this->m_FibresVolume = NULL;

	//mark lines of action as invalid
	this->ProximalLineOfAction[0] = 0.0;
#pragma warning (suppress: 4723)
	this->DistalLineOfAction[0] = 1.0 / this->ProximalLineOfAction[0];
	this->ProximalLineOfAction[0] = this->DistalLineOfAction[0];

	//"approximately 15 − 20% of the entire fascicle length is included in each  of the proximal and distal regions"
	//see Lee et al., 2014
	this->DistalThreshold = this->ProximalThreshold = 0.15;
}

#ifndef vtkDEL
#define vtkDEL(x) if ((x) != NULL) { (x)->Delete(); (x) = NULL; }
#endif

vtkLHPMuscleFibresAnalysisFilter::~vtkLHPMuscleFibresAnalysisFilter()
{
	vtkDEL(this->m_FibresSegmentsCount);
	vtkDEL(this->m_FibresSegmentsLengths);
	vtkDEL(this->m_FibresLengths);
	vtkDEL(this->m_FibresParameterization);
	vtkDEL(this->m_FibresCRParameterization);
	vtkDEL(this->m_FibresTangentVectors);
	vtkDEL(this->m_FibresProximalAvgTgVec);
	vtkDEL(this->m_FibresDistalAvgTgVec);
	vtkDEL(this->m_FibresProximalPA);
	vtkDEL(this->m_FibresDistalPA);
	vtkDEL(this->m_FibresRadii);
	vtkDEL(this->m_FibresVolume);
}

/**
Detect changes.
*/
/* virtual*/ void vtkLHPMuscleFibresAnalysisFilter::MustUpdateCheck(vtkPolyData* inputPoly)
{	 
	if (inputPoly->GetMTime() > this->m_LastExecuteInputTimeStamp) {
		this->m_CachedDataFlags = COMPUTE_NOTHING;	//nothing is cached
	}
	else
	{
		//if UseCentripetalCatmullRom value has changed =>
		//all TANGENTVECTORS, LINEOFACTION, PENNATIONANGLE, PCSA must be recalculated
		if ((m_SettingsChangedFlags & USECENTRIPETALCATMULLROM_CHANGED) != 0)
			this->m_CachedDataFlags &= ~(COMPUTE_FIBRES_TANGENTVECTORS | COMPUTE_FIBRES_PROXIMAL_TANGENTVECTORS |
			COMPUTE_FIBRES_DISTAL_TANGENTVECTORS | COMPUTE_PROXIMAL_LINEOFACTION | COMPUTE_DISTAL_LINEOFACTION |
			COMPUTE_FIBRES_PROXIMAL_PENNATIONANGLE | COMPUTE_FIBRES_DISTAL_PENNATIONANGLE | COMPUTE_PCSA);

		//if ProximalThreshold value has changed => 
		if ((m_SettingsChangedFlags & PROXIMALTHRESHOLD_CHANGED) != 0)
			this->m_CachedDataFlags &= ~(COMPUTE_FIBRES_PROXIMAL_TANGENTVECTORS |
			COMPUTE_PROXIMAL_LINEOFACTION | COMPUTE_FIBRES_PROXIMAL_PENNATIONANGLE | COMPUTE_PCSA);

		//if DistalThreshold value has changed => 
		if ((m_SettingsChangedFlags & DISTALTHRESHOLD_CHANGED) != 0)
			this->m_CachedDataFlags &= ~(COMPUTE_FIBRES_DISTAL_TANGENTVECTORS |
			COMPUTE_DISTAL_LINEOFACTION | COMPUTE_FIBRES_DISTAL_PENNATIONANGLE | COMPUTE_PCSA);

		if ((m_SettingsChangedFlags & VOLUME_CALCULATION_METHOD_CHANGED) != 0)
			this->m_CachedDataFlags &= ~(COMPUTE_FIBRES_VOLUME |
			COMPUTE_TOTAL_FIBRES_VOLUME | COMPUTE_PCSA);
	}

	m_SettingsChangedFlags = 0;	//the change accepted

	//if line of actions have been specified, do not calculate them
	if (this->ComputeProximalLineOfAction == 0 && std::isfinite(this->ProximalLineOfAction[0]))
		this->m_CachedDataFlags |= COMPUTE_PROXIMAL_LINEOFACTION;

	if (this->ComputeDistalLineOfAction == 0 && std::isfinite(this->DistalLineOfAction[0]))
		this->m_CachedDataFlags |= COMPUTE_DISTAL_LINEOFACTION;


	//Call inputPoly->BuildCells, if it has not been already built 
	//Cells are required to be built prior to calling time-efficient routines such that GetCellPoints
	//However, calling BuildCells recreates the structure always no matter if it exists or not
	//Since there is no way to find out if the structure has been built from external classes,
	//we need to call one of safe routines that automatically builds cells if they are not present
	inputPoly->GetCellType(0);

	this->m_LastExecuteInputTimeStamp = inputPoly->GetMTime();
}

#if VTK_MAJOR_VERSION < 5
#define VTK6_RESULT_OK
#define VTK6_RESULT_FAIL
/**
Executes the data operation.

@param [in,out]	output	If non-null, the output.
*/
/*virtual*/ void vtkLHPMuscleFibresAnalysisFilter::ExecuteData(vtkDataObject *output)
{
	vtkPolyData* inputPoly = GetInput();
	vtkPolyData* outputPoly = vtkPolyData::SafeDownCast(output);
#else

#define VTK6_RESULT_OK		1
#define VTK6_RESULT_FAIL	0

// This is the method of the algorithm in which the algorithm should fill in the output ports
/*virtual*/ int vtkLHPMuscleFibresAnalysisFilter::RequestData(vtkInformation* request,
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

	//Pass the geometry and topology from the input to output 
	outputPoly->CopyStructure(inputPoly);
	
	//Pass the point data and cell data from the input
	vtkPointData* outPD = outputPoly->GetPointData();
	vtkPointData* inPD = inputPoly->GetPointData();
	if (inPD != NULL && outPD != NULL) {		
		outPD->PassData(inPD);
	}

	vtkCellData* outCD = outputPoly->GetCellData();
	vtkCellData* inCD = inputPoly->GetCellData();
	if (inCD != NULL && outCD != NULL) {
		outCD->PassData(inCD);
	}

	//now, detect what must be recalculated
	MustUpdateCheck(inputPoly);

	if (this->ComputeFibresLengths != 0) {
		outCD->AddArray(this->GetFibresLengths());
	}

	if (this->ComputeFibresParameterization != 0) {
		outPD->AddArray(this->GetFibresParameterization());
	}

	if (this->ComputeFibresTangentVector != 0) {
		outPD->AddArray(this->GetFibresTangentVectors());
	}

	if (this->ComputeFibresRadii != 0) {
		outPD->AddArray(this->GetFibresRadii());
	}

	if (this->ComputeProximalLineOfAction != 0) {
		GetProximalLineOfAction();	//force its calculation
	}

	if (this->ComputeDistalLineOfAction != 0) {
		GetDistalLineOfAction();	//force its calculation
	}

	if (this->ComputeProximalPennationAngle != 0) {
		outCD->AddArray(this->GetFibresProximalPA());
	}

	if (this->ComputeDistalPennationAngle != 0) {
		outCD->AddArray(this->GetFibresDistalPA());
	}

	if (this->ComputeFibresVolumes != 0) {
		outCD->AddArray(this->GetFibresVolumes());
	}

	if (this->ComputeTotalFibresVolume != 0) {
		this->GetTotalFibresVolume();	//force its calculation
	}

	if (this->ComputePCSA != 0) {
		this->GetPCSA();	//force its calculation
	}

	return VTK6_RESULT_OK;
}

#if defined(_MSC_VER) && _MSC_VER >= 1600
#pragma region Get the cached data
#endif

//Gets the radius of every fibre vertex, i.e., radius varies along the fibre.
/*virtual*/ vtkDoubleArray* vtkLHPMuscleFibresAnalysisFilter::GetFibresRadii()
{
	if ((this->m_CachedDataFlags & COMPUTE_FIBRES_RADII) == 0)
	{
		if (this->m_FibresRadii == NULL)
		{
			this->m_FibresRadii = vtkDoubleArray::New();
			this->m_FibresRadii->SetName(this->FibresRadiiFieldName);
		}

		this->CalculateFibresRadii(vtkPolyData::SafeDownCast(this->GetInput()), this->GetFibresTangentVectors(), this->m_FibresRadii);

		this->m_CachedDataFlags |= COMPUTE_FIBRES_RADII;
	}

	return this->m_FibresRadii;
}

//Gets the number of segments of every input fibre.
/*virtual*/ vtkIntArray* vtkLHPMuscleFibresAnalysisFilter::GetFibresSegmentsCount()
{
	if ((this->m_CachedDataFlags & COMPUTE_FIBRES_SEGMENTS_COUNT) == 0)
	{
		if (this->m_FibresSegmentsCount == NULL)
		{
			this->m_FibresSegmentsCount = vtkIntArray::New();
			this->m_FibresSegmentsCount->SetName(this->FibresSegmentsCountFieldName);
		}

		vtkLHPMuscleFibresMath::CountFibresSegments(vtkPolyData::SafeDownCast(this->GetInput()), this->m_FibresSegmentsCount);

		this->m_CachedDataFlags |= COMPUTE_FIBRES_SEGMENTS_COUNT;
	}

	return this->m_FibresSegmentsCount;
}

//Gets the lengths of segments of every input fibre.
/*virtual*/ vtkDoubleArray* vtkLHPMuscleFibresAnalysisFilter::GetFibresSegmentsLengths()
{
	if ((this->m_CachedDataFlags & COMPUTE_FIBRES_SEGMENTS_LENGTHS) == 0)
	{
		if (this->m_FibresSegmentsLengths == NULL)
		{
			this->m_FibresSegmentsLengths = vtkDoubleArray::New();
			this->m_FibresSegmentsLengths->SetName(this->FibresSegmentsLengthsFieldName);
		}

		vtkLHPMuscleFibresMath::CalculateFibresSegmentsLengths(vtkPolyData::SafeDownCast(this->GetInput()), this->m_FibresSegmentsLengths);

		this->m_CachedDataFlags |= COMPUTE_FIBRES_SEGMENTS_LENGTHS;
	}

	return this->m_FibresSegmentsLengths;
}

//Gets the lengths of segments of every input fibre.
/*virtual*/ vtkDoubleArray* vtkLHPMuscleFibresAnalysisFilter::GetFibresLengths()
{
	if ((this->m_CachedDataFlags & COMPUTE_FIBRES_LENGTH) == 0)
	{
		if (this->m_FibresLengths == NULL)
		{
			this->m_FibresLengths = vtkDoubleArray::New();
			this->m_FibresLengths->SetName(this->FibresLengthsFieldName);
		}

		vtkLHPMuscleFibresMath::CalculateFibresLengths(this->GetFibresSegmentsCount(),
			this->GetFibresSegmentsLengths(), this->m_FibresLengths);

		this->m_CachedDataFlags |= COMPUTE_FIBRES_LENGTH;
	}

	return this->m_FibresLengths;
}


//Gets the arc-length (Chordal) parameterization of every fibre normalized to 0-1.
/*virtual*/ vtkDoubleArray* vtkLHPMuscleFibresAnalysisFilter::GetFibresChordalParameterization()
{
	if ((this->m_CachedDataFlags & COMPUTE_FIBRES_CHORDAL_PARAMETERIZATION) == 0)
	{
		if (this->m_FibresParameterization == NULL)
		{
			this->m_FibresParameterization = vtkDoubleArray::New();
			this->m_FibresParameterization->SetName(this->FibresChordalParameterizationFieldName);
		}

		vtkLHPMuscleFibresMath::CalculateFibresParameterization(vtkPolyData::SafeDownCast(this->GetInput()),
			this->GetFibresSegmentsLengths(), false, this->m_FibresParameterization);
		vtkLHPMuscleFibresMath::NormalizeFibresParameterization(vtkPolyData::SafeDownCast(this->GetInput()),
			this->m_FibresParameterization, this->m_FibresParameterization);

		this->m_CachedDataFlags |= COMPUTE_FIBRES_CHORDAL_PARAMETERIZATION;
	}

	return this->m_FibresParameterization;
}

//Gets the centripetal parameterization of every fibre normalized to 0-1. 
/*virtual*/ vtkDoubleArray* vtkLHPMuscleFibresAnalysisFilter::GetFibresCentripetalParameterization()
{
	if ((this->m_CachedDataFlags & COMPUTE_FIBRES_CENTRIPETAL_PARAMETERIZATION) == 0)
	{
		if (this->m_FibresCRParameterization == NULL)
		{
			this->m_FibresCRParameterization = vtkDoubleArray::New();
			this->m_FibresCRParameterization->SetName(this->FibresCentripetalParameterizationFieldName);
		}

		vtkLHPMuscleFibresMath::CalculateFibresParameterization(vtkPolyData::SafeDownCast(this->GetInput()),
			this->GetFibresSegmentsLengths(), true, this->m_FibresCRParameterization);
		vtkLHPMuscleFibresMath::NormalizeFibresParameterization(vtkPolyData::SafeDownCast(this->GetInput()),
			this->m_FibresCRParameterization, this->m_FibresCRParameterization);

		this->m_CachedDataFlags |= COMPUTE_FIBRES_CENTRIPETAL_PARAMETERIZATION;
	}

	return this->m_FibresCRParameterization;
}

//Gets the tangent vector of every fibre point.
/*virtual*/ vtkDoubleArray* vtkLHPMuscleFibresAnalysisFilter::GetFibresTangentVectors()
{
	if ((this->m_CachedDataFlags & COMPUTE_FIBRES_TANGENTVECTORS) == 0)
	{
		if (this->m_FibresTangentVectors == NULL)
		{
			this->m_FibresTangentVectors = vtkDoubleArray::New();
			this->m_FibresTangentVectors->SetNumberOfComponents(3);	//vector
			this->m_FibresTangentVectors->SetName(vtkLHPMuscleFibresAnalysisFilter::FibresTangentVectorsFieldName);
		}

		this->CalculateFibresTangentVectors(vtkPolyData::SafeDownCast(this->GetInput()), this->GetFibresParameterization(), this->m_FibresTangentVectors);

		this->m_CachedDataFlags |= COMPUTE_FIBRES_TANGENTVECTORS;
	}

	return this->m_FibresTangentVectors;
}

//Gets the averaged proximal tangent vector of every fibre point.
/*virtual*/ vtkDoubleArray* vtkLHPMuscleFibresAnalysisFilter::GetFibresProximalTangentVectors()
{
	if ((this->m_CachedDataFlags & COMPUTE_FIBRES_PROXIMAL_TANGENTVECTORS) == 0)
	{
		if (this->m_FibresProximalAvgTgVec == NULL)
		{
			this->m_FibresProximalAvgTgVec = vtkDoubleArray::New();
			this->m_FibresProximalAvgTgVec->SetNumberOfComponents(3);	//vector
			this->m_FibresProximalAvgTgVec->SetName(vtkLHPMuscleFibresAnalysisFilter::FibresProximalTangentVectorsFieldName);
		}

		this->CalculateFibresAvgTgVec(vtkPolyData::SafeDownCast(this->GetInput()), this->GetFibresTangentVectors(),
			this->GetFibresChordalParameterization(), true, this->m_FibresProximalAvgTgVec);

		this->m_CachedDataFlags |= COMPUTE_FIBRES_PROXIMAL_TANGENTVECTORS;
	}

	return this->m_FibresProximalAvgTgVec;
}

//Gets the averaged distal tangent vector of every fibre point.
/*virtual*/ vtkDoubleArray* vtkLHPMuscleFibresAnalysisFilter::GetFibresDistalTangentVectors()
{
	if ((this->m_CachedDataFlags & COMPUTE_FIBRES_DISTAL_TANGENTVECTORS) == 0)
	{
		if (this->m_FibresDistalAvgTgVec == NULL)
		{
			this->m_FibresDistalAvgTgVec = vtkDoubleArray::New();
			this->m_FibresDistalAvgTgVec->SetNumberOfComponents(3);	//vector
			this->m_FibresDistalAvgTgVec->SetName(vtkLHPMuscleFibresAnalysisFilter::FibresDistalTangentVectorsFieldName);
		}

		this->CalculateFibresAvgTgVec(vtkPolyData::SafeDownCast(this->GetInput()), this->GetFibresTangentVectors(),
			this->GetFibresChordalParameterization(), false, this->m_FibresDistalAvgTgVec);

		this->m_CachedDataFlags |= COMPUTE_FIBRES_DISTAL_TANGENTVECTORS;
	}

	return this->m_FibresDistalAvgTgVec;
}

//Gets the proximal line of action to which proximal PA are calculated 
/*virtual*/ const double* vtkLHPMuscleFibresAnalysisFilter::GetProximalLineOfAction()
{
	if ((this->m_CachedDataFlags & COMPUTE_PROXIMAL_LINEOFACTION) == 0)
	{
		this->CalculateLineOfAction(vtkPolyData::SafeDownCast(this->GetInput()), this->GetFibresProximalTangentVectors(),
			true, this->ProximalLineOfAction, this->ProximalLineOfActionPos
#ifdef _DEBUG
			, this->ProximalLineOfActionLR
#endif
			);

		this->m_CachedDataFlags |= COMPUTE_PROXIMAL_LINEOFACTION;
	}

	return this->ProximalLineOfAction;
}

//Gets the distal line of action to which distal PA are calculated 
/*virtual*/ const double* vtkLHPMuscleFibresAnalysisFilter::GetDistalLineOfAction()
{
	if ((this->m_CachedDataFlags & COMPUTE_DISTAL_LINEOFACTION) == 0)
	{
		this->CalculateLineOfAction(vtkPolyData::SafeDownCast(this->GetInput()), this->GetFibresDistalTangentVectors(),
			false, this->DistalLineOfAction, this->DistalLineOfActionPos
#ifdef _DEBUG
			, this->DistalLineOfActionLR
#endif
			);

		this->m_CachedDataFlags |= COMPUTE_DISTAL_LINEOFACTION;
	}

	return this->DistalLineOfAction;
}

//Gets the proximal line of action to which proximal PA are calculated 
/*virtual*/ void vtkLHPMuscleFibresAnalysisFilter::GetProximalLineOfAction(double out[3])
{
	const double* app = GetProximalLineOfAction();
	for (int k = 0; k < 3; k++){
		out[k] = app[k];
	}
}

//Gets the distal line of action to which distal PA are calculated 
/*virtual*/ void vtkLHPMuscleFibresAnalysisFilter::GetDistalLineOfAction(double out[3])
{
	const double* app = GetDistalLineOfAction();
	for (int k = 0; k < 3; k++){
		out[k] = app[k];
	}
}

//Gets the position of the proximal line of action., i.e., point through which it comes. 
/*virtual*/ const double* vtkLHPMuscleFibresAnalysisFilter::GetProximalLineOfActionPos()
{
	GetProximalLineOfAction();	//update the values
	return this->ProximalLineOfActionPos;
}

//Gets the position of the distal line of action., i.e., point through which it comes. 
/*virtual*/ const double* vtkLHPMuscleFibresAnalysisFilter::GetDistalLineOfActionPos()
{
	GetDistalLineOfAction();	//update the values
	return this->DistalLineOfActionPos;
}

//Gets the position of the proximal line of action., i.e., point through which it comes. 
/*virtual*/ void vtkLHPMuscleFibresAnalysisFilter::GetProximalLineOfActionPos(double out[3])
{
	const double* app = GetProximalLineOfActionPos();
	for (int k = 0; k < 3; k++){
		out[k] = app[k];
	}
}

//Gets the position of the distal line of action., i.e., point through which it comes. 
/*virtual*/ void vtkLHPMuscleFibresAnalysisFilter::GetDistalLineOfActionPos(double out[3])
{
	const double* app = GetDistalLineOfActionPos();
	for (int k = 0; k < 3; k++){
		out[k] = app[k];
	}
}

#ifdef _DEBUG
//Gets the proximal line of action to which proximal PA are calculated.
/*virtual*/ const double* vtkLHPMuscleFibresAnalysisFilter::GetProximalLineOfActionLR()
{
	GetProximalLineOfAction();	//update the values
	return this->ProximalLineOfActionLR;
}

//Gets the distal line of action to which distal PA are calculated.
/*virtual*/ const double* vtkLHPMuscleFibresAnalysisFilter::GetDistalLineOfActionLR()
{
	GetDistalLineOfAction();	//update the values
	return this->DistalLineOfActionLR;
}

//Gets the proximal line of action to which proximal PA are calculated.
/*virtual*/ void vtkLHPMuscleFibresAnalysisFilter::GetProximalLineOfActionLR(double out[3])
{
	const double* app = GetProximalLineOfActionLR();
	for (int k = 0; k < 3; k++){
		out[k] = app[k];
	}
}

//Gets the distal line of action to which distal PA are calculated.
/*virtual*/ void vtkLHPMuscleFibresAnalysisFilter::GetDistalLineOfActionLR(double out[3])
{
	const double* app = GetDistalLineOfActionLR();
	for (int k = 0; k < 3; k++){
		out[k] = app[k];
	}
}
#endif

//Gets the proximal pennation angle for every fibre.
/*virtual*/ vtkDoubleArray* vtkLHPMuscleFibresAnalysisFilter::GetFibresProximalPA()
{
	if ((this->m_CachedDataFlags & COMPUTE_FIBRES_PROXIMAL_PENNATIONANGLE) == 0)
	{
		if (this->m_FibresProximalPA == NULL)
		{
			this->m_FibresProximalPA = vtkDoubleArray::New();
			this->m_FibresProximalPA->SetName(vtkLHPMuscleFibresAnalysisFilter::FibresProximalPAFieldName);
		}

		this->CalculateFibresPA(this->GetFibresProximalTangentVectors(),
			this->GetProximalLineOfAction(), this->m_FibresProximalPA);

		this->m_CachedDataFlags |= COMPUTE_FIBRES_PROXIMAL_PENNATIONANGLE;
	}

	return this->m_FibresProximalPA;
}

//Gets the distal pennation angle for every fibre. 
/*virtual*/ vtkDoubleArray* vtkLHPMuscleFibresAnalysisFilter::GetFibresDistalPA()
{
	if ((this->m_CachedDataFlags & COMPUTE_FIBRES_DISTAL_PENNATIONANGLE) == 0)
	{
		if (this->m_FibresDistalPA == NULL)
		{
			this->m_FibresDistalPA = vtkDoubleArray::New();
			this->m_FibresDistalPA->SetName(vtkLHPMuscleFibresAnalysisFilter::FibresDistalPAFieldName);
		}

		this->CalculateFibresPA(this->GetFibresDistalTangentVectors(),
			this->GetDistalLineOfAction(), this->m_FibresDistalPA);

		this->m_CachedDataFlags |= COMPUTE_FIBRES_DISTAL_PENNATIONANGLE;
	}

	return this->m_FibresDistalPA;
}

//Gets the volume of every fibre vertex.
/*virtual*/ vtkDoubleArray* vtkLHPMuscleFibresAnalysisFilter::GetFibresVolumes()
{
	if ((this->m_CachedDataFlags & COMPUTE_FIBRES_VOLUME) == 0)
	{
		if (this->m_FibresVolume == NULL)
		{
			this->m_FibresVolume = vtkDoubleArray::New();
			this->m_FibresVolume->SetName(vtkLHPMuscleFibresAnalysisFilter::FibresVolumeFieldName);
		}

		switch (this->FibresVolumeCalculationMethod)
		{
		case VOLUME_CALCULATION_LEE14:
			this->CalculateFibresVolumeLee14(vtkPolyData::SafeDownCast(this->GetInput()),
				this->GetFibresSegmentsLengths(), this->m_FibresVolume);
			break;

		case VOLUME_CALCULATION_CYLINDRICAL:
			this->CalculateFibresVolumeCylindrical(vtkPolyData::SafeDownCast(this->GetInput()),
				this->GetFibresSegmentsLengths(), this->GetFibresRadii(), this->m_FibresVolume);
			break;

		case VOLUME_CALCULATION_RAVICHANDIRAN:
			this->CalculateFibresVolumeRavichandrian(vtkPolyData::SafeDownCast(this->GetInput()), this->GetFibresLengths(),
				this->GetFibresRadii(), this->m_FibresVolume);
			break;

		default:
			assert(false);
		}

		this->m_CachedDataFlags |= COMPUTE_FIBRES_VOLUME;
	}

	return this->m_FibresVolume;
}

//Gets the total volume of all fibres.
/*virtual*/ double vtkLHPMuscleFibresAnalysisFilter::GetTotalFibresVolume()
{
	if ((this->m_CachedDataFlags & COMPUTE_TOTAL_FIBRES_VOLUME) == 0)
	{
		this->m_TotalVolume = vtkLHPMuscleFibresMath::Sum(this->GetFibresVolumes());

		this->m_CachedDataFlags |= COMPUTE_TOTAL_FIBRES_VOLUME;
	}

	return this->m_TotalVolume;
}

//Gets the physiological cross-section area.
/*virtual*/ double vtkLHPMuscleFibresAnalysisFilter::GetPCSA()
{
	if ((this->m_CachedDataFlags & COMPUTE_PCSA) == 0)
	{
		//according Ravichandiran et al. 2009, PCSA = (muscle volume * cos (average PA) / average nFBL
		//where average nFBL = averaged normalized fibre length
		//and average PA = sum(PAdistal + PAproximal) / 2*nFibres


		double averagePA = (vtkLHPMuscleFibresMath::Sum(this->GetFibresProximalPA())
			+ vtkLHPMuscleFibresMath::Sum(this->GetFibresDistalPA()))
			/ (2 * vtkPolyData::SafeDownCast(this->GetInput())->GetNumberOfCells());

		this->m_PCSA = this->GetTotalFibresVolume() * cos(averagePA) / vtkLHPMuscleFibresMath::Avg(this->GetFibresLengths());

		this->m_CachedDataFlags |= COMPUTE_PCSA;
	}

	return this->m_PCSA;
}

#if defined(_MSC_VER) && _MSC_VER >= 1600
#pragma endregion Get the cached data
#endif

//Computes the volume of every fibre using the method described by Lee et al. 2014.
void  vtkLHPMuscleFibresAnalysisFilter::CalculateFibresVolumeLee14(const vtkPolyData* poly,
	const vtkDoubleArray* segLengths, vtkDoubleArray* outVolume)
{
	assert(false);

	//TODO: Not implemented
}

//Computes the volume of every fibre using the own proprietal method.
void  vtkLHPMuscleFibresAnalysisFilter::CalculateFibresVolumeCylindrical(const vtkPolyData* poly,
	const vtkDoubleArray* segLengths, const vtkDoubleArray* fibRadii, vtkDoubleArray* outVolume)
{
	//The volume of a fibre is given as sum of truncated cones having the height of segment length and 
	//radii as measured in the vertices of the fibre

	int nFibres = const_cast<vtkPolyData*>(poly)->GetNumberOfCells();
	const double* pR = const_cast<vtkDoubleArray*>(fibRadii)->GetPointer(0);
	const double* pLen = const_cast<vtkDoubleArray*>(segLengths)->GetPointer(0);
	double* pVol = outVolume->WritePointer(0, nFibres);

	double pithird = vtkMath::Pi() / 3;

	for (int i = 0; i < nFibres; i++) //cannot be easily parallelized
	{
		vtkIdType nPts, *pPts;
		const_cast<vtkPolyData*>(poly)->GetCellPoints(i, nPts, pPts);

		//calculate the volume
		//V of the truncated cone is 1/3 * PI * (r1^2 + r2^2 + r1*r2)* h
		double totVol = 0.0;
		for (int j = 1; j < nPts; j++) {
			totVol += (
				pR[pPts[j - 1]] * pR[pPts[j - 1]] +
				pR[pPts[j - 1]] * pR[pPts[j]] +
				pR[pPts[j]] * pR[pPts[j]]
				)*pLen[j - 1];
		}

		//calculate volume
		pVol[i] = totVol * pithird;

		pLen += (nPts - 1); //move to next Fibre
	}
}


//Computes the volume of every fibre using the method described by Ravichandrian et al. 2009.
void  vtkLHPMuscleFibresAnalysisFilter::CalculateFibresVolumeRavichandrian(const vtkPolyData* poly,
	const vtkDoubleArray* fibLengths, const vtkDoubleArray* fibRadii, vtkDoubleArray* outVolume)
{
	//Ravichandrian et al. 2009: To calculate muscle volume, the ﬁbre bundles are treated as 3D cylinders.
	//The fibre length FBL was used as the height of each cylinder. The radius of the cylinder(R) is calculated 
	//for each ﬁbre bundle as half the mean distance between its b-spline curve and the b-spline curve of 
	//the nearest neighbouring fibre bundle.

	int nFibres = const_cast<vtkPolyData*>(poly)->GetNumberOfCells();
	const double* pR = const_cast<vtkDoubleArray*>(fibRadii)->GetPointer(0);
	const double* pLen = const_cast<vtkDoubleArray*>(fibLengths)->GetPointer(0);
	double* pVol = outVolume->WritePointer(0, nFibres);

#pragma omp parallel for shared(pVol)
	for (int i = 0; i < nFibres; i++)
	{
		vtkIdType nPts, *pPts;
		const_cast<vtkPolyData*>(poly)->GetCellPoints(i, nPts, pPts);

		//calculate mean R
		double meanR = 0.0;
		for (int j = 0; j < nPts; j++) {
			meanR += pR[pPts[j]];
		}

		meanR /= nPts;

		//calculate volume
		pVol[i] = pLen[i] * vtkMath::Pi()*meanR*meanR;
	}
}



//Computes the tangent vector at every vertex of every fibre and stores them into outTgv.
void vtkLHPMuscleFibresAnalysisFilter::CalculateFibresTangentVectors(
	const vtkPolyData* poly, const vtkDoubleArray* params, vtkDoubleArray* outTgv)
{
	int nFibres = const_cast<vtkPolyData*>(poly)->GetNumberOfCells();
	const double* pPars = const_cast<vtkDoubleArray*>(params)->GetPointer(0);

	VCoord* pTg = (VCoord*)this->m_FibresTangentVectors->WritePointer(0,
		3 * const_cast<vtkPolyData*>(poly)->GetNumberOfPoints());

#pragma omp parallel for shared(pTg)
	for (int i = 0; i < nFibres; i++)
	{
		vtkIdType nPts, *pPts;
		const_cast<vtkPolyData*>(poly)->GetCellPoints(i, nPts, pPts);

		//Let the fibre be formed by points Q0, Q1, ... Qn-1 having Catmull-Rom parameterization Q0_t, Q1_t, ...Qn-1_t
		//Cubic Catmul-Rom curve C is defined by 4 points P0..P3 and their knots (parameterization) t0..t3
		//where knots ti+1 are defined as ti+1 = ti + ||Pi+1 - Pi||^alpha 
		//whereby alpha = 1/2 for centripetal CR, 1 for chordal
		//see: http://en.wikipedia.org/wiki/Centripetal_Catmull%E2%80%93Rom_spline
		//Note: the curve is valid only from P1 to P2 (points P0 and P3 are used to define it), i.e.,
		//the valid range of the parameter t is <t1...t2>
		//
		//
		//C'(t1) = (P1 - P0) / (t1 - t0) - (P2 - P0) / (t2 - t0) + (P2 - P1) / (t2 - t1)
		//C'(t2) = (P2 - P1) / (t2 - t1) - (P3 - P1) / (t3 - t1) + (P3 - P2) / (t3 - t2)
		//
		//	Both formulas personally verified in Matlab
		//
		//For Q0 the tangent vector must be approximated
		//For Q1, it is calculated using C'(t1) formula for Q0, Q1, Q2, Q3 being P0, .., P3
		//For Q2, it is calculated using C'(t1) formula for Q1, Q2, Q3, Q4 being P0, .., P3
		// ...
		//For Qn-2, it is calculated using C'(t2) formula for Qn-4, Qn-3, Qn-2, Qn-1 being P0, .., P3
		//For Qn-1 the tangent vector must be approximated


		double P[4][3];	//buffer for the points P0..P3
		double t[4];		//buffer for the knots values t0..t3

		const_cast<vtkPolyData*>(poly)->GetPoint(pPts[0], P[1]); t[1] = pPars[pPts[0]];	//P1 is the first point of the fibre
		const_cast<vtkPolyData*>(poly)->GetPoint(pPts[1], P[2]); t[2] = pPars[pPts[1]];	//P2 is the second point of the fibre

		//construct an artificial point P0 to define the Catmull-Rom cubic spline in the
		//in the opposite direction from Q1 as Q0A = 2*Q0 - Q1
		for (int k = 0; k < 3; k++) {
			P[0][k] = 2 * P[1][k] - P[2][k];
		}

		t[0] = t[1] - t[2];

		for (int j = 0; j < nPts - 1; j++)
		{
			//load the fourth point			
			int iv0 = (j + 0) % 4;
			int iv1 = (j + 1) % 4;
			int iv2 = (j + 2) % 4;
			int iv3 = (j + 3) % 4;

			if (j + 2 < nPts)
			{
				const_cast<vtkPolyData*>(poly)->GetPoint(pPts[j + 2], P[iv3]);
				t[iv3] = pPars[pPts[j + 2]];
			}
			else
			{
				//construct an artificial point P0 to define the Catmull-Rom cubic spline in the
				//in the opposite direction from Q2 as Q3A = 2*Q2 + Q1
				for (int k = 0; k < 3; k++) {
					P[iv3][k] = 2 * P[iv2][k] - P[iv1][k];
				}

				t[iv3] = 2 * t[iv2] - t[iv1];

				//and calculate the tangent vector at t2, i.e., at the last vertex of the fibre:
				//C'(t2) = (P2 - P1) / (t2 - t1) - (P3 - P1) / (t3 - t1) + (P3 - P2) / (t3 - t2)
				int index = pPts[j + 1];
				for (int k = 0; k < 3; k++) {
					pTg[index][k] = (P[iv2][k] - P[iv1][k]) / (t[iv2] - t[iv1]) -
						(P[iv3][k] - P[iv1][k]) / (t[iv3] - t[iv1]) +
						(P[iv3][k] - P[iv2][k]) / (t[iv3] - t[iv2]);
				}
			}

			//calculate the tangent vector at t0, i.e., the current point
			//C'(t1) = (P1 - P0) / (t1 - t0) - (P2 - P0) / (t2 - t0) + (P2 - P1) / (t2 - t1)
			int index = pPts[j];
			for (int k = 0; k < 3; k++) {
				pTg[index][k] = (P[iv1][k] - P[iv0][k]) / (t[iv1] - t[iv0]) -
					(P[iv2][k] - P[iv0][k]) / (t[iv2] - t[iv0]) +
					(P[iv2][k] - P[iv1][k]) / (t[iv2] - t[iv1]);
			}
		} //end for
	}

	//done.
}

//Computes the averaged tangent vectors for every fibre either for its proximal or distal part.
void vtkLHPMuscleFibresAnalysisFilter::CalculateFibresAvgTgVec(const vtkPolyData* poly,
	const vtkDoubleArray* tgvecs, const vtkDoubleArray* chordalParams, bool proximal, vtkDoubleArray* outAvgTgv)
{
	//Calculate the averaged tangent vectors for proximal and distal orientations for every fibre 
	//whereby a local proximal or distal region  of a fibre is considered to be 15-20% (this->ProximalThreshold) 
	//of the entire fibre length	- see Lee et al., 2014

	//common initialization	
	int nFibres = const_cast<vtkPolyData*>(poly)->GetNumberOfCells();

	double thr[2];		//threshold in par
	int rng[3];				//ranges for the internal for-cycle to speed-up the processing	

	if (proximal)
	{
		thr[0] = this->ProximalThreshold;
		thr[1] = 0.0;

		rng[1] = +1;
	}
	else
	{
		thr[0] = 1.0;
		thr[1] = 1.0 - this->DistalThreshold;

		rng[1] = -1;
	}

	rng[1 - rng[1]] = 0;


	VCoord* pAvgTgs = (VCoord*)outAvgTgv->WritePointer(0, 3 * nFibres);
	memset(pAvgTgs, 0, sizeof(VCoord)*nFibres);

	const double* pPars = const_cast<vtkDoubleArray*>(chordalParams)->GetPointer(0);
	const VCoord* pTgs = (VCoord*)const_cast<vtkDoubleArray*>(tgvecs)->GetPointer(0);
#pragma omp parallel for shared(pAvgTgs) firstprivate(rng)
	for (int i = 0; i < nFibres; i++)
	{
		vtkIdType nPts, *pPts;
		const_cast<vtkPolyData*>(poly)->GetCellPoints(i, nPts, pPts);

		rng[1 + rng[1]] = nPts - 1;

		//get the region
		int index = rng[0];
		while (index != rng[2])
		{
			double t = pPars[pPts[index]];
			if (t > thr[0] || t < thr[1])
				break;	//we are outside the region

			//get the tangent vector
			const VCoord& tgvec = pTgs[pPts[index]];
			for (int k = 0; k < 3; k++){
				pAvgTgs[i][k] += tgvec[k];
			}

			index += rng[1];
		} //end while

		index = abs(index - rng[0]);
		for (int k = 0; k < 3; k++){
			pAvgTgs[i][k] /= index;
		}
	}
}

//------------------------------------------------------------------------
//Computes centroid of points
void vtkLHPMuscleFibresAnalysisFilter::ComputeCentroid(const vtkPoints* points, VCoord centroid)
//------------------------------------------------------------------------
{
	int N = const_cast<vtkPoints*>(points)->GetNumberOfPoints();
	centroid[0] = centroid[1] = centroid[2] = 0.0;
	for (int i = 0; i < N; i++)
	{
		const double* pcoords = const_cast<vtkPoints*>(points)->GetPoint(i);
		for (int j = 0; j < 3; j++){
			centroid[j] += pcoords[j];
		}
	}

	for (int j = 0; j < 3; j++){
		centroid[j] /= N;
	}
}

//------------------------------------------------------------------------
// Computes three eigen vectors of the given point set.
//Centroid may be NULL, if it should be calculated automatically, otherwise, it must be centroid of points.
//The returned eigen vectors are ordered by their lengths. The longest (i.e., the principal axis) is denoted by the first one.
void vtkLHPMuscleFibresAnalysisFilter::ComputeEigenVects(const vtkPoints* points, const VCoord centroid, VCoord eigenvects[3])
//------------------------------------------------------------------------
{
	double origin[3];
	if (centroid == NULL) {
		ComputeCentroid(points, origin);
		centroid = origin;
	}

	//compute the covariance matrix, which is a symmetric matrix
	double A[3][3], eigenvals[3], eigenvectsT[3][3];

	//fill the matrix
	memset(A, 0, sizeof(A));
	int N = const_cast<vtkPoints*>(points)->GetNumberOfPoints();
	for (int k = 0; k < N; k++)
	{
		const double* pcoords = const_cast<vtkPoints*>(points)->GetPoint(k);

		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++) {
				A[i][j] += (pcoords[i] - centroid[i])*(pcoords[j] - centroid[j]);
			}
		}
	}

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++){
			A[i][j] /= (N - 1);
		}
	}

	//compute eigen vectors, the principal axis is the first one
	double *ATemp[3], *V[3];
	for (int i = 0; i < 3; i++)
	{
		ATemp[i] = A[i];
		V[i] = eigenvectsT[i];
	}

	vtkMath::Jacobi(ATemp, eigenvals, V);

	//copy the result
	//N.B. Jacobi returns vectors in columns!
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			eigenvects[i][j] = eigenvectsT[j][i];
		}
	}
}

//Computes the line of action for proximal or distal part.
void vtkLHPMuscleFibresAnalysisFilter::CalculateLineOfAction(const vtkPolyData* poly,
	const vtkDoubleArray* avgtgvecs, bool proximal, double la[3], double laPos[3]
#ifdef _DEBUG
	, double laLR[3]
#endif
	)
{
	//the Line of Action is calculated as described in Lee et al., 2014
	//http://dx.doi.org/10.1080/10255842.2014.917294

	//It requires to 1) calculate the principal axes from the attachment points of fibres (end-points) 
	//to automatically determine if the muscle is pennate or not
	//2) if the muscle is pennate, the lines of action are the principal axes
	//3) if the muscle is not pennate, calculate the averaged tangent vectors for proximal 
	// and distal orientations for every fibre and then the averaged vectors for the entire muscle
	//The lines of action are these averaged vectors.


	//get the principal axes of the end-point	
	int nFibres = const_cast<vtkPolyData*>(poly)->GetNumberOfCells();

	vtkPoints* points = vtkPoints::New();
	for (int i = 0; i < nFibres; i++)
	{
		vtkIdType nPts, *pPts;
		const_cast<vtkPolyData*>(poly)->GetCellPoints(i, nPts, pPts);

		points->InsertNextPoint(const_cast<vtkPolyData*>(poly)->GetPoint(pPts[proximal ? 0 : nPts - 1]));
	}

	VCoord eigenvects[3];
	ComputeCentroid(points, laPos);
	ComputeEigenVects(points, laPos, eigenvects);

	//the principal axis is defined by the centroid and the first eigen vectors, i.e., eigenvects[0]
	vtkMath::Normalize(eigenvects[0]);	//make sure it is unit vector

#ifdef _DEBUG	
	for (int k = 0; k < 3; k++) {
		laLR[k] = eigenvects[0][k];
	}
#endif

	//evaluate the quality of the fit (r2)
	double nom = 0.0, denom = 0.0;
	for (int i = 0; i < nFibres; i++)
	{
		double* pt = points->GetPoint(i);

		//calculate the distance between pt (Q) and the principal axes
		//L(t) = P + d*t, where P is centroid and d is unit vector (eigenvects[0])		
		VCoord difPQ, ptProj;
		for (int k = 0; k < 3; k++)
		{
			difPQ[k] = pt[k] - laPos[k];
			denom += difPQ[k] * difPQ[k];
		}

		double ti = vtkMath::Dot(eigenvects[0], difPQ);
		for (int k = 0; k < 3; k++) {
			ptProj[k] = laPos[k] + eigenvects[0][k] * ti;
		}

		nom += vtkMath::Distance2BetweenPoints(pt, ptProj);
	}

	points->Delete();

	double r2 = 1.0 - nom / denom;

	//according to experiments by Lee et al., r2 = 0.92 was the minimal value found for a pennate muscle, 
	//but there was also 0.96 non-pennate muscle => it is difficult to guess the best threshold, Lee used 0.9
	if (r2 > 0.95) //5% reliability
	{
		for (int k = 0; k < 3; k++) {
			la[k] = eigenvects[0][k];
		}

		return;
	}

	//if the muscle is not pennate, calculate the averaged tangent vectors for proximal
	//and distal orientations for every fibre whereby a local proximal or distal region 
	//of a fibre is considered to be 15-20% (this->ProximalThreshold) of the entire fibre length
	const VCoord* pAvgTgs = (VCoord*)const_cast<vtkDoubleArray*>(avgtgvecs)->GetPointer(0);
	memset(la, 0, sizeof(VCoord));

	for (int i = 0; i < nFibres; i++)
	{
		for (int k = 0; k < 3; k++){
			la[k] += pAvgTgs[i][k];
		}
	}

	//and now compute the averaged vectors for the entire muscle
	for (int k = 0; k < 3; k++){
		la[k] /= nFibres;
	}

	vtkMath::Normalize(la);	//make it unit vector
}

//Computes the pennation angle for every fibre.
void vtkLHPMuscleFibresAnalysisFilter::CalculateFibresPA(const vtkDoubleArray* avgtgvecs,
	const double lineOfAction[3], vtkDoubleArray* outPA
	)
{
	//common initialization
	int nFibres = const_cast<vtkDoubleArray*>(avgtgvecs)->GetNumberOfTuples();
	const VCoord* pAvgTgVecs = (VCoord*)const_cast<vtkDoubleArray*>(avgtgvecs)->GetPointer(0);
	double* pPAF = outPA->WritePointer(0, nFibres);

	double normLA = 1.0 / vtkMath::Norm(lineOfAction);
#pragma omp parallel for shared(pPAF)
	for (int i = 0; i < nFibres; i++){
		pPAF[i] = acos(normLA * vtkMath::Dot(lineOfAction, pAvgTgVecs[i]) / vtkMath::Norm(pAvgTgVecs[i]));
	}
}

//Computes the radius of a fibre at each of its vertex.
void vtkLHPMuscleFibresAnalysisFilter::CalculateFibresRadii(const vtkPolyData* poly,
	const vtkDoubleArray* tgvecs, vtkDoubleArray* outRadii)
{
	int nFibres = const_cast<vtkPolyData*>(poly)->GetNumberOfCells();
	const VCoord* pTgs = (VCoord*)const_cast<vtkDoubleArray*>(tgvecs)->GetPointer(0);

	double* pR = outRadii->WritePointer(0, const_cast<vtkDoubleArray*>(tgvecs)->GetNumberOfTuples());

	//brute-force approach: for each fibre vertex, construct a plane passing through this vertex and being
	//perpendicular to the fibre and then use this plane to cut every other fibre and indentify the point of 
	//intersection that is nearest to the fibre vertex being processed - the distance between these two points 
	//is the radius we are searching for
#pragma omp parallel for shared(pR)
	for (int i = 0; i < nFibres; i++)
	{
		vtkIdType nPts, *pPts;
		const_cast<vtkPolyData*>(poly)->GetCellPoints(i, nPts, pPts);

		double lastR = 0.0;
		for (int j = 0; j < nPts; j++)
		{
			double pt[3];
			const_cast<vtkPolyData*>(poly)->GetPoint(pPts[j], pt);

			//the plane is defined by the point pt and normal *pTgs		
			pR[pPts[j]] = GetNearestFibreDistance(poly, pt, pTgs[pPts[j]], i) / 2;

			if (pR[pPts[j]] == 0.0)
			{
				//there is no other fibre cut by the given plane => use the same radius as was the previous one
				pR[pPts[j]] = lastR;
			}
			else
			{
				if (lastR == 0.0)
				{
					//vertices 0 .. j-1 have currently zero radius => replace it by *pR
					//pR points to 
					for (int k = 0; k < j; k++) {
						pR[pPts[k]] = pR[pPts[j]];
					}
				}

				lastR = pR[pPts[j]];
			}
		}

		if (lastR == 0.0)
		{
			//singular case, as a last resort specify some constant
			for (int j = 0; j < nPts; j++) {
				pR[pPts[j]] = 0.1;	//some constant
			}
		}
	}
}

//Gets the distance between the given point and its nearest fibre.
double vtkLHPMuscleFibresAnalysisFilter::GetNearestFibreDistance(const vtkPolyData* poly,
	const VCoord pt, const VCoord normal, int excludeFibreId)
{
	//The calculation is done on a plane defined by the given point(pt) and normal(normal) and
	//	the fibre with excludeFibreId identifier is excluded from the calculation.
	//	The method returns 0.0, if there is no fibre intersected by the plane.
	//	excludeFibreId may by - 1, if every fibres should be tested.

	double dblRet = 0.0;

	int nFibres = const_cast<vtkPolyData*>(poly)->GetNumberOfCells();
	for (int i = 0; i < nFibres; i++)
	{
		if (i == excludeFibreId)
			continue;	//skip this fibre

		double dist = Distance2BetweenFibreAndPoint(poly, i, pt);
		if (dblRet == 0.0 || dist < dblRet) {
			dblRet = dist;
		}
	}

	return sqrt(dblRet);
}

//Gets the square distance between the given point and the fibre.
double vtkLHPMuscleFibresAnalysisFilter::Distance2BetweenFibreAndPoint(const vtkPolyData* poly,
	int fibreId, const VCoord pt)
{
	vtkIdType nPts, *pPts;
	const_cast<vtkPolyData*>(poly)->GetCellPoints(fibreId, nPts, pPts);

	double AB[2][3];	//buffer containing both end-points of the segment
	const_cast<vtkPolyData*>(poly)->GetPoint(pPts[0], AB[0]);

	//distance to the first point of the curve
	double dblRet = vtkMath::Distance2BetweenPoints(pt, AB[0]);

	int iFirstPt = 0;
	for (int j = 1; j < nPts; j++)
	{
		const_cast<vtkPolyData*>(poly)->GetPoint(pPts[j], AB[1 - iFirstPt]);	//get the next point		

		//calculate the distance of the point pt to line segment AB
		//Line: L(t) = A + (B-A)*t, Q = projection of P onto the line, i.e.,
		//dist(QA) = (P-A)*(B-A) / |B - A| => Q = A + (B-A)*distQA =>
		//Q = A + (B-A)*(B-A)*(P-A) / |B - A|

		VCoord difPA, difBA, ptProj;
		for (int k = 0; k < 3; k++)
		{
			difPA[k] = pt[k] - AB[0][k];
			difBA[k] = AB[1][k] - AB[0][k];
		}

		double t = vtkMath::Dot(difPA, difBA);
		if (t < 0.0)
			continue; //no intersection here

		t /= vtkMath::Dot(difBA, difBA);
		if (t > 1.0)
			continue; //no intersection here

		for (int k = 0; k < 3; k++) {
			ptProj[k] = AB[0][k] + difBA[k] * t;
		}

		double dist = vtkMath::Distance2BetweenPoints(pt, ptProj);
		if (dist < dblRet) {
			dblRet = dist;
		}

		iFirstPt = 1 - iFirstPt;
	}

	//and finaly test the distance to the last point
	double dist = vtkMath::Distance2BetweenPoints(pt, AB[1 - iFirstPt]);
	if (dist < dblRet) {
		dblRet = dist;
	}

	return dblRet;
}

