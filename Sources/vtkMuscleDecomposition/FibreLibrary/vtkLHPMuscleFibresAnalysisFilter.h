#pragma once

#ifndef __vtkLHPMuscleFibresAnalysisFilter_h
#define __vtkLHPMuscleFibresAnalysisFilter_h

#include "vtkObject.h"
#if VTK_MAJOR_VERSION < 5
#include "vtkPolyDataToPolyDataFilter.h"
#define vtkPolyDataAlgorithm vtkPolyDataToPolyDataFilter
#else
#include "vtkPolyDataAlgorithm.h"
#endif

class vtkDoubleArray;
class vtkIntArray;
class vtkPoints;

/**
 Muscle fibres analysis filter
 
 Performs analysis of the muscle fibres passed into input. The requirements are:
 1) Fibres are represented by poly-lines,  2) One poly-line is stored as one cell,
 i.e., the cell data structure of the polydata object is:

 M_0, i_{0,0}, i_{0,1}, ... i_{0,M0-1},
 M_1, i_{1,0}, i_{1,1}, ... i_{1,M1-1},
 ...
 M_{N-1}, i_{N-1,0}, i_{N-1,1}, ... i_{N-1,M1-1},

 where N is the number fibres, M_f is the number of vertices forming the poly-line representing the f-th fibre,
 i_{f, v} is the index of v-th vertex of f-th fibre.

 3) No vertex can be used repeatedly, and 4) every vertex belongs to some fibre.
 
 The output polydata contains the same data as the input plus in addition there are several fields:

 CELL DATA FIELDS:
 Fibres Lengths					- SCALAR field -	N x double value
 Fibres Proximal Pennate Angle (PA)	- SCALAR field -	N x double value
 Fibres Distal Pennate Angle (PA)		- SCALAR field -	N x double value
 Fibres Volumes					- SCALAR field -	N x double value

 POINT DATA FIELDS:
 Fibres Parameterization		- SCALAR field -	M_0 x double value, M_1 x double value, ...
 Fibres Tangent Vector			- VECTOR field -	M_0 x xyz double values, M_1 x xyz double values, ...
 
 The construction of above given fields depend on the settings of the filter.
 The filter also can provide additional output information (if required):

 TotalFibresVolume	- the sum of fibres volumes
 Physiological Cross-Section Area (PCSA) -	calculated as Kajeandra Ravichandiran, Mayoorendra Ravichandiran, Michele L. Oliver, 
									Karan S. Singh, Nancy H. McKee, Anne M.R. Agur, Determining physiological cross-sectional 
									area of extensor carpi radialis longus and brevis as a whole and by regions using 3D 
									computer muscle models created from digitized fiber bundle data, 
									Computer Methods and Programs in Biomedicine, Volume 95, Issue 3, September 2009, 
									Pages 203-212, ISSN 0169-2607, http://dx.doi.org/10.1016/j.cmpb.2009.03.002.

 Proximal and Distal Lines of Action - if not given as input, these are calculated using the automatic technique described in:
									Lee, Dongwoon and Li, Zhi and Sohail, Qazi Zain and Jackson, Ken and Fiume, Eugene and Agur, Anne,
									A three-dimensional approach to pennation angle estimation for human skeletal muscle,
									Computer Methods in Biomechanics and Biomedical Engineering, 2014,
									http://dx.doi.org/10.1080/10255842.2014.917294
 */
class vtkLHPMuscleFibresAnalysisFilter : public vtkPolyDataAlgorithm
{
public:
	static const char* FibresLengthsFieldName;
	static const char* FibresProximalPAFieldName;
	static const char* FibresDistalPAFieldName;

	static const char* FibresChordalParameterizationFieldName;
	static const char* FibresCentripetalParameterizationFieldName;
	static const char* FibresTangentVectorsFieldName;

	static const char* FibresRadiiFieldName;
	static const char* FibresVolumeFieldName;

protected:
	static const char* FibresSegmentsCountFieldName;
	static const char* FibresSegmentsLengthsFieldName;
	static const char* FibresProximalTangentVectorsFieldName;
	static const char* FibresDistalTangentVectorsFieldName;

protected:
	enum COMPUTE_FLAGS
	{
		COMPUTE_NOTHING = 0x0,						//no change at all, just skip the processing

		COMPUTE_FIBRES_LENGTH = 0x1,												//fibres lengths must be calculated
		COMPUTE_FIBRES_CHORDAL_PARAMETERIZATION = 0x2,			//fibres chordal parameterization must be calculated		
		COMPUTE_FIBRES_CENTRIPETAL_PARAMETERIZATION = 0x4,	//fibres centripetal parameterization must be calculated
		COMPUTE_FIBRES_TANGENTVECTORS = 0x8,								//fibres tangent vector must be calculated
		COMPUTE_FIBRES_PROXIMAL_TANGENTVECTORS = 0x10,			//the averaged tangent vectors of the proximal region for each fibre
		COMPUTE_FIBRES_DISTAL_TANGENTVECTORS = 0x20,				//the averaged tangent vectors of the distal region for each fibre 
		COMPUTE_PROXIMAL_LINEOFACTION = 0x40,								//the proximal line of action must be calculated
		COMPUTE_DISTAL_LINEOFACTION = 0x80,									//the distal line of action must be calculated
		COMPUTE_FIBRES_PROXIMAL_PENNATIONANGLE = 0x100,			//fibres proximal pennate angle must be calculated
		COMPUTE_FIBRES_DISTAL_PENNATIONANGLE = 0x200,				//fibres distal pennate angle must be calculated		
		COMPUTE_FIBRES_SEGMENTS_LENGTHS = 0x400,						//fibres segments lengths must be calculated
		COMPUTE_FIBRES_SEGMENTS_COUNT = 0x800,							//fibres number of segments must be calculated
		COMPUTE_FIBRES_RADII = 0x1000,											//sizes (radius) of every fibre at every vertex must be calculated		
		COMPUTE_FIBRES_VOLUME = 0x2000,											//fibres volumes must be calculated
		COMPUTE_TOTAL_FIBRES_VOLUME = 0x4000,								//fibres volumes must be calculated
		COMPUTE_PCSA = 0x8000,															//PCSA must be calculated	

		COMPUTE_EVERYTHING = -1 ,					//typically during serialization, or if the input has changed
	};

	enum SETTINGS_CHANGED_FLAGS
	{
		NOTHINGCHANGED = 0x0,
		USECENTRIPETALCATMULLROM_CHANGED  = 0x1,
		PROXIMALTHRESHOLD_CHANGED = 0x2,
		DISTALTHRESHOLD_CHANGED = 0x4,
		VOLUME_CALCULATION_METHOD_CHANGED = 0x8,
	};

public:
	enum VOLUME_CALCULATION_METHODS
	{
		VOLUME_CALCULATION_RAVICHANDIRAN,		//Length*PI*avgR^2
		VOLUME_CALCULATION_LEE14,						//area of Voronoi
		VOLUME_CALCULATION_CYLINDRICAL,			//My own
	};

#ifndef vtkSetMacroEx
#define vtkSetMacroEx(name,type, modifiedPar) \
virtual void Set##name (type _arg) \
			{ \
	vtkDebugMacro(<< this->GetClassName() << " (" << this << "): setting " #name " to " << _arg); \
	if (this->name != _arg) \
			{ \
		this->name = _arg; \
		this->Modified(modifiedPar); \
			} \
			} 
#endif

	//define vector
	typedef double VCoord[3];

public:
	vtkTypeMacro(vtkLHPMuscleFibresAnalysisFilter, vtkPolyDataAlgorithm);

	/** Constructs a new analysis filter that calculates nothing by the default. */
	static vtkLHPMuscleFibresAnalysisFilter *New();

	/** Turn on/off the computation of fibres length */
	vtkBooleanMacro(ComputeFibresLengths, int);
	vtkSetMacro(ComputeFibresLengths, int);
	vtkGetMacro(ComputeFibresLengths, int);

	/** Turn on/off the computation of the parameterization of fibres.
	It will be either chordal or centripetal depending on the UseCentripetalCatmullRom option. */
	vtkBooleanMacro(ComputeFibresParameterization, int);
	vtkSetMacro(ComputeFibresParameterization, int);
	vtkGetMacro(ComputeFibresParameterization, int);

	/** Turn on/off the using of centripal Catmull-Rom parameterization of fibres
	If off, the default chordal Catmull-Rom parameterization is used. */
	vtkBooleanMacro(UseCentripetalCatmullRom, int);
	vtkSetMacroEx(UseCentripetalCatmullRom, int, USECENTRIPETALCATMULLROM_CHANGED);
	vtkGetMacro(UseCentripetalCatmullRom, int);
	
	/** Turn on/off the computation of  tangent vectors for each node of every fibre */
	vtkBooleanMacro(ComputeFibresTangentVector, int);
	vtkSetMacro(ComputeFibresTangentVector, int);
	vtkGetMacro(ComputeFibresTangentVector, int);

	/** Turn on/off the computation of the proximal line of action */
	vtkBooleanMacro(ComputeProximalLineOfAction, int);
	vtkSetMacro(ComputeProximalLineOfAction, int);
	vtkGetMacro(ComputeProximalLineOfAction, int);

	/** Turn on/off the computation of the distal line of action */
	vtkBooleanMacro(ComputeDistalLineOfAction, int);
	vtkSetMacro(ComputeDistalLineOfAction, int);
	vtkGetMacro(ComputeDistalLineOfAction, int);

	/** Sets of gets the fraction  of the fibre used for calculating the proximal line of action (0..1) */
	vtkSetMacroEx(ProximalThreshold, double, PROXIMALTHRESHOLD_CHANGED);
	vtkGetMacro(ProximalThreshold, double);

	/** Sets of gets the fraction  of the fibre used for calculating the distal line of action (0..1) */
	vtkSetMacroEx(DistalThreshold, double, DISTALTHRESHOLD_CHANGED);
	vtkGetMacro(DistalThreshold, double);

	/** Turn on/off the computation of the volumes of fibres */
	vtkBooleanMacro(ComputeFibresVolumes, int);
	vtkSetMacro(ComputeFibresVolumes, int);
	vtkGetMacro(ComputeFibresVolumes, int);

	/** Turn on/off the computation of the total volume of the muscle via summation of fibres volumes */
	vtkBooleanMacro(ComputeTotalFibresVolume, int);
	vtkSetMacro(ComputeTotalFibresVolume, int);
	vtkGetMacro(ComputeTotalFibresVolume, int);

	/** Turn on/off the computation of the proximal PA */
	vtkBooleanMacro(ComputeProximalPennationAngle, int);	
	vtkSetMacro(ComputeProximalPennationAngle, int);
	vtkGetMacro(ComputeProximalPennationAngle, int);

	/** Turn on/off the computation of the distal PA */
	vtkBooleanMacro(ComputeDistalPennationAngle, int);
	vtkSetMacro(ComputeDistalPennationAngle, int);
	vtkGetMacro(ComputeDistalPennationAngle, int);

	/** Turn on/off the computation of fibres radii */
	vtkBooleanMacro(ComputeFibresRadii, int);
	vtkSetMacro(ComputeFibresRadii, int);
	vtkGetMacro(ComputeFibresRadii, int);

	/** Gets/Sets the volume calculation method to one of the supported VOLUME_CALCULATION_METHODS*/
	vtkSetMacroEx(FibresVolumeCalculationMethod, int, VOLUME_CALCULATION_METHOD_CHANGED);
	vtkGetMacro(FibresVolumeCalculationMethod, int);

	/** Turn on/off the computation of PCSA */
	vtkBooleanMacro(ComputePCSA, int);
	vtkSetMacro(ComputePCSA, int);
	vtkGetMacro(ComputePCSA, int);

	/** Sets the proximal line of action to which proximal PA are calculated */	
	vtkSetVector3Macro(ProximalLineOfAction, double);
	
	/** Gets the proximal line of action to which proximal PA are calculated */
	virtual const double* GetProximalLineOfAction();
	virtual void GetProximalLineOfAction(double out[3]);

	/** Sets the proximal line of action to which proximal PA are calculated */
	vtkSetVector3Macro(DistalLineOfAction, double);

	/** Gets the distal line of action to which distal PA are calculated */
	virtual const double* GetDistalLineOfAction();
	virtual void GetDistalLineOfAction(double out[3]);

	/** Gets the position of the proximal line of action., i.e., point through which it comes. */
	virtual const double* GetProximalLineOfActionPos();
	virtual void GetProximalLineOfActionPos(double out[3]);

	/** Gets the position of the distal line of action., i.e., point through which it comes.*/
	virtual const double* GetDistalLineOfActionPos();
	virtual void GetDistalLineOfActionPos(double out[3]);

#ifdef _DEBUG
	/** Gets the proximal line of action to which proximal PA are calculated.
	This line of action is calculated from the attachment points of fibres. */
	virtual const double* GetProximalLineOfActionLR();
	virtual void GetProximalLineOfActionLR(double out[3]);	

	/** Gets the distal line of action to which distal PA are calculated
	This line of action is calculated from the attachment points of fibres. */
	virtual const double* GetDistalLineOfActionLR();
	virtual void GetDistalLineOfActionLR(double out[3]);	
#endif
	
	/** Gets the number of segments of every input fibre. */
	virtual vtkIntArray* GetFibresSegmentsCount();

	/** Gets the lengths of segments of every input fibre. */
	virtual vtkDoubleArray* GetFibresSegmentsLengths();

	/** Gets the lengths of every input fibre. */
	virtual vtkDoubleArray* GetFibresLengths();

	/** Gets the parameterization of every fibre normalized to 0-1. 
	It is chordal or centripetal parameterization depending on UseCentripetalCatmullRom parameter. */
	inline virtual vtkDoubleArray* GetFibresParameterization() {
		return this->UseCentripetalCatmullRom != 0 ? GetFibresCentripetalParameterization() : GetFibresChordalParameterization();
	}

	/** Gets the arc-length (Chordal) parameterization of every fibre normalized to 0-1. */
	virtual vtkDoubleArray* GetFibresChordalParameterization();

	/** Gets the centripetal parameterization of every fibre normalized to 0-1. */
	virtual vtkDoubleArray* GetFibresCentripetalParameterization();

	/** Gets the tangent vector of every fibre point. */
	virtual vtkDoubleArray* GetFibresTangentVectors();
	
	/** Gets the averaged proximal tangent vector of every fibre point. */
	virtual vtkDoubleArray* GetFibresProximalTangentVectors();

	/** Gets the averaged distal tangent vector of every fibre point. */
	virtual vtkDoubleArray* GetFibresDistalTangentVectors();

	/** Gets the proximal pennation angle for every fibre. */
	virtual vtkDoubleArray* GetFibresProximalPA();

	/** Gets the distal pennation angle for every fibre. */
	virtual vtkDoubleArray* GetFibresDistalPA();
	
	/** Gets the radius of every fibre vertex, i.e., radius varies along the fibre. */
	virtual vtkDoubleArray* GetFibresRadii();

	/** Gets the volume of every fibre vertex. */
	virtual vtkDoubleArray* GetFibresVolumes();

	/** Gets the total volume of all fibres. */
	virtual double GetTotalFibresVolume();

	/** Gets the physiological cross-section area. */
	virtual double GetPCSA();
	
	/*virtual*/ void Modified() {
		Superclass::Modified();
	}

	/** Makes the object dirty */
	virtual void Modified(int whatChanged) 
	{
		this->m_SettingsChangedFlags |= whatChanged;
		this->Modified();
	}

protected:	
	vtkLHPMuscleFibresAnalysisFilter();	
	~vtkLHPMuscleFibresAnalysisFilter();

protected:
#if VTK_MAJOR_VERSION < 5
	/**
	This method is the one that should be used by subclasses, right now the
	default implementation is to call the backwards compatibility method */
	/*virtual*/void ExecuteData(vtkDataObject *output);
#else
	// Description:	
	// This is the method of the algorithm in which the algorithm should fill in the output ports
	/*virtual*/ int RequestData(vtkInformation* request,
		vtkInformationVector** inputVector, vtkInformationVector* outputVector);
#endif

	/** Detects the changes and specifies what must be recomputed */
	virtual void MustUpdateCheck(vtkPolyData* inputPoly);

	/** Computes the arc-length parameterization of input fibres and stores them into m_FibresParameterization.
	N.B. CalculateFibresLengths must be called prior this method! */
	void CalculateFibresParameterization();

	/** Computes the Catmull-Rom parameterization of input fibres and stores them into m_FibresCatmullRomParameterization.
	N.B. CalculateFibresSegmentsLengths must be called prior to calling this method*/
	void CalculateFibresCatmullRomParameterization();

	/** Computes the tangent vector at every vertex of every fibre.
	\param poly input fibres (usually, this->GetInput())
	\param params parameterization at every vertex of every fibre.
	\param outTgv output buffer with tangent vectors*/		
	void CalculateFibresTangentVectors(const vtkPolyData* poly, const vtkDoubleArray* params, vtkDoubleArray* outTgv);

	/** Computes the averaged tangent vectors for every fibre either for its proximal or distal part.
	\param poly input fibres (usually, this->GetInput())
	\param tgvecs tangent vectors at every vertex of every fibre
	\param chordalParams chordal parameterization at every vertex of every fibre.
	\param proximal if true, the proximal region of a fibre is used, otherwise the distal one.		
	\param outAvgTgv output buffer with the averaged tangent vectors, one vector per fibre.*/	
	void CalculateFibresAvgTgVec(const vtkPolyData* poly, const vtkDoubleArray* tgvecs,
		const vtkDoubleArray* chordalParams, bool proximal, vtkDoubleArray* outAvgTgv);		

	/** Computes the line of action for proximal or distal part.
	\param poly input fibres (usually, this->GetInput())
	\param avgtgvecs proximal or distal tangent vectors for every fibre
	\param proximal if true, the proximal line of action is calculated, otherwise distal one
	\param la line of action orientation
	\param laPos one point through which the line of action passes
	\param laLR line of action calculated from the attachment points; la == laLR for pennate muscles
	*/
	void CalculateLineOfAction(const vtkPolyData* poly, const vtkDoubleArray* avgtgvecs, 
		bool proximal, double la[3], double laPos[3]
#ifdef _DEBUG
		, double laLR[3]
#endif
		);

	/** Computes the pennation angle for every fibre.
	\param avgtgvecs proximal or distal tangent vectors for every fibre
	\param lineOfAction proximal or distal line of action
	\param outPA output buffer for pennation angles for every fibre*/
	void CalculateFibresPA(const vtkDoubleArray* avgtgvecs, const double lineOfAction[3], vtkDoubleArray* outPA);

	/** Computes the radius of a fibre at each of its vertex.
	\param poly input fibres (usually, this->GetInput())
	\param tgvecs tangent vectors at every vertex of every fibre
	\param outRadii output buffer for the radius of every fibre vertex*/
	void CalculateFibresRadii(const vtkPolyData* poly, const vtkDoubleArray* tgvecs, vtkDoubleArray* outRadii);

	/** Computes the volume of every fibre using the method described by Lee et al. 2014.
	\param poly input fibres (usually, this->GetInput())
	\param segLengths segment lengths of every fibre
	\param outVolume output buffer for the volume of every fibre*/
	void CalculateFibresVolumeLee14(const vtkPolyData* poly, const vtkDoubleArray* segLengths, vtkDoubleArray* outVolume);

	/** Computes the volume of every fibre using the own proprietal method.
	\param poly input fibres (usually, this->GetInput())
	\param segLengths segment lengths of every fibre
	\param fibRadii the radii at every vertex of every fibre
	\param outVolume output buffer for the volume of every fibre*/
	void CalculateFibresVolumeCylindrical(const vtkPolyData* poly, const vtkDoubleArray* segLengths,
		const vtkDoubleArray* fibRadii, vtkDoubleArray* outVolume);


	/** Computes the volume of every fibre using the method described by Ravichandrian et al. 2009.
	\param poly input fibres (usually, this->GetInput())
	\param fibLengths total lengths of every fibre
	\param fibRadii the radii at every vertex of every fibre
	\param outVolume output buffer for the volume of every fibre*/
	void CalculateFibresVolumeRavichandrian(const vtkPolyData* poly, const vtkDoubleArray* fibLengths, 
		const vtkDoubleArray* fibRadii, vtkDoubleArray* outVolume);

protected:
	/** Computes centroid of points */
	void ComputeCentroid(const vtkPoints* points, VCoord centroid);	//TODO: refactor this, duplication of code - see vtkMAFMuscleDecomposition

	/** Computes three eigen vectors of the given point set.
	Centroid may be NULL, if it should be calculated automatically, otherwise, it must be centroid of points.
	The returned eigen vectors are ordered by their lengths. The longest (i.e., the principal axis) is denoted by the first one.*/
	void ComputeEigenVects(const vtkPoints* points, const VCoord centroid, VCoord eigenvects[3]); 	//TODO: refactor this, duplication of code - see vtkMAFMuscleDecomposition	

	/** Gets the distance between the given point and its nearest fibre.
	The calculation is done on a plane defined by the given point (pt) and normal (normal) and
	the fibre with excludeFibreId identifier is excluded from the calculation. 
	The method returns 0.0, if there is no fibre intersected by the plane.
	excludeFibreId may by -1, if every fibres should be tested. */
	double GetNearestFibreDistance(const vtkPolyData* poly, const VCoord pt, 
		const VCoord normal, int excludeFibreId);
	
	/** Gets the square distance between the given point and the fibre. */
	double Distance2BetweenFibreAndPoint(const vtkPolyData* poly, int fibreId, const VCoord pt);
protected:	
	unsigned long m_LastExecuteInputTimeStamp;	//<timestamp for the last successful execute	
	int m_SettingsChangedFlags;									//<contains what has changed since the last execute	
	int m_CachedDataFlags;											//<what is currently cached (updated by MustUpdateCheck)

	int ComputeFibresLengths;						///<non-zero to compute fibres length
	int ComputeFibresParameterization;	///<non-zero to compute 0-1 arc-length parameterization for each fibre	
	int ComputeFibresTangentVector;			///<non-zero to compute tangent vectors for each node of every fibre

	int ComputeProximalLineOfAction;		///<non-zero to estimate proximal line of action from the fibres shape
	int ComputeDistalLineOfAction;			///<non-zero to estimate distal line of action from the fibres shape

	int ComputeProximalPennationAngle;	///<non-zero to calculate pennation angle (PA) of every fibre to the proximal line of action
	int ComputeDistalPennationAngle;		///<non-zero to calculate pennation angle (PA) of every fibre to the distal line of action

	int ComputeFibresRadii;						///<non-zero to calculate the radii of every fibre (at every its vertex)

	int ComputeFibresVolumes;					///<non-zero to calculate volumes of fibres (fascicles) in the muscle
	int ComputeTotalFibresVolume;			///<non-zero to calculate the volume of all fibres (fascicles) in the muscle (sum of Fibres Volumes)


	int ComputePCSA;									///<non-zero to calculate the physiological cross-section area (PCSA)

	int UseCentripetalCatmullRom;		///<if non-zero, the centripetal Catmull-Rom parameterization will be used, 
																	///otherwise chordal parameterization will be used	

	int FibresVolumeCalculationMethod;	///<one of CALCULATE_VOLUME_METHODS
	
	double ProximalThreshold;				///<Fraction  of the fibre used for calculating the proximal line of action (0..1)
	double DistalThreshold;					///<Fraction  of the fibre used for calculating the distal line of action (0..1)

	//NOTE: these arrays are not to be accessed directly but via appropriate Get methods that ensures its validity
	vtkIntArray* m_FibresSegmentsCount;							///<Cached Fibres number of segments each fibre consists of
	vtkDoubleArray* m_FibresSegmentsLengths;		///<Cached Fibres lengths	

	vtkDoubleArray* m_FibresLengths;					///<Cached Fibres lengths	- fibre = poly-line	
	vtkDoubleArray* m_FibresParameterization;	///<Cached Fibres parameterization
	vtkDoubleArray* m_FibresCRParameterization;///<Cached Fibres parameterization for Catmull-Rom
	vtkDoubleArray* m_FibresTangentVectors;			///<Cached Fibres tangent vectors for every point	
	vtkDoubleArray* m_FibresProximalAvgTgVec;	///<Cached the averaged tangent vectors for proximal and distal region of every fibre
	vtkDoubleArray* m_FibresDistalAvgTgVec;		///<Cached the averaged tangent vectors for proximal and distal region of every fibre
	vtkDoubleArray* m_FibresProximalPA;			///<Cached the proximal pennate angle (PA) for every fibre
	vtkDoubleArray* m_FibresDistalPA;				///<Cached the distal pennate angle (PA) for every fibre
	vtkDoubleArray* m_FibresRadii;					///<Cached the radii for every fibre
	vtkDoubleArray* m_FibresVolume;					///<Cached the volume for every fibre
	double m_TotalVolume;										///<Cached the total volume
	double m_PCSA;													///<Cached the PCSA

	double ProximalLineOfAction[3];		///<proximal line of action, either calculated or specified
	double DistalLineOfAction[3];			///<distal line of action, either calculated or specified

	double ProximalLineOfActionPos[3];	///<position of the proximal line of action calculated from the attachment points of fibres
	double DistalLineOfActionPos[3];		///<position of the distal line of action calculated from the attachment points of fibres

#ifdef _DEBUG
	double ProximalLineOfActionLR[3];		///<proximal line of action calculated from the attachment points of fibres
	double DistalLineOfActionLR[3];			///<distal line of action calculated from the attachment points of fibres	
#endif
	
private:
	vtkLHPMuscleFibresAnalysisFilter(const vtkLHPMuscleFibresAnalysisFilter&);  // Not implemented.
	void operator = (const vtkLHPMuscleFibresAnalysisFilter&);					// Not implemented.  
};


#if VTK_MAJOR_VERSION < 5
#undef vtkPolyDataAlgorithm
#endif

#endif //__vtkLHPMuscleFibresAnalysisFilter_h
