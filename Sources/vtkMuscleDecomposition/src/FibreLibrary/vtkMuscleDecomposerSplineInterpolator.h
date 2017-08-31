#pragma once
#include "vtkMuscleDecomposer.h"

class vtkMuscleDecomposer::SplineInterpolator
{
	public:
		SplineInterpolator(int splineType);
		~SplineInterpolator();

		void Interpolate(vtkMuscleDecomposer::Fiber* fiber, std::vector<vtkMuscleDecomposer::Point*> & points);

	private:
		int splineType;

		vtkMuscleDecomposer::Fiber* fiber;   // Pointer to fiber data
		std::vector<vtkMuscleDecomposer::Point*> pointCloud; // Point data
		int n; // number of points

		std::vector<double> subdiagonal;			// Diagonal under main diagonal (for all axes)
		std::vector<double> superdiagonal;		// Diagonal above main diagonal (for all axes)
		std::vector<double> diagonal[3];			// Main diagonal (for 3 axes - modified by calculations)
		std::vector<double> result[3];			// Resulting vector, containing S_i..S_n (for 3 axes)
		std::vector<double> rightHand[3];		// Right hand vector (for 3 axes)

		std::vector<double> derivations[3];		// Derivations for Constrained Splines

		void BuildParameter(bool useCentripetal);	// Builds independent variable (parameter) t
		void BuildEquation();						// Creates sub, sup and diagonal vectors and righthand vectors
		void SolveEquation();						// Solves equation and associates weights to points
		void ComputeWeights();						// Computes weights for all points using eq. result

		void ComputeFibreDerivations();				// Computes derivations for contrained splines
		void ComputeConstrainedWeights();			// Computes weights for constrained splines

		void ComputeCatmullRom();					// Prepares links for Catmull-Rom splines
	public:
		int GetInterpolationType();

		
};

