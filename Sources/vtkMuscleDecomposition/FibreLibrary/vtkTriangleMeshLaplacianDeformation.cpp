/*=========================================================================
Program: Musculoskeletal Modeling (VPHOP WP10)
Module: vtkTriangleMeshLaplacianDeformation.cpp

Authors: Josef Kohout
==========================================================================
Copyright (c) 2014 University of West Bohemia (www.zcu.cz)
See the COPYINGS file for license details
=========================================================================
*/
#include "vtkTriangleMeshLaplacianDeformation.h"
#include "vtkObjectFactory.h"
#include "vtkPolyData.h"
#include "vtkCellArray.h"
#include "vtkIdList.h"
#include "vtkMath.h"

extern "C"
{
#include "csparse.h"
}


vtkStandardNewMacro(vtkTriangleMeshLaplacianDeformation);


vtkTriangleMeshLaplacianDeformation::vtkTriangleMeshLaplacianDeformation()
{
	this->FixedPointsIds = NULL;
	this->MovedPointsIds = NULL;
	this->MovedPointsCoords = NULL;
}


vtkTriangleMeshLaplacianDeformation::~vtkTriangleMeshLaplacianDeformation()
{
	if (this->FixedPointsIds != NULL)
	{
		this->FixedPointsIds->Delete();
		this->FixedPointsIds = NULL;
	}

	if (this->MovedPointsIds != NULL)
	{
		this->MovedPointsIds->Delete();
		this->MovedPointsIds = NULL;
	}

	if (this->MovedPointsCoords != NULL)
	{
		this->MovedPointsCoords->Delete();
		this->MovedPointsCoords = NULL;
	}
}


//Constructs sparse Laplacian matrix for the input triangular mesh
void vtkTriangleMeshLaplacianDeformation::ConstructMatrixL(
	vtkPolyData* inputPoly, std::vector< RowMap >& matrixL) const
{
	//Force calling inputPoly->BuildCells, if it has not been already built 	
	//to be able to call GetCellPoints
	inputPoly->GetCellType(0);

	int nPoints = inputPoly->GetNumberOfPoints();
	matrixL.resize(nPoints);
	
	int nCells = inputPoly->GetNumberOfCells();
	//#pragma omp parallel for shared(cpu_sparse_matrix, nCells, surface, processed, nPoints)
	for (int i = 0; i < nCells; i++)
	{
		vtkIdType nPts, *pPts;
		inputPoly->GetCellPoints(i, nPts, pPts);

		//get coordinates we will need for our calculation			
		VCoord Coords[3];

		for (int j = 0; j < 3; j++) {
			inputPoly->GetPoint(pPts[j], Coords[j].xyz);
		}

		//compute cotans and add them to cpu_sparse_matrix
		for (int j = 0; j < 3; j++)
		{
			//cotg(alfa) = cos(alfa)/sin(alfa)
			//cos(alfa) = u*v /(|u|*|v|), where u and v are vectors of the triangle
			//sin(alfa) = |u x v| / (|u|*|v|)
			//thus cotg(alfa) = u*v / |u x v|

			vtkIdType iv1 = (j + 1) % 3;
			vtkIdType iv2 = (j + 2) % 3;

			double u[3], v[3];
			for (int k = 0; k < 3; k++)
			{
				u[k] = Coords[iv1].xyz[k] - Coords[j].xyz[k];
				v[k] = Coords[iv2].xyz[k] - Coords[j].xyz[k];
			}

			double uxv[3];
			vtkMath::Cross(u, v, uxv);
			double cotan_2 = 0.5*vtkMath::Dot(u, v) / vtkMath::Norm(uxv);

			//the calculated cotan is used in iv1-iv2 edge				
			//#pragma omp critical 
			{
				matrixL[pPts[iv1]][pPts[iv2]] += cotan_2;
				matrixL[pPts[iv2]][pPts[iv1]] += cotan_2;
			}
		}
	} //end for cells	

	//normalization and adding 1.0 on diagonals
	//#pragma omp parallel for shared(cpu_sparse_matrix)
	for (int i = 0; i < nPoints; i++)
	{
		double dblSum = 0.0;
		for (RowMap::const_iterator it = matrixL[i].begin();
			it != matrixL[i].end(); it++) {
			dblSum += (*it).second;
		}

		for (RowMap::iterator it = matrixL[i].begin();
			it != matrixL[i].end(); it++)
		{
			(*it).second /= -dblSum;	//weight
		}

		matrixL[i][i] = 1.0;	//add diagonal			
	}
}

//Constructs the vector B = A^T*b where b is n+m+o x 1 vector of right side
void vtkTriangleMeshLaplacianDeformation::ConstructVectorB(vtkPolyData* inputPoly,
	const std::vector< RowMap >& matrixL, std::vector< VCoord >& vectorB) const
{
	int nPoints = inputPoly->GetNumberOfPoints();
	vectorB.resize(nPoints);	//zeroes	

	//B'_j = A^T_j,i*b_i =
	//L^T_j,i*d_i	[i = 0 .. n - 1] + 	F^T_j,i*x0_i	[i = 0 .. f - 1] +	M^T_j,i*xm_i	[i = 0 .. m - 1] =
	//L_i,j*d_i		[i = 0 .. n - 1] + F_i,j*x0_i	[i = 0 .. f - 1] +	M_i,j*xm_i	[i = 0 .. m - 1] =>

	//#pragma omp parallel for shared(cpu_sparse_matrix)
	for (int i = 0; i < nPoints; i++)
	{
		//calculate d_i = (V_i + sum(L_i,j*V_j))
		VCoord di = { 0, 0, 0 };

		//process V_j
		for (RowMap::const_iterator it = matrixL[i].begin();
			it != matrixL[i].end(); it++)
		{
			VCoord Vj;
			inputPoly->GetPoint((*it).first, Vj.xyz);
			for (int k = 0; k < 3; k++) {
				di.xyz[k] += Vj.xyz[k] * (*it).second;
			}
		}

		//calculate b'j += L_i,j*d_i
		for (RowMap::const_iterator it = matrixL[i].begin();
			it != matrixL[i].end(); it++)
		{
			for (int k = 0; k < 3; k++) {
				vectorB[(*it).first].xyz[k] += (*it).second * di.xyz[k];
			}
		}		
	}

	//append the fixed points
	//F_i,j != 0 only for one element
	if (this->FixedPointsIds != NULL)
	{
		int nFixedPoints = this->FixedPointsIds->GetNumberOfIds();
		vtkIdType* pFIds = this->FixedPointsIds->GetPointer(0);
		for (int i = 0; i < nFixedPoints; i++)
		{
			VCoord x0i;
			inputPoly->GetPoint(pFIds[i], x0i.xyz);

			for (int k = 0; k < 3; k++) {
				vectorB[pFIds[i]].xyz[k] += x0i.xyz[k];
			}
		}
	}

	//append the moved points
	//F_i,j != 0 only for one element
	if (this->MovedPointsIds != NULL && this->MovedPointsCoords != NULL)
	{
		int nMovedPoints = this->MovedPointsIds->GetNumberOfIds();
		vtkIdType* pMIds = this->MovedPointsIds->GetPointer(0);
		for (int i = 0; i < nMovedPoints; i++)
		{
			VCoord xmi;
			this->MovedPointsCoords->GetPoint(i, xmi.xyz);

			for (int k = 0; k < 3; k++) {
				vectorB[pMIds[i]].xyz[k] += xmi.xyz[k];
			}
		}
	}
}

// Constructs sparse matrix A^T*A
//N.B. this method assumes matrixL is correctly filled by ConstructMatrixL
void vtkTriangleMeshLaplacianDeformation::ConstructMatrixAtA(
	vtkPolyData* inputPoly, const std::vector< RowMap >& matrixL, std::vector< RowMap >& matrixQ) const
{
	int nPoints = inputPoly->GetNumberOfPoints();
	matrixQ.resize(nPoints);

	//		|	L	|	 (n x n) with 1 at diagonals and non-zero values only at i, j for vertices V_i, V_j connected by an edge
	//A =	|	F	|	 (m x n) with 1 at only one column at each row	
	//		|	M	|	 (o x n) with 1 at only one column at each row	

	//Q_i,j = sum(L^T_i,k*L_k,j) + sum(F^T_i,k*F^T_k,j) + sum(M^T_i,k*M^T_k,j) =
	//sum(L_k,i*L_k,j) + sum(F_k,i*F_k,j) + sum(M_k,i*M_k,j) => the matrix is symmetric
	for (int k = 0; k < nPoints; k++)
	{
		for (RowMap::const_iterator it = matrixL[k].begin();
			it != matrixL[k].end(); it++)
		{
			//L_k,i is (*it).second, i is (*it).first
			for (RowMap::const_iterator it2 = matrixL[k].begin();
				it2 != matrixL[k].end(); it2++)
			{
				//L_k,j is (*it2).second, j is (*it2).first
				double Lkij = (*it).second * (*it2).second;
				matrixQ[(*it).first][(*it2).first] += Lkij;
				matrixQ[(*it2).first][(*it).first] += Lkij;
			}
		}
	}

	//append the fixed points
	//F_k,i*F_k,j != 0 only for i == j because in one row of F there is only one non-zero value
	//and this non-zero value is in the column corresponding to the index of k-th fixed point
	if (this->FixedPointsIds != NULL)
	{
		int nFixedPoints = this->FixedPointsIds->GetNumberOfIds();
		vtkIdType* pFIds = this->FixedPointsIds->GetPointer(0);
		for (int i = 0; i < nFixedPoints; i++) {
			matrixQ[pFIds[i]][pFIds[i]] += 1;
		}
	}

	//append the moved points in a similar fashion as the fixed points	
	if (this->MovedPointsIds != NULL && this->MovedPointsCoords != NULL)
	{
		int nMovedPoints = this->MovedPointsIds->GetNumberOfIds();
		vtkIdType* pMIds = this->MovedPointsIds->GetPointer(0);
		for (int i = 0; i < nMovedPoints; i++) {
			matrixQ[pMIds[i]][pMIds[i]] += 1;
		}
	}
}

#if VTK_MAJOR_VERSION < 5
#define VTK6_RESULT_OK
#define VTK6_RESULT_FAIL
/**
Executes the data operation.

@param [in,out]	output	If non-null, the output.
*/
/*virtual*/ void vtkTriangleMeshLaplacianDeformation::ExecuteData(vtkDataObject *output)
{
	vtkPolyData* inputPoly = GetInput();
	vtkPolyData* outputPoly = vtkPolyData::SafeDownCast(output);
#else

#define VTK6_RESULT_OK		1
#define VTK6_RESULT_FAIL	0

// This is the method of the algorithm in which the algorithm should fill in the output ports
/*virtual*/ int vtkTriangleMeshLaplacianDeformation::RequestData(vtkInformation* request,
	vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
	vtkPolyData* inputPoly = vtkPolyData::GetData(inputVector[0]);
	vtkPolyData* outputPoly = vtkPolyData::GetData(outputVector);
#endif

	//check whether output is valid	
	if (inputPoly == NULL || inputPoly->GetPoints() == NULL ||
		inputPoly->GetPoints()->GetNumberOfPoints() == 0)
	{
		vtkErrorMacro(<< "Invalid input for vtkTriangleMeshLaplacianDeformation.");
		return VTK6_RESULT_FAIL;   //we have no valid input
	}
	
	if (outputPoly == NULL)
	{
		vtkErrorMacro(<< "Invalid output for vtkTriangleMeshLaplacianDeformation.");
		return VTK6_RESULT_FAIL;   //we have no valid output
	}	

	//copy input to the output
	outputPoly->ShallowCopy(inputPoly);

	if ((this->FixedPointsIds == NULL || this->FixedPointsIds->GetNumberOfIds() == 0) &&
		(this->MovedPointsIds == NULL || this->MovedPointsIds->GetNumberOfIds() == 0))
	{
		//no change so we are ready			
		return VTK6_RESULT_OK;
	}

	//we need to construct a system of linear equations A*x = b:
	//	|	L	|		| d  |
	//	|	F	| * X = | x0 |
	//	|	M	|		| xm |
	//where L is n x n Laplacian submatrix containing the shape constraints such 
	//that wj*Vj - Vi = di whereas Vj are vertices adjacent to Vi,
	//i.e., Lij = 1, if i == j, Lij = -wij, if there is an edge between Vi and Vj, Lij = 0 otherwise
	//F contains only fixed points constraints, i.e., Fij = 1, if Vj is fixed to x0j
	//M contains only moved points constraints, i.e., Mij = 1, if Vj has moved to xmj
		
	//A*x = b is over-constrained and, therefore, it will be solved by least-squares method, 
	//which means we convert the equation into A^T*A * X = A^T*b -> Q * X = B

	//first, construct the matrix L (sparse one)
	std::vector< RowMap > cpu_sparse_matrixL;
	ConstructMatrixL(inputPoly, cpu_sparse_matrixL);

	//next, construct the vector B
	std::vector< VCoord > B;
	ConstructVectorB(inputPoly, cpu_sparse_matrixL, B);
	
	//finally construct A^T*A
	std::vector< RowMap > matrixQ;
	ConstructMatrixAtA(inputPoly, cpu_sparse_matrixL, matrixQ);
	
	//determine the number of non-zero values
	size_t nonzeroes = 0, nPoints = matrixQ.size();
	for (size_t i = 0; i < nPoints; i++) {
		nonzeroes += matrixQ[i].size();
	}
	
	//process each dimension (x,y,z) independently
	double* pValues = new double[nPoints];
	for (int dim = 0; dim < 3; dim++)
	{
		//copy the vector of right side
		for (size_t i = 0; i < nPoints; i++) {
			pValues[i] = B[i].xyz[dim];
		}

		//convert matrixQ to CS sparse matrix
		//allocate a sparse matrix (in compressed column format)
		cs* A = cs_spalloc((int)nPoints, (int)nPoints, (int)nonzeroes, 1, 0);

		//since cs is column ordered, let us find the number of nnz in every column
		memset(A->p, 0, (nPoints + 1) * sizeof(A->p[0]));

		//#pragma omp parallel for shared(cpu_sparse_matrix, A)
		for (int i = 0; i < nPoints; i++)
		{
			for (RowMap::const_iterator it = matrixQ[i].begin();
				it != matrixQ[i].end(); it++) {
				A->p[(*it).first + 1]++;
			}
		}

		//calculate column pointers (0, n1, n1+n2, n1+n2+n3, ...)
		for (int j = 1; j <= nPoints; j++) {
			A->p[j] += A->p[j - 1];
		}

		//store values
		for (int i = 0; i < nPoints; i++)
		{
			for (RowMap::const_iterator it = matrixQ[i].begin();
				it != matrixQ[i].end(); it++)
			{
				int index = A->p[(*it).first];	//column index ptr
				A->i[index] = i;								//row index
				A->x[index] = (*it).second;			//value
				A->p[(*it).first]++;						//increase column index
			}
		}

		//column pointers must be shifted down (n1, n1+n2, n1+n2+n3, ...) -> (0, n1, n1+n2, ...)
		for (size_t j = nPoints - 1; j > 0; j--) {
			A->p[j] = A->p[j - 1];
		}

		A->p[0] = 0;

		int result = cs_lusol(A, pValues, 0, 1e-8);

		cs_spfree(A);	//free the matrix

		//store the vector of right side
		for (size_t i = 0; i < nPoints; i++) {
			B[i].xyz[dim] = pValues[i];
		}
	}

	delete[] pValues;

	//set the points to the output
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	points->SetNumberOfPoints(nPoints);
	for (size_t i = 0; i < nPoints; i++) {
		points->SetPoint(i, B[i].xyz);
	}

	outputPoly->SetPoints(points);	
	return VTK6_RESULT_OK;
}