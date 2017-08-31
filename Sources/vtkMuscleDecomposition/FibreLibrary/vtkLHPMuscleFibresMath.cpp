#include "vtkLHPMuscleFibresMath.h"
#include "vtkPolyData.h"
#include "vtkCellArray.h"
#include "vtkMath.h"
#include "vtkDoubleArray.h"
#include "vtkIntArray.h"


/** Counts the segments of every fibre in poly and stores this information into out */
/*static*/ void vtkLHPMuscleFibresMath::CountFibresSegments(const vtkPolyData* poly, vtkIntArray* out)
{
	int nFibres = const_cast<vtkPolyData*>(poly)->GetNumberOfCells();
	CountFibresSegments(poly, out->WritePointer(0, nFibres));
}

/** Computes the lengths of each segment of the fibres in poly and stores them into out. */
/*static*/ void vtkLHPMuscleFibresMath::CalculateFibresSegmentsLengths(const vtkPolyData* poly, vtkDoubleArray* out)
{
	int nPoints = const_cast<vtkPolyData*>(poly)->GetNumberOfPoints();
	int nFibres = const_cast<vtkPolyData*>(poly)->GetNumberOfCells();
	
	CalculateFibresSegmentsLengths(poly, out->WritePointer(0, nPoints - nFibres));
}

/** Computes the overall lengths of input fibres and stores them into out.	*/
/*static*/ void vtkLHPMuscleFibresMath::CalculateFibresLengths(const vtkIntArray* segmentsCount, const vtkDoubleArray* segmentsLengths, vtkDoubleArray* out)
{
	int nFibres = const_cast<vtkIntArray*>(segmentsCount)->GetNumberOfTuples();
	CalculateFibresLengths(nFibres, const_cast<vtkIntArray*>(segmentsCount)->GetPointer(0),
		const_cast<vtkDoubleArray*>(segmentsLengths)->GetPointer(0), out->WritePointer(0, nFibres));
}

/** Computes either chordal or centripetal parameterization of the input fibres and stores it into out.
Note: the parameterization of a fibre is from 0 to t_m-1. */
/*static*/ void vtkLHPMuscleFibresMath::CalculateFibresParameterization(
	const vtkPolyData* poly, const vtkDoubleArray* segmentsLengths,
	bool centripetal, vtkDoubleArray* out)
{
	CalculateFibresParameterization(poly, const_cast<vtkDoubleArray*>(segmentsLengths)->GetPointer(0),
		centripetal, out->WritePointer(0, const_cast<vtkPolyData*>(poly)->GetNumberOfPoints()));
}

/** Normalizes the given parameterization of fibres into 0 - 1.
It is safe to call this method with out being parms for in - place normalization. * /
/*static*/ void vtkLHPMuscleFibresMath::NormalizeFibresParameterization(
	const vtkPolyData* poly,	const vtkDoubleArray* params, vtkDoubleArray* out)
{
	NormalizeFibresParameterization(poly, const_cast<vtkDoubleArray*>(params)->GetPointer(0),
		out->WritePointer(0, const_cast<vtkDoubleArray*>(params)->GetNumberOfTuples()));
}

/** Counts the segments of every fibre in poly and stores this information into pnVerts.
N.B. the buffer must be big enough to hold poly->GetNumberOfCells() items. */
/*static*/ void vtkLHPMuscleFibresMath::CountFibresSegments(const vtkPolyData* poly, int* pnVerts)
{
	int nFibres = const_cast<vtkPolyData*>(poly)->GetNumberOfCells();	
	for (int i = 0; i < nFibres; i++)
	{
		vtkIdType nPts, *pPts;
		const_cast<vtkPolyData*>(poly)->GetCellPoints(i, nPts, pPts);

		*pnVerts = nPts - 1;	//number of segments is number of vertices - 1
		pnVerts++;
	}
}

/** Computes the lengths of each segment of the fibres in poly and stores them into pdblLen.
N.B. the buffer must be big enough to hold poly->GetNumberOfPoints() - GetNumberOfCells() items. */
/*static*/ void vtkLHPMuscleFibresMath::CalculateFibresSegmentsLengths(const vtkPolyData* poly, double* pdblLen)
{
	int nFibres = const_cast<vtkPolyData*>(poly)->GetNumberOfCells();
	int nSegments = const_cast<vtkPolyData*>(poly)->GetNumberOfPoints() - nFibres;
	
	for (int i = 0; i < nFibres; i++)
	{
		vtkIdType nPts, *pPts;
		const_cast<vtkPolyData*>(poly)->GetCellPoints(i, nPts, pPts);

		double AB[2][3];	//buffer containing both end-points of the segment
		const_cast<vtkPolyData*>(poly)->GetPoint(pPts[0], AB[0]);

		int iFirstPt = 0;
		for (int j = 1; j < nPts; j++)
		{
			const_cast<vtkPolyData*>(poly)->GetPoint(pPts[j], AB[1 - iFirstPt]);	//get the next point		

			*pdblLen = sqrt(vtkMath::Distance2BetweenPoints(AB[0], AB[1]));
			pdblLen++;

			iFirstPt = 1 - iFirstPt;
		}
	}
}

/** Computes the overall lengths of input fibres and stores them into pdblTotLen.
N.B. the output buffer must be capable to hold :nFibres: items. */
/*static*/ void vtkLHPMuscleFibresMath::CalculateFibresLengths(int nFibres, const int* pdblSegCount,
	const double* pdblSegLen, double* pdblTotLen)
{	
	for (int i = 0; i < nFibres; i++)
	{
		double dblLen = 0.0;
		for (int j = 0; j < (*pdblSegCount); j++)
		{
			dblLen += *pdblSegLen;
			pdblSegLen++;
		}

		*pdblTotLen = dblLen;
		pdblTotLen++;
		pdblSegCount++;
	}
}

/** Computes either chordal or centripetal parameterization of the input fibres and stores it into out.
N.B. the output buffer must be capable holding nPoints (poly->GetNumberOfPoints()) parameters. */
/*static*/ void vtkLHPMuscleFibresMath::CalculateFibresParameterization(
	const vtkPolyData* poly, const double* pdblSegLen, bool centripetal, double* pPars)
{
	int nFibres = const_cast<vtkPolyData*>(poly)->GetNumberOfCells();
	for (int i = 0; i < nFibres; i++)	//cannot be parallelized
	{
		vtkIdType nPts, *pPts;
		const_cast<vtkPolyData*>(poly)->GetCellPoints(i, nPts, pPts);

		pPars[pPts[0]] = 0.0;

		double dblLen = 0.0;
		for (int j = 1; j < nPts; j++)
		{
			if (!centripetal)
				dblLen += *pdblSegLen;
			else
				dblLen += sqrt(*pdblSegLen);

			pdblSegLen++;

			pPars[pPts[j]] = dblLen;
		}
	}
}

/** Normalizes the given parameterization of fibres into 0 - 1.
It is safe to call this method with out being parms for in-place normalization.
N.B. he output buffer must be capable holding nPoints (poly->GetNumberOfPoints()) parameters. */
/*static*/ void vtkLHPMuscleFibresMath::NormalizeFibresParameterization(
	const vtkPolyData* poly, const double* params, double* out)
{
	int nFibres = const_cast<vtkPolyData*>(poly)->GetNumberOfCells();

#pragma omp parallel for shared(out)
	for (int i = 0; i < nFibres; i++)
	{
		vtkIdType nPts, *pPts;
		const_cast<vtkPolyData*>(poly)->GetCellPoints(i, nPts, pPts);

		double dblTotLen = params[pPts[nPts - 1]];
		out[pPts[nPts - 1]] = 1.0;

		for (int j = 1; j < nPts - 1; j++) {
			out[pPts[j]] = params[pPts[j]] / dblTotLen;
		}
	}
}

/*static*/ double vtkLHPMuscleFibresMath::Sum(const vtkDoubleArray* inArray)
{
	return Sum(const_cast<vtkDoubleArray*>(inArray)->GetPointer(0),
		const_cast<vtkDoubleArray*>(inArray)->GetNumberOfTuples());
}

/*static*/ double vtkLHPMuscleFibresMath::Avg(const vtkDoubleArray* inArray)
{
	return Avg(const_cast<vtkDoubleArray*>(inArray)->GetPointer(0),
		const_cast<vtkDoubleArray*>(inArray)->GetNumberOfTuples());
}

/*static*/ double vtkLHPMuscleFibresMath::Min(const vtkDoubleArray* inArray)
{
	return Min(const_cast<vtkDoubleArray*>(inArray)->GetPointer(0),
		const_cast<vtkDoubleArray*>(inArray)->GetNumberOfTuples());
}

/*static*/ double vtkLHPMuscleFibresMath::Max(const vtkDoubleArray* inArray)
{
	return Max(const_cast<vtkDoubleArray*>(inArray)->GetPointer(0),
		const_cast<vtkDoubleArray*>(inArray)->GetNumberOfTuples());
}

/*static*/ double vtkLHPMuscleFibresMath::Dev(const vtkDoubleArray* inArray)
{
	return Dev(const_cast<vtkDoubleArray*>(inArray)->GetPointer(0),
		const_cast<vtkDoubleArray*>(inArray)->GetNumberOfTuples());
}

/*static*/ double vtkLHPMuscleFibresMath::Dev(const double* values, int count)
{
	if (count == 1)
		return 0.0;

	double mean = Avg(values, count);
	double dev = 0.0;

	for (int i = 0; i < count; i++) {
		dev += (values[i] - mean)*(values[i] - mean);
	}
		
	return sqrt(dev / (count - 1));
}
