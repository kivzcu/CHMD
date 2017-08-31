#pragma once

void RegisterFibresAndMesh(vtkPolyData * fibresData, vtkPolyData * muscleData
#ifdef _DEBUG
	, const char* path
#endif
	);

/** Constructs the Catmull-Rom surface patch interpolating the input fibres
ensuring that the orientation of its cells and its normals is consistent
with the orientation of the cells (and normals) of the muscle surface. */
void ConstructFibresPatch(vtkPolyData* fibresData, vtkPolyData* muscleData, vtkPolyData* out);

/** Moves the fibres patch in the given direction to make sure that 
the transformed patch does not intersect the muscle. */
void TransformFibresPatch(vtkPolyData* fibresPoly, vtkPolyData* muscleData, 
	const double dir[3], vtkPolyData* out);

/** Returns true if the orientation of the normals of the fibres patch surface 
is consistent with the orientation of the normals of the muscle surface. */
bool FibresPatchNormalsOK(vtkPolyData* fibresPoly, vtkPolyData* musclePoly);

/** Computes the normals at points of the given polydata.
N.B. the caller is responsible for deleting the returned object.*/
vtkDataArray* ComputePointNormals(vtkPolyData* polyData);

/** Computes the average vector (normal) of the given vector (normal) field. */
void GetAvgNormal(vtkDataArray* normals, double avgNormal[3]);
