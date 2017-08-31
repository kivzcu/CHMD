#ifndef __Fiber_h
#define __Fiber_h

#include "vtkstd/vector"
#include "vtkstd/algorithm"

class Fiber
{
public:
	vtkstd::vector<int>  PointIDs;
	int Id;
	bool BBoundary;
	Fiber(void);
	~Fiber(void);
};
#endif