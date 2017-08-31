#pragma once

#include "vtkRenderer.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkProperty.h"
#include <vector>

class vtkMapperActor
{
public:
	vtkMapperActor(vtkPolyData* data, int R, int G, int B, int A, 
		bool wireframe, vtkRenderer* ren1, int groupId = -1, bool show = true);
	~vtkMapperActor(void);

	static void DeleteMAs(int groupId);	//deletes MAs of the specified group
	static void DeleteAllMAs();					//deletes all MAs

	static void ShowHideMAs(int groupId, bool show);	//Shows or hides MAs of the specified group	

	inline vtkActor* GetActor() {
		return actor;
	}

private:
	static std::vector<vtkMapperActor*> AllMAs;

	vtkRenderer* renderer;
	vtkActor* actor;
	
	int groupId;
};

