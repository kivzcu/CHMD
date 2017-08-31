#include "vtkMapperActor.h"
#include "vtkRenderer.h"

#include "memdebug.h"

vtkMapperActor::vtkMapperActor(vtkPolyData* data, int R, int G, int B, int A, 
	bool wireframe, vtkRenderer* ren1, int groupId, bool show)
{
	vtkPolyDataMapper* mapper = vtkPolyDataMapper::New();
#if VTK_MAJOR_VERSION < 5
	mapper->SetInput(data);
#else
	mapper->SetInputData(data);
#endif
	mapper->SetScalarVisibility(0);

	this->actor = vtkActor::New();
	this->actor->SetMapper(mapper);
	this->actor->GetProperty()->SetDiffuseColor(R / 255.0, G / 255.0, B/255.0);
	this->actor->GetProperty()->SetOpacity(A / 255.0);
	if (wireframe)
	{
		this->actor->GetProperty()->SetRepresentationToWireframe();
	}

	this->groupId = groupId;
	
	this->actor->SetVisibility(show ? 1 : 0);
	ren1->AddActor( this->actor );	

	this->renderer = ren1;
	this->renderer->Register(NULL);
	AllMAs.push_back(this);

	mapper->Delete();
}

vtkMapperActor::~vtkMapperActor(void)
{	
	//this->renderer->RemoveActor(this->actor);
	this->actor->SetVisibility(0);
	this->renderer->UnRegister(NULL);
	this->actor->Delete();
}

std::vector<vtkMapperActor*> vtkMapperActor::AllMAs;

void vtkMapperActor::DeleteAllMAs()
{
	for (int i = 0; i < (int)AllMAs.size(); i++)
	{
		delete AllMAs[i];
	}
	AllMAs.clear();
}

void vtkMapperActor::DeleteMAs(int groupId)
{
	for (int i = 0; i < (int)AllMAs.size(); i++)
	{
		if (AllMAs[i]->groupId == groupId){
			delete AllMAs[i];

			AllMAs.erase(AllMAs.begin() + i);
		}
	}	
}

void vtkMapperActor::ShowHideMAs(int groupId, bool show)
{
	for (int i = 0; i < (int)AllMAs.size(); i++)
	{
		if (AllMAs[i]->groupId == groupId){
			AllMAs[i]->actor->SetVisibility(show ? 1 : 0);
		}
	}
}
