#ifndef __vtkLHPSetGet_h_
#define __vtkLHPSetGet_h_

#include "vtkSetGet.h"

//
// Set built-in type.  Creates member Set"name"() (e.g., SetVisibility());
// flagsvar is a member attribute with the flags, flagvalue is the value to be set / reset
#define vtkSetFlagsBodyMacro(name, value, flagsvar, flag) \
{ \
	vtkDebugMacro(<< this->GetClassName() << " (" << this << "): setting " #name " to " << value); \
	value *= flag; \
	if ((flagsvar & flag) != value) \
	{\
		flagsvar = (flagsvar & ~flag) | value; \
		this->Modified(); \
	} \
}

//
// Set built-in type.  Creates member Set"name"() (e.g., SetVisibility());
// flagsvar is a member attribute with the flags, flagvalue is the value to be set / reset
#define vtkSetFlagsMacro(name, flagsvar, flag) \
virtual void Set##name (int _arg) \
	{ \
		vtkSetFlagsBodyMacro(name, _arg, flagsvar, flag) \
	} \


//
// Get built-in type.  Creates member Get"name"() (e.g., GetVisibility());
//
#define vtkGetFlagsMacro(name,flagsvar,flag) \
virtual int Get##name () { \
  vtkDebugMacro(<< this->GetClassName() << " (" << this << "): returning " << #name " of " << this->name ); \
  return (flagsvar & flag) / flag; \
}


#endif // __vtkLHPSetGet_h_