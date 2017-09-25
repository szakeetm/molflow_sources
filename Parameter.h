#ifndef _PARAMETERH_
#define _PARAMETERH_

#include "Distributions.h"
#include <string>
#include <vector>

class Parameter:public Distribution2D {
public:
	std::string name;
	bool fromCatalog;
	
	Parameter();
};

/*
class StringClass:public Distribution2D {
public:
	std::string name;
	bool fromCatalog;
	StringClass();
};*/
#endif