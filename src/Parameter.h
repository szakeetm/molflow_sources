
#pragma once

#include "Distributions.h"
#include <string>
#include <vector>

class Parameter:public Distribution2D {
public:
    Parameter() {
        fromCatalog=false;
        logXinterp = false;
        logYinterp = false;
    }
	std::string name;
	bool fromCatalog=false;
};