#pragma once
#include <string>

#if defined(MOLFLOW)
//Hard-coded identifiers, update these on new release and rebuild solution
//---------------------------------------------------
static const std::string appName = "Molflow";
static const int appVersionId = 2907; //Compared with available updates. Global variable, so rebuild whole solution if changed.
static const std::string appVersionName = "2.9.7 beta";

static const std::string appTitle = "Molflow+ " + appVersionName 
#if defined(_DEBUG)
+ " [DEBUG mode]"
#endif
#if defined (MEASURE_FORCES)
+ " with force measurement"
#endif
+ " (" __DATE__ ")";
#endif