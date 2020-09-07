#pragma once
#include <string>

#if defined(MOLFLOW)
//Hard-coded identifiers, update these on new release and rebuild solution
//---------------------------------------------------
static const std::string appName = "Molflow";
static const int appVersionId = 2802; //Compared with available updates. Global variable, so rebuild whole solution if changed.
static const std::string appVersionName = "2.8.2 beta";
//---------------------------------------------------
#if defined(_DEBUG)
static const std::string appTitle = "Molflow+ debug version (Compiled " __DATE__ " " __TIME__ ")";
#else
static const std::string appTitle = "Molflow+ " + appVersionName + " (" __DATE__ ")";
#endif
#endif