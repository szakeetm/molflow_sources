#pragma once
#include <string>

#if defined(MOLFLOW)
//Hard-coded identifiers, update these on new release and rebuild solution
//---------------------------------------------------
static const std::string appName = "Molflow";
static const int appVersionId = 0000; //Compared with available updates. Global variable, so rebuild whole solution if changed.
static const std::string appVersionName = "2.9.9 beta";

static const std::string appTitle = "Molflow+ " + appVersionName
+ " (" __DATE__ ")"
#if defined(_DEBUG)
+" [DEBUG mode]"
#endif
;
#endif