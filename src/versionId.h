#pragma once
#include <string>

#if defined(MOLFLOW)
//Hard-coded identifiers, update these on new release
//---------------------------------------------------
static const std::string appName = "Molflow";
static const int appVersionId = 2923; //Compared with available updates 
static const std::string appVersionName = "2.9.23 beta";

static const std::string appTitle = "Molflow+ " + appVersionName
+ " (" __DATE__ ")"
#if defined(DEBUG)
+" [DEBUG mode]"
#endif
;
#endif