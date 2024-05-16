#pragma once
#include <string>

#if defined(MOLFLOW)
//Hard-coded identifiers, update these on new release
//---------------------------------------------------
static const std::string appName = "MolFlow";
static const int appVersionId = 2924; //Compared with available updates 
static const std::string appVersionName = "2.9.24";

static const std::string appTitle = "MolFlow " + appVersionName
+ " (" __DATE__ ")"
#if defined(DEBUG)
+" [DEBUG mode]"
#endif
;
#endif