#include <stdlib.h>
#include <string>
#include <stdio.h>
#include <iostream>
#include <Windows.h>
#include <io.h>
#include <direct.h>

std::string exec(char* cmd);
int FileExists(char *fileName);

int main(int argc,char* argv[]) {
	//__debugbreak();
	char key;
	std::string result;
	if (argc!=2) {
		std::cout<<"No input file argument.";
		ShowWindow( GetConsoleWindow(), SW_RESTORE );
		std::cin>>key;
		return 0;
	}
	char command[4096];
	char fileName[1024];
	char fileNameWith7z[1024];
	char fileNameGeometry[1024];
	memcpy(fileName,argv[1],strlen(argv[1])*sizeof(char));
	fileName[strlen(argv[1])*sizeof(char)]='\0';
	std::cout<<"\nargv0: "<<argv[0];
	std::cout<<"\nargv1: "<<argv[1];
	sprintf_s(fileNameWith7z,"%s7z",fileName);
	if (!FileExists("7za.exe")) {
		printf("\n7za.exe not found. Cannot compress.\n");
			std::cin>>key;
			return 0;
	}
	char *dir;
	dir = strrchr(fileName,'\\');
	memcpy(fileNameGeometry,fileName,sizeof(char)*(dir-fileName));
	fileNameGeometry[dir-fileName]=NULL;
	sprintf_s(fileNameGeometry,"%s\\Geometry.geo",fileNameGeometry);
	sprintf_s(command,"move \"%s\" \"%s\"",fileName,fileNameGeometry);
	result=exec(command);
	char CWD [MAX_PATH];
	_getcwd( CWD, MAX_PATH );
	sprintf_s(command,"cmd /C \"pushd \"%s\"&&7za.exe u -t7z \"%s\" \"%s\"&&popd\"",CWD,fileNameWith7z,fileNameGeometry);
	std::cout<<"\nCommand: "<<command<<"\n\nStarting compression...\nYou can continue using Molflow while compressing.\n";
	result=exec(command);
	size_t found;
	//printf("\nresult: %s\n",result);
	found=result.find("Everything is Ok");
	if (found!=std::string::npos) {
		printf("\nCompression seems legit. Deleting GEO file.");
		remove(fileNameGeometry);
		return 0;
	}
	//printf("\nresult: %s\n",result);
	ShowWindow( GetConsoleWindow(), SW_RESTORE );
	sprintf_s(command,"move \"%s\" \"%s\"",fileNameGeometry,fileName);
	result=exec(command);
	printf("\nSomething went wrong during the compression, read above. GEO file kept."
		"\nType any letter and press Enter to exit\n");
	std::cin>>key;
	return 0;
}

std::string exec(char* cmd) {
    FILE* pipe = _popen(cmd, "r");
    if (!pipe) return "ERROR";
    char buffer[128];
    std::string result = "";
    while(!feof(pipe)) {
        if(fgets(buffer, 128, pipe) != NULL)
                result += buffer;
		        printf(buffer);
    }
	result=result+'0';
    _pclose(pipe);
    return result;
}

int FileExists(char *fileName) {

  struct _finddata_t seqfile;
  intptr_t h;

  if( (h=_findfirst( fileName , &seqfile )) != -1L ) {
	  _findclose(h);
	}

  return (h != -1L);
}