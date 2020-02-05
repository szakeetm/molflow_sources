//
// Created by pbahr on 03/02/2020.
// With several additions from Ingo Wald's GDT helper files.
//

#ifndef MOLFLOW_PROJ_HELPER_OUTPUT_H
#define MOLFLOW_PROJ_HELPER_OUTPUT_H

#define MF_TERMINAL_RED "\033[1;31m"
#define MF_TERMINAL_GREEN "\033[1;32m"
#define MF_TERMINAL_YELLOW "\033[1;33m"
#define MF_TERMINAL_BLUE "\033[1;34m"
#define MF_TERMINAL_RESET "\033[0m"
#define MF_TERMINAL_DEFAULT MF_TERMINAL_RESET
#define MF_TERMINAL_BOLD "\033[1;1m"

#ifndef PRINT
# define PRINT(var) std::cout << #var << "=" << var << std::endl;
# define PING std::cout << __FILE__ << "::" << __LINE__ << ": " << __FUNCTION__ << std::endl;
#endif

#endif //MOLFLOW_PROJ_HELPER_OUTPUT_H
