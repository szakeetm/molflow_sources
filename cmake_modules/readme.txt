Molflow is separated into three standalone CMake projects:

- molflow_core: static library with shared gui/cli code
- molflowCLI (molflow_cli folder): CLI, depends on molflow_core
- molflow (molflow_gui folder): GUI, depends on molflow_core and molflow_cli