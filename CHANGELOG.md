# Changelog

## v2.9.10 beta

### Feature

- Auto-updater can now execute post-install scripts
   - Gives executable permission to updates (Mac and Linux)

### Interface

- Don't allow 0K temperature

## v2.9.9 beta (2023-04-13)

### Feature
- Immediate ZIP compression without intermediate XML
   - This solves hang when saving in cloud-synced folders
- New Clipper2 polygon library providing drastically faster texture meshing
- Facet collapse can leave max. vertices per facet (default: 100)
- PugiXML "compact" mode reduces memory consumption when opening/saving large XML files

### Interface
- Fully rewritten OpenGL rendering, much faster graphics allowing up to 1M facets
   - Using triangles for volume and direction textures
   - Using vertex arrays for facet wireframe, including deduplication of shared edges
   - Using gridlines (instead of cell quads) for texture grid
- Texture mesh (now renamed "grid") rendering drastically faster
- Coordinate axes in lower left corner
- File save progress displayed in console
- Color-highlight facets only if profile plotter is open
- Support for progressbars with wide/multiline status
- Added few missing options to Global Settings
- Added IMGui test window (Test/ImGui test suite/Demo window in menu)
- Light-rotation (ALT + right-drag) speed independent of "Angle step" setting

### Bugfix

- Fixed a few OpenGL crashes related to camera views
- Fixed crash in certain cases when collapsing
- Fixed convergence plotter interface elements becoming black rectangles
- Fixed crash at "auto-create difference of 2 facets" function
- Fixed crash at 0-depth geometries at certain camera angles
- "Front" and "Back" facet render modes were inversed
- Don't ask twice for overwriting files
- MPI version compiles again
- MPI builds throw error if no MPI is present (previously silent fail)

## v2.9.8 (2022.12.05)
### Feature
- Added option to measure molecular forces and torque (See blog post on website)

### Bugfix
- Fixed very serious memory leak when collapsing a geometry
- Fixed 100% CPU usage when simulation paused
- Fixed app updater checking for debian versions on fedora

## v2.9.7 (2022.10.13)
### Bugfix
- Simulation is properly reset on facet change even if scene auto-update is off
- Fixed crash on texture change while simulation was running
- Fedora builds: updater no longer downloads Debian updates
- Proper reset allows to save geometry after isolated vertices are removed
- Correct angle map status refresh when simulation starts

## v2.9.6 (2022.08.17)
### Bugfix
- Collinearity condition fixed (no more false NULL facets)
- Create facet with convex hull fixed, no more duplicate vertices
- Load and save old TXT format including facet hits
- Fix crash when file loading is unsuccessful

## v2.9.5 (2022.06.23)

### Feature
- Experimental windows based on IMGui
- Formulas are case insensitive
- Allow to import SYN version 12 (Upcoming Synrad version)
- Null facet detection criteria matched to Synrad

### Bugfix
- Rare crash when collapsing facets or vertices
- Linux+Mac: Fixed install location with CMake
- Correct saving of part of the geometry (File->Save selected facets)
- Inserting STL files doesn't crash Molflow anymore
- CLI: crash when no input file is given
- CLI: no more creation if intermediate result files
- CLI: error handling for file copy

## v2.9.4 (2022.04.14)
### Feature
- MPI simulations (see documentation)

### Change
- Dramatically accelerated collapse operations (vertices/facets)
- Load files via drag and drop from file browser
- Added CMake install option
- CLI option for an auto generated test case (pipe) based on oblique prisms
- Facet highlighting is disabled when more than 500 facets are selected to increase rendering performance
- Time-dependent mode: particle reflection takes into account sojourn time

### Bugfix
- Crash on creating a new facet by difference
- Various synchronization bugs between recording and generating angle map states
- Not copying last point before 'latestMoment' of user-defined outgassing resulted in wrong TD results
- Convergence data would not be plotted when loading from a XML file
- Arrow heads of coordinate system are properly rendered
- CLI could endlessly run when neither a valid time or desorption limit has been set
- Thread status shows properly on all supported OS
- CLI will properly save changes done via parameter changes to the output file
- Fixed incorrect generation of angle maps

## v2.9.3 (2022.01.31)
### Change
- Changing parameters (CLI) will now give errors or warnings for edge cases
- Profile plotter now gives a unified warning for each plot concerned
- Particle Logger logs for all threads/particles again
### Bugfix
- Fixed CSV export for CLI misfunctioning, due to a bug with fmtlib
- Multiple profiles wouldn't plot correctly depending on display mode
- Desorptions would not get recorded for particle logger
- Fixed a rare case of heap corruptions when loading XML files with multiple structures

## v2.9.2 (2021.11.22)
### Feature
- User notification on repetitive update fails
### Bugfix
- Fixed empty file filters in Import Desorption file dialog on Windows
- Fixed Molflow preventing the machine going into standby mode
- Fixed a rare crash concerning files using the default parameter catalog
- Fixed occasions where exponential values would be truncated in the Profile Plotter
- Fixed desorptions not getting handled properly for line visualization
- Fixed a crash after successively loading new geometries
- Removed thread information from Global Settings where no details can be given
### Development
- Changes surrounding the App Updater

## v2.9.1 (2021.11.08)
### Feature
- Added two new export commands for the CLI application
    - `--writeFacetDetails` and `--writeFacetQuantities` (see Documentation)
- Smart Selection now has almost no wait time
### Bugfix
- Fixed simulations on Windows sometimes feeling unresponsive
- Texture values for small/cut-off texture cells are now properly scaled
- Fixed various application crashes
- Fixed cases where the wrong output file format was used for CLI save files
- Fixed dynamic desorptions leading to faulty outgassing maps when imported from Synrad
- Fixed cases where the simulation time would increment by faulty values

## v2.9.0 (2021.09.30)
### Feature
- Added CLI application (see [documentation](https://molflow.web.cern.ch/node/381))
- Improved simulation performance with major rework of the simulation framework
- Support dynamic outgassing from non-square facets
- Added MacOS ARM support for M1 chips
- Added support for multibody STL files
### Minor changes
- Some usability changes for several tool windows
- Removed texture limits from output file
- Renamed "overwrite" to "save/overwrite" for "View" and "Selection" menus
### Developers
- Added google-test support
### Bugfix
- Fixed proper distribution of desorptions on all CPU cores when a limit is set
- Log/Log interpolation used when importing desorption yield for dynamic desorption
- Fixed double to double comparison in Facetdetails
- Keep child windows on screen when size of main window is reduced
- Fix crash when bounding box of the geometry had 0-length side
- Fixed some bugs when files contained "." characters or had no extension
- Fixed some issues concerning import and export from STL files