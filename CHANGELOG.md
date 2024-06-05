# Changelog

## v2.9.24 (2024-05-16)

### Bugfix

- Fix mixed up hit colors
- Record on transparent facets only if they are before a volume event's location
  - Change in case of background collisions, radioactive decay and time-dep. simulations with cutoff
- Facet mirroring by equation: normalize coefficients
- Ship runtime DLLs on Windows (fix from 2.9.23)

### License and versioning

- Back to GPL3 license (open-source)
- Beta designation removed, 2.9 is now considered stable

## v2.9.22 beta (2024-04-03)

### Feature

- Experimental support for hard-sphere collisions with static background gas
  - See website blog post for more info 

### Bugfix

- Don't update manual texture limit textboxes on every second
- Fixed Facet Coordinates window crash when inserting vertex with no facet sel.
- Crash (mostly on macOS) when issuing an "unselect all" command

### Fedora Linux

- Removed outdated shared libraries, MolFlow+ is now built on AlmaLinux 9

### License

- Removed GNU Scientific Library dependencies, moving to closed source
  - Discussion in progress about the long term. Probably the physics (molflowCLI)
  - will remain open-source, and the GUI closed-source

## v2.9.21 beta (2024-02-19)

### Feature

- molflowCLI: allow 'specific outgassing' parameter change

### Interface

- Histogram settings panel shows memory estimate of histograms
- Show actual desorption limit in Global Settings on change button
- Sanitize desorption limit user input before applying
- molflowCLI: feedback (success/error) for parameter changes
- Easier to understand Texture Scaling window controls
- Select by texture type: allow 0 as min. ratio

### Bugfix

- MolFlow runs on OpenGL1.1 environments again (i.e. Remote Desktop servers)
- Selection sometimes not finding facets on click fixed
- Inserting geometries: formula and selection group offsetting working again
- Stop simulation correctly if desorption limit reached almost instantly after start
- molflowCLI: allow to run simulations for very small desorption limits
- molflowCLI: --setParamsByFile with outgassing: fixed factor of 10 error
- crash after Copy Mirror command
- every second line empty on Texture Plotter export

## v2.9.18 beta (2024-01-08)

- Hotfix after 2.9.17: formula parsing works again

## v2.9.17 beta (2023-12-14)

### Feature

- Cross section viewer (find in View menu)
- Last camera view restored when loading a file (if saved with 2.9.17+)
- Texture scaling settings (auto/manual limits) saved and loaded again
- Security: hard-code updater URL to avoid tampering with updater_config.xml
- Security (macOS, Linux): disable custom post-update scripts,
  restrict to give executable permission by code
- Particle logger: all processes add to log (previously only first)

### Interface

- Formula editor: allow to batch copy values for all moments
- Auto-update formulas moved from Global Settings to Formula Editor
- Easier to understand molflowCLI status messages and help
- Allow "-inf" and "+inf" texture limits if no textures (or no moments)

### Bugfix

- Histogram plotter had no values for global histograms and for facet-based flight time
- Angle map CSV import now allows 64-bit unsigned integers (previously max 2^31-1)
- Full screen mode works again
- Moments load progressbar off-by-one error
- STL loading: better handle scientific numbers, give meaningful error messages
- Clear formulas on File/New geometry

## v2.9.16 beta (2023-09-18)

### Interface

- Message boxes (error and confirmation) now accept Enter and Esc
- Global settings right pane now shows that settings are file-based

### Bugfix

- Desorption limit working again (GUI & CLI)
- No more false warnings for app updater
- Low-flux mode and gas mass saved with file, not with app settings

## v2.9.15 beta (2023-08-21)

### Feature

- Facet temperatures can be time-dependent
    - Currently only for non-outgassing facets
- Lowered XML file memory usage on load by ~25%
- XML format change: refer to parameters by name, not by ID
    - Backwards compatible with previous versions, except for temperature

### Interface

- Time-dependent sticking, opacity and temp. now shown in Facet details
    - In blue color, value for current moment, as defined in Time Settings
- Total outgassing now calculated by simulation model (as opposed to GUI)
    - Needs manual recalc. button push if model not synced
- Model sanity check now also done by simulation model

### Bugfix

- Fixed hang when exiting Molflow with simulation running

## v2.9.14 beta (2023-08-08)

### Feature

- Load status window for simulation-interface operations that take >0.5s
- Possibility to abort simulation start/stop/reload if taking too long

### Interface

- Resizable Global Settings window (for long thread list)
- Approximation for Sim.Manager and thread memory usage

### Bugfix

- Formula parsing fixed for consecutive multiplications/divisions
- Convergence data saved again
- Prevent crash if convergence data is expected but loaded file hasn't

## v2.9.13 beta (2023-07-21)

### Feature

- Support binary STL files (single body, load only)
- Formulas can refer to formulas defined above
- Progress printed for file saving (GUI and CLI)
- Google Analytics replaced with Matomo (GDPR)
- CLI: --notInteractive flag to suppress percentage update (log files)

### Bugfix

- Include missing vcomp140.dll for OpenMP parallelization on Windows
- Fixed incorrect filename check saying xml not supported

## v2.9.12 beta (2023-06-01)

### Bugfix

- Hotfix for bug causing updater to crash. Update manually if updater fails.

## v2.9.11 beta (2023-05-23)

### Feature

- Faster geometry initialization (after load, move, scale etc.)

### Bugfix

- Restored loading of textures
- Prevented some race conditions

## v2.9.10 beta (2023-05-15)

### Feature

- Parallelized many slow parts of the code:
    - Facet meshing
    - Facet and vertex search on screen (when selecting)
    - OpenGL geometry building (geometry initialization)
    - Temporary counters creation for subprocesses (load, reset)
- Speedup: only reload simulation to threads if it changed (not on resume/reset)
- Auto-updater can now execute post-install scripts
   - Gives executable permission to updates (Mac and Linux)
- Removed: legacy .ASE file loader

### Interface

- New unified progressbars for GUI and CLI (ready for file load, save coming soon)
- Don't accept 0 Kelvin temperature
- Deduplicate shared edges of selected facets (faster rendering)
- Show progress message when exiting MolFlow to identify blocks

### Bugfix

- Fixed all facets rendering once per structure
- Fixed Global Settings impossible to minimize

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