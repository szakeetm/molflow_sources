# Changelog

## v2.9.4 (TBA)
### Change
- ...
- ...
- Added CMake install option
### Bugfix
- ...

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