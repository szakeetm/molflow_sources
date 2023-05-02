# Changelog

## v2.8.12

###Feature
- Immediate ZIP compression without intermediate XML

###Interface
- Support for progressbars with wide/multiline status
- Allow computer to use lock / screensaver while MolFlow is running
- Don't allow 0K temperature

### Bugfix

- Fixed convergence plotter interface elements becoming black rectangles

## v2.8.11 (2023-02-09)

### Bugfix

- Fixed crash when inserting STL geometries

### Interface

- Display collapse dialog after STL insert

## v2.8.10 (2023-01-11)

### Bugfix

- Fixed memory leak when collapsing
- Fixed error messages on dynamic desorption import

## v2.8.9 (2022-11-14)

### Feature

- Ported near-instant STL collapse and Smart Select from 2.9 beta
- Support multibody STL files
- When deleting facets, teleports are also renumbered

## v2.8.8 (2022-07-12)

### Bugfix
- Collinearity condition fixed (no more false NULL facets)
- Create facet with convex hull fixed, no more duplicate vertices