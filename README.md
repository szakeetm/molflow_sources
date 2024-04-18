# MolFlow+

A Monte Carlo simulator for Ultra High Vacuum systems

**Authors:** Roberto KERSEVAN, Marton ADY, Tymoteusz MROCZKOWSKI, Pascal Rene BAEHR, Jean-Luc PONS  
**Website:** https://cern.ch/molflow  
**Copyright:** CERN (2024)  
**License:** GNU GPL v3. [Details](https://molflow.docs.cern.ch/about/#license)

<img src="https://molflow.web.cern.ch/sites/molflow.web.cern.ch/files/pictures/2018-10-09%2016_14_20-PowerPoint%20Slide%20Show%20%20-%20%20Presentation1.png" alt="Molflow image" width="800"/>

# Cloning the project (all OS)

* Clone the Molflow project with `git clone`
* Go into the `molflow` directory and initialize the `src_shared` submodule
* Update the `src_shared` submodule

```
git clone https://gitlab.cern.ch/molflow_synrad/molflow.git
cd molflow
git submodule init
git submodule update
```
  
# Building

Molflow uses `cmake` for its build system. On Windows it comes with *Visual Studio 2022* or it has to be built/downloaded manually from [cmake's download page](https://cmake.org/download/).
On Linux and macOS `cmake` can be installed with package managers.

## Windows

In Windows, you need to clone the project (see below) and open in Visual Studio Community 2022.

With Visual Studio's CMake support, you can now open the main folder containing *CMakeLists.txt* (File->Open Folder) and build the project directly (Build / Build All in menu or F7)

## Prepare for building - Debian Linux (like Ubuntu)

### Ubuntu 20

```
sudo apt install build-essential git cmake gcc g++
sudo apt install libsdl2-dev libpng-dev libgtk-3-dev libcurl4-gnutls-dev libatlas-base-dev p7zip
```
## Prepare for building - Fedora Linux (like CentOS)

### CentOS 9 stream

Set up build tools:

```
dnf install -y git cmake gcc g++
```

Install required packages to compile molflow:

```
dnf install -y zlib-devel libpng-devel gtk3-devel libcurl-devel SDL2-devel p7zip
```
Some packages might require enabling the `crb` repo:
```
dnf config-manager --set-enabled crb
/usr/bin/crb enable
dnf install epel-release
```
In our test setup (Windows Subsystem for Linux), we used X-Win32 window server, which required setting the display and installing mesa drivers:
```
export DISPLAY=0:0
dnf install mesa-dri-drivers
```

### CentOS 8

```
dnf config-manager --set-enabled PowerTools
dnf in -y https://dl.fedoraproject.org/pub/epel/epel-release-latest-8.noarch.rpm
dnf in -y libglvnd-opengl libpng15 SDL2 p7zip
dnf group install "Development Tools"
```

### CentOS 7

```
sudo yum install epel-release
sudo yum install p7zip SDL2 libglvnd-opengl
```

If your CMake is outdated, install a newer version:
```
wget https://github.com/Kitware/CMake/releases/download/v3.18.3/cmake-3.18.3.tar.gz
tar -zxvf cmake-3.18.3.tar.gz
cd cmake-3.18.3
./bootstrap
make
sudo make install
```

## Prepare for building - macOS

* Use Homebrew to install build tools, like g++-8, the SDL2 library, libpng, curl, p7zip  

The procedure looks as follows:
1. Install command line tools
  - `xcode-select --install`
2. Install homebrew
  - from https://brew.sh/
  - `/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"`
3. Install cmake and the necessary dependencies
  - `brew install cmake libpng sdl2 p7zip libomp libcurl`

## Manual build with CMake (Linux/MacOS)

* Go to the project root folder `cd molflow`
* Create a build directory `mkdir build`
* Enter the build directory `cd build`
* Init the cmake project inside the new directory `cmake ..`
* Start the build with `make`
```
cd molflow
mkdir build
cd build
cmake ..
make -j8
```

## CMake configuration flags

* `CMAKE_BUILD_TYPE`
  * by default `Release`
  * can set to `Debug`
  * these two options are case sensitive!
* `USE_TESTS`
  * by default `OFF`
  * set to `ON` to build `testsuite.exe` that runs molflowCLI unit tests
* `USE_MPI`
  * by default `OFF`
  * use `ON`, most likely with interface disabled (below) for MPI builds
* `NO_INTERFACE`
  * by default `OFF`
  * use `ON` to build only `molflowCLI` and not the GUI `molflow` (useful for clusters, headless servers, MPI builds)
* `CMAKE_C_COMPILER` and `CMAKE_CXX_COMPILER`
  * Allows to use not the default compiler. For example, on macOS, by default, Apple CLang is used, but you could use Homebrew gcc:
  ```
  cmake .. -DCMAKE_C_COMPILER=/usr/local/bin/gcc-13 -DCMAKE_CXX_COMPILER=/usr/local/bin/g++-13
  ```

Usage example (for an MPI build without the GUI):

```
cmake .. -DUSE_MPI=ON -DNO_INTERFACE=ON
```

## CMake / make installation

Given a default Cmake build procedure `mkdir build && cd build && cmake ..` Molflow can be installed into the users home folder by `cmake --install .` --> ~/molflow/

To change the directory, the installation prefix can be adjusted `cmake --install . --prefix ~/Apps/` --> ~/Apps/molflow

For a make installation, after the standard CMake procedure, by default Molflow will also be installed in the users home folder `make install` --> ~/molflow/

The installation path can be changed by adding an installation prefix to the CMake build command `cmake -DCMAKE_INSTALL_PREFIX:PATH=~/Apps/ ..` and `make install` --> ~/Apps/molflow

## Build example on Mac

1) We install Homebrew and dependencies:

```
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
brew install cmake libpng sdl2 p7zip libomp
```

2) then we clone and build the repo:

```
git checkout https://gitlab.cern.ch/molflow_synrad/molflow.git
cd molflow
git submodule init
git submodule update
mkdir build
cd build
cmake ..
make -j8
```

On Linux, the dependency part is different (using `apt` or `yum`), but the second part is the same.

## Branches

* For MolFlow, `master` is the 2.9 beta, and `master_2.8` is the 2.8 public
* For SynRad, `master_cli` is the 1.5 beta, and `master` is the 1.4 public

## Building for MPI

- Prepare the environment with MPI and a compatible GCC version (GCC_VERSION >= 8)
- Make sure a recent version of CMake is available (>= 3.12)
    - Otherwise setup cmake from `https://github.com/Kitware/CMake/releases/`
    - `wget https://github.com/Kitware/CMake/releases/download/v3.20.1/cmake-3.20.1-linux-x86_64.tar.gz` 
    - `tar -xvf cmake-3.20.1-linux-x86_64.tar.gz.1`
    - `cd cmake-3.20.1-linux-x86_64` 
- Install Molflow by cloning the Repository
    - `git clone https://gitlab.cern.ch/molflow_synrad/molflow.git`
    - `cd molflow`
- Build molflow with the installed cmake (or the pre-installed cmake)
    - we require headless (NO_INTERFACE) and MPI (USE_MPI) for the CMake build
    - `mkdir build && cd build` 
    - `cmake -DCMAKE_C_COMPILER=/path_to_custom_gcc -DCMAKE_CXX_COMPILER=/path_to_custom_g++ -DNO_INTERFACE=ON -DUSE_MPI=ON ..`
- Use as explained by the MPI service of choice, e.g. with `mpirun`
    `mpirun -n 64 ./molflowCLI -f TestCases/06-dynamic_desorption_from_synrad.xml -t 180 --reset`

## Build using vcpkg

Detailed instructions coming soon, after testing.

Working as of 2024.03.07:

- `git clone https://github.com/microsoft/vcpkg.git`
- `cd vcpkg`
- `git checkout 2024.01.12` - This vcpkg version, with its packages, works on all platforms (avoids SDL 2.30+ causing failed start on Windows Remote Desktop)
- `./bootstrap-vcpkg.sh` or `./bootstrap-vcpkg.bat`
- `./vcpkg integrate install` - note the toolchain location in the command's output
- On all platforms: `./vcpkg install cereal cimg curl fmt libpng zlib pugixml sdl2`
- On Fedora (omit SDL2, as it has a problem): `./vcpkg install cereal cimg curl fmt libpng zlib pugixml`
  - Install `SDL2-devel` with yum or dnf instead
- Navigate to molflow repo, and from a `build` or similar dir:
- `cmake .. "-DCMAKE_TOOLCHAIN_FILE=/path/to/vcpkg/scripts/buildsystems/vcpkg.cmake"` (toolchain location as you noted down)
- `make`

# Running

## Windows

Run `molflow.exe` (in the `bin/release` directory if you built it yourself)

## Linux

* Make `molflow`, `molflowCLI` and `compress` executable:  
  `chmod +x molflow molflowCLI compress`  
  (No need if you built it from source)
* Run `./molflow`  

Detailed instructions: 

* [For Debian](https://molflow.web.cern.ch/node/296)
* [For Fedora](https://molflow.web.cern.ch/node/302)

## macOS

* Use Homebrew to install dependencies, like `brew install sdl2 libpng gsl gcc p7zip libomp`
* Make `molflow`, `molflowCLI` and `compress` executable  
  `chmod +x molflow molflowCLI compress`  
  (No need if you built it from source)
* Run `./molflow`

[Detailed instructions (macOS)](https://molflow.web.cern.ch/node/294)

# Repository snapshots
Commits are constantly pushed to this primary repo, and some of them might break - temporarily - the build scripts. If you want to fork Molflow, it is recommended that you download a [tag snapshot](https://gitlab.cern.ch/molflow_synrad/molflow/-/tags) of a guaranteed-to-work state. Usually these snapshots are made at every public release of Molflow.
