# MolFlow+
A Monte Carlo simulator for Ultra High Vacuum systems

**Authors:** Marton ADY, Pascal Rene BAEHR, Roberto KERSEVAN, Jean-Luc PONS  
**Website:** https://cern.ch/molflow  
**Copyright:** CERN (2021)  
**License:** GNU GPLv2 or later

<img src="https://molflow.web.cern.ch/sites/molflow.web.cern.ch/files/pictures/2018-10-09%2016_14_20-PowerPoint%20Slide%20Show%20%20-%20%20Presentation1.png" alt="Molflow image" width="800"/>


  
# Building
Molflow uses cmake for its build system. On Windows it comes with *Visual Studio 2019* or it has to be build/downloaded manually from [cmake's download page](https://cmake.org/download/).
On Linux and macOS it is part of most package managers.

## Windows

In Windows, you need to clone the project (see below) and open in Visual Studio Community 2019.

With *Visual Studio 2019*'s CMake support, you can now open the main folder containing *CMakeLists.txt* (File->Open Folder) and build the project directly

## Prepare for building - Debian Linux, like Ubuntu

- Install the necessary packages with the apt package manager, for Ubuntu 20:
```
sudo apt install p7zip-full
sudo apt install build-essential git cmake
sudo apt install libsdl2-dev libpng-dev libgtk-3-dev libgsl-dev libcurl4-gnutls-dev gsl-bin libatlas-base-dev p7zip
```
## Prepare for building - Fedora Linux, like Cent OS

- Install the necessary packages with the package manager

yum (e.g. CentOS 7):
```
sudo yum install epel-release
sudo yum install p7zip SDL2 gsl libglvnd-opengl
```

dnf (e.g. CentOS 8):
```
dnf config-manager --set-enabled PowerTools
dnf in -y https://dl.fedoraproject.org/pub/epel/epel-release-latest-8.noarch.rpm
dnf in -y libglvnd-opengl libpng15 SDL2 p7zip
dnf group install "Development Tools"
```
You might have to create a symlink for the GSL library:
```
ln -s /usr/lib64/libgsl.so.23 /usr/lib64/libgsl.so.0
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

* Use Homebrew to install build tools, like g++-8, the SDL2 library, libpng, gsl, curl, p7zip  

The procedure looks as follows:
1. Install command line tools
  - `xcode-select --install`
2. Install homebrew
  - from https://brew.sh/
  - `/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"`
3. Install cmake and the necessary dependencies
  - `brew install cmake libpng gsl sdl2 p7zip libomp`

## Cloning the project (all OS)
* Clone the Molflow project with `git clone`
* Go into the `Molflow` directory and:
  * `git submodule init`
  * `git submodule update`
```
git clone https://gitlab.cern.ch/molflow_synrad/molflow.git
cd molflow
git submodule init
git submodule update
```

## Manual build with CMake (Linux/MacOS)

* Go back to the project root folder `cd ..`
* Create a build directory
* Init the cmake project inside the new directory
* Start the build with `make`
```
cd ..
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=RELEASE -DUSE_TESTS=OFF ..
make
```

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
    - `~/cmake-3.20.1-linux-x86_64/bin/cmake -DCMAKE_C_COMPILER=/path_to_custom_gcc -DCMAKE_CXX_COMPILER=/path_to_custom_g++ -DNO_INTERFACE=ON -DUSE_MPI=ON ..`
- Use as explained by the MPI service of choice, e.g. with `mpirun`
    `mpirun -n 64 ./molflowCLI -f TestCases/06-dynamic_desorption_from_synrad.xml -t 180 --reset`

# Running
## Windows
Use the shortcut (that changes the working directory and launches *molflow.exe*) in *bin\win\release*
## Linux (Debian)
* Install dependencies with the *apt* package manager, like *libsdl2-2.0*, *gsl-bin*, *libatlas-base-dev*  
* In the *release/bin* folder, make *molflow* and *compress* executable
* Run *molflow*  

[Detailed instructions here](https://molflow.web.cern.ch/node/296)
## Linux (Fedora-based, like Cent OS)
* Make *launch_molflow.sh* executable
* Run *launch_molflow.sh* (It adds the lib folder to the library search path and launches molflow)

[Detailed instructions here](https://molflow.web.cern.ch/node/302)
## macOS
* Use Homebrew to install dependencies, like *sdl2*, *libpng*, *gsl*, *gcc*  
* In the *release/bin* folder, make *molflow* and *compress* executable
* Run *molflow*  

[Detailed instructions here](https://molflow.web.cern.ch/node/294)

# Repository snapshots
Commits are constantly pushed to this primary repo, and some of them might break - temporarily - the build scripts. If you want to fork Molflow, it is recommended that you download a [snapshot](https://molflow.web.cern.ch/content/developers) of a guaranteed-to-work state. Usually these snapshots are made at every public release of Molflow.
