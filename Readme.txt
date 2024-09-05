This repository contains the 2D GRLMBRT code developed in the context of my PhD thesis, "General Relativistic Lattice Boltzmann Method for Radiation Transport", and my paper of the same title.
It is a standalone radiation transport code to simulate the behavior of radiation (photons and neutrinos) in the presence of a matter fluid in curved spacetime.
The code is not coupled to a matter fluid code yet, making all interactions between matter and radiation static.
The master branch is parallelized via openMP.
Tim Brohn developed the MPI branch in his Bachelor thesis, "MPI parallelization analysis of lattice Boltzmann methods for radiative transport in computational astrophysics", which allows for multi-node simulations.



How To use:
The code was developed on Linux. I have not tested it on other systems and do not guarantee it will work there.

Download with submodules:
-git clone --recurse-submodules <repository-url>

Or clone first and then init submodels:
-git clone <repository-url>
-git submodule init
-git submodule update

Dont forget to create the output dir and change the global output path in ControlFlow.hh

The following libraries are necessary to compile the code:
openmp:     -sudo apt-get install libomp-dev

To use the code, include the 'src/Radiation.h' in your .cpp project.
This will include all needed header files.
Examples of code usage can be found in exe/main.cpp.

To build the code, add your .cpp file in the 'CMakeLists.txt' as an executable,
and add the c++ flags to it (analogous to 'main.cpp'):
-add_executable(myCode.out exe/myCode.cpp ${srcs})
-target_compile_options(myCode.out PUBLIC -O3 -ffast-math)
-target_link_libraries(myCode.out OpenMP::OpenMP_CXX)

Then run cmake inside the build folder, and afterward the makefile:
-cd build
-cmake ..
-make <executable name>

The src folder contains all the source code, including the submodule eigen.
The exe folder contains all .cpp files that get compiled into executable .out files.
Here, you can find simple programs that contain tests for the individual modules of the code.
Have a look at these to understand the general code structure.
