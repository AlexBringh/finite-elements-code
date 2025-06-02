Elasto-Plastic Finite Element code developed for Master Thesis at HVL "Computational Plasticity with Work Hardening and the Flow Rule using the Finite Element Method", June 2nd, 2025 by Alexander B. Ringheim

Code is not tested against real results and is not to be used outside of evaluation for master thesis untill further notice. The repository is open so that the censors may evaluate the work.

Building with CMake:


Configuration of build:

Build with all modules
cmake -DBUILD_ALL=ON ..
cmake --build

Build with specific module
cmake -DBUILD_ALL=OFF -DBUILD_WITH_<MODULE-NAME>=ON ..


Building the project after configurations are set:
cd build/
cmake --build