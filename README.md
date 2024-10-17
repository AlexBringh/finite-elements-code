

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