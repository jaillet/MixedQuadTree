mkdir build
cd build
rm -R *
cmake ../src -D_CMAKE_BUILD_TYPE=Debug
cmake .
make