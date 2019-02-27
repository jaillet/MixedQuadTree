cd build
rm -R *
cmake ../src -D_CMAKE_DEBUG=Debug
cmake .
make