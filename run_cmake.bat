@echo off

rem set CMAKE_PREFIX_PATH=C:/lapack-3.12.0
set CMAKE_PREFIX_PATH=C:/openblas-0.3.30
mkdir build
cd build
cmake .. -G "MinGW Makefiles" -D CMAKE_INSTALL_PREFIX=C:/calc_coupling -LH
rem cmake --build . -j
rem cmake --install .
cd ..
