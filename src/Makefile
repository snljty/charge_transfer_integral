# Makefile for calc_coupling.exe
# libraries lapack and blas are required to be linked.
# for Windows version, they are already contained in current directory, and you can also 
# download them from http://icl.cs.utk.edu/lapack-for-windows/lapack/ yourself.
# for Linux version, you may use apt or yum to install from repository, or download source code
# from http://www.netlib.org/lapack/ and compile yourself.
# This program is not parallized, partly due to lapack and blas are not multithread libraries, 
# and the author of this program is not familiar with OpenBLAS or MKL, as well as parallel programming.

SHELL = cmd
FC := gfortran
FLINKER := $(FC)
INSTALL_PATH := ../bin
LAPACKPATH := C:/lapack-3.12.0
LIBPATH := -L $(LAPACKPATH)/lib

.PHONY: all
all: calc_coupling.exe

.PHONY: install
install: all
	dir /ad /b "$(INSTALL_PATH)" 1> NUL 2> NUL || mkdir "$(INSTALL_PATH)"
	copy /y calc_coupling.exe "$(INSTALL_PATH)"

calc_coupling.exe: calc_coupling.obj
	@echo Linking $@ against $^ ...
	$(FLINKER) -o $@ $^ $(LIBPATH) -l lapack -l blas -static -s

%.obj: %.f90
	@echo Compiling $@ ...
	$(FC) -o $@ -c $< -O2 -s -ffpe-summary=none

.PHONY: clean
clean:
	-del /q calc_coupling.obj 2> NUL

.PHONY: veryclean
veryclean: clean
	-del /q calc_coupling.exe 2> NUL

