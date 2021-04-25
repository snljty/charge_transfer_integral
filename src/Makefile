# Makefile for calc_coupling.exe
# libraries lapack and blas are required to be linked.
# for Windows version, they are already contained in current directory, and you can also 
# download them from http://icl.cs.utk.edu/lapack-for-windows/lapack/ yourself.
# for Linux version, you may use apt or yum to install from repository, or download source code
# from http://www.netlib.org/lapack/ and compile yourself.
# This program is not parallized, partly due to lapack and blas are not multithread libraries, 
# and the author of this program is not familiar with OpenBLAS or MKL, as well as parallel programming.

SHELL := cmd
FC := gfortran
FLINKER := $(FC)
LIBPATH := .
# LIBPATH = -L C:\\PATH\\TO\\YOUR\\LAPACK

.PHONY: all clean_tmp clean

all: calc_coupling.exe

calc_coupling.exe: calc_coupling.o
	@echo Linking $@ against $^ ...
	$(FLINKER) -o $@ $^ $(LIBPATH) -l lapack -l blas -static

calc_coupling.o: calc_coupling.f95
	@echo Compiling $@ ...
	$(FC) -o $@ -c $<

clean_tmp:
	-del /q calc_coupling.o 2> NUL

clean: clean_tmp
	-del /q calc_coupling.exe 2> NUL

