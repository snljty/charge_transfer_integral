# Makefile for calc_coupling_gau_semi_empirical_use_out

SHELL = cmd
FC := gfortran
FLINKER := $(FC)
INSTALL_PATH := ../bin

.PHONY: all
all: calc_coupling_gau_semi_empirical_use_out.exe

.PHONY: install
install: all
	dir /ad /b "$(INSTALL_PATH)" 1> NUL 2> NUL || mkdir "$(INSTALL_PATH)"
	copy /y calc_coupling_gau_semi_empirical_use_out.exe "$(INSTALL_PATH)"

calc_coupling_gau_semi_empirical_use_out.exe: calc_coupling_gau_semi_empirical_use_out.obj
	@echo Linking $@ against $^ ...
	$(FLINKER) -o $@ $^ -static -s

calc_coupling_gau_semi_empirical_use_out.obj: calc_coupling_gau_semi_empirical_use_out.f90
	@echo Compiling $@ ...
	$(FC) -o $@ -c $< -O2 -s

.PHONY: clean
clean:
	-del /q calc_coupling_gau_semi_empirical_use_out.obj 2> NUL

.PHONY: veryclean
veryclean: clean
	-del /q calc_coupling_gau_semi_empirical_use_out.exe 2> NUL

