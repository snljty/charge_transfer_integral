![neko](neko.jpg)

# calc_coupling

## Calculates charge transfer integrals from files contain molecule wavefunction information.



## A brief introduction of the algorithm

Below the symbol "@" stands for matrix multiplication, and a vector is transferred into a matrix with size 1 in a dimension that suits the prerequisite of matrix multiplication.

For dimer, we have F @ C = S @ C @ E, where F is the Fock matrix , C is the coefficient matrix, S is the overlap matrix of the basis functions, and E is a diagonal matrix with energy of each orbitals on its diagonal, hence F = S @ C @ E @ C^(-1), where C^(-1) is the inverse matrix of C (in this code it is called C_inv).

Let block diagonal matrix C_sep be a matrix with the same shape as C, and the upper left block of C_sep be the coefficient matrix of monomer1, the lower right block of C_sep be the coefficient matrix of monomer2, and other two off-diagonal blocks be zero. 

For the computation details of the algorithm to obtain the coupling, just read the subroutine calc_coupling, and remember that it is done by the site-energy correction method, index1 and index2 stands for the indices of orbitals in monomer1 and monomer2, individually, and J_eff12 stands for the transfer integral between these two orbitals. The unit for all real(kind=8) numbers in this subroutine is a.u. (a.k.a. Hartree).



## Application

All files that both contains the information needed (the information of the basis functions, the coefficients of basis functions in molecule orbitals, and energies of molecule orbitals) and supported by **Multiwfn** should be supported by this program.  For example, fch/fchk files of **Gaussian** (note that traditional semi-empirical methods such as PM7 are NOT supported!) and molden files of **ORCA**. For **Gaussian** Users, key-work "NoSymmetry" must present in all of your inputs files to avoid the basis functions being rotated or translated.

If linear-dependant basis functions are found in any of these calculations, this program will NOT be available due to two reasons: some quantum software may reduce the linear-dependant basis functions so that the matrix of coefficients is a rectangle instead of a square, or even if you prevent the program from doing that, a square matrix with linear-dependant columns is singular, in both cases the matrix of coefficients is not inversible.

Phase-matching method of orbitals is currently NOT implied currently, hence ONLY the absolute value of the transfer integral calculated by this program is currently meaningful.

## Before use

The libraries **lapack** and **blas** (http://www.netlib.org/lapack/) is required to calculate C_inv from C. For windows users, you can download from http://icl.cs.utk.edu/lapack-for-windows/lapack/ or use the libraries provided in this folder. For Linux Users, you can get these libraries through yum/apt, or get the source code from http://www.netlib.org/lapack/ and compile yourself.

This executable is compiled by gfortran from TDM-gcc-9.2.0-x64 under Windows with lapack-3.9.0 compiled with the same compilers and cmake-3.15.2-x64 for Windows. If you want to compile yourself, you can change the Makefile to make sure the LIBPATH contains the path to your **lapack** and **blas** libraries, e.g., liblapack.a and libblas.a, the choices of FC (The Fortran compiler that supports at least Fortran 95 standard) and FLINKER (The Fortran Linker) are available, and others to fit the requirement of your system, then simply type "make" to generate the executable.

The program **Multiwfn** need to be properly installed before. You can obtain it from http://sobereva.com/multiwfn/ . Please set the environmental variable "Multiwfnpath" to the directory containing Multiwfn executable, and add this directory to environmental variable "PATH". Some settings in "settings.ini" for Multiwfn may also be required to  be changed before running this program.

## Usage

Use "Calc_coupling.exe --help" to see datails.

You can assign files names for dimer and monomers through command arguments, or by interactive mode. For command arguments, the argument after "--dimer" is the file for dimer (e.g. dimer.fch, dimer.molden), the argument after "--monomer1" is for the first monomer and the argument after "--monomer2" is for the second monomer.

When the matrix of coefficients is singular, you can provide a plain text file containing the Fock Matrix of the Dimer. This file should have the first line as comment, followed by the lower-triangle format Fock matrix for dimer. The "lower-triangle format" can be found in examples\dimer_Fock_Gaussian.txt. For Gaussian users, add a Link-1 task after the input file for dimer with Guess=Read IOp(5/33=3) and the same checkpoint path as above. Then search for the last "Fock matrix (alpha)" or "Fock matrix (beta)" from the Gaussian output file and save it to a new plain text file. Then rerun this program with command argument "-dF your_file_name_for_Fock_matrix_for_dimer".

All needed arguments that are not passed through command arguments, will be asked by interactive mode, please read the hint information on the screen.

## Contact

For bug report and further discussion, please send emails to snljty@sina.com .