#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([calc_coupling], [1.0])

# Checks for programs.

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AC_CONFIG_SRCDIR([src/calc_coupling.f90])
AM_INIT_AUTOMAKE(subdir-objects)
AC_LANG(Fortran)
AC_PROG_FC
AC_CHECK_LIB([blas], [dgemm], [], [AC_MSG_ERROR([BLAS not found])])
AC_CHECK_LIB([lapack], [dsyev], [], [AC_MSG_ERROR([LAPACK not found])])
LDFLAGS="-static"
AC_SUBST(LDFLAGS)

AC_OUTPUT([Makefile])

