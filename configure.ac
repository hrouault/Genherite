# genherite - Evolutionary algorithm for the prediciton of genetic networks.
#    
# Copyright (C) 2005-2011 Hervé Rouault <rouault@lps.ens.fr>
# 
# This file is part of Genherite.
# 
# Genherite is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Genherite is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Imogene; see the file COPYING  If not, see
# <http://www.gnu.org/licenses/>.
#
#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.61])
AC_INIT(Genherite, m4_esyscmd([tools/git-version-gen .tarball-version]),
    [rouault@lps.ens.fr])
AM_INIT_AUTOMAKE()

AC_CONFIG_SRCDIR([src/genherite.cpp])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])
m4_include([m4/ax_blas.m4])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AM_PROG_CC_C_O

LT_INIT([dlopen])

# Checks for libraries.
AC_CHECK_LIB([m],[cos],,[AC_MSG_ERROR([math functions required.])])

##AX_BLAS( [ AC_MSG_NOTICE( [ blas libs variable : x$BLAS_LIBS BLAS_LIBS ] ) ], )
AX_BLAS(LIBS="$BLAS_LIBS $LIBS $FLIBS", AC_CHECK_LIB([gslcblas],
        [cblas_dgemm], , AC_MSG_ERROR([blas library required.])))

#AX_BLAS([],[AC_MSG_ERROR([blas library required.])])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm],,[AC_MSG_ERROR([libgsl required.])])

AC_CHECK_LIB(pthread,pthread_create)

# Checks for header files.
AC_CHECK_HEADERS([gsl/gsl_rng.h],[],[AC_MSG_ERROR([gsl headers required.])])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL

# Checks for library functions.
AC_CHECK_FUNCS([floor pow sqrt], [],
               [AC_MSG_ERROR([floor,pow and sqrt functions required])])

# Check for python
AM_PATH_PYTHON([2.5],,[AC_MSG_ERROR([no python version >=2.5 detected])])

AC_CONFIG_FILES([Makefile
                 doc/Makefile
                 src/Makefile
                 cython_wrap/Makefile])

AC_OUTPUT
