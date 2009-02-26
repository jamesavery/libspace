# This file is part of The Generalized Sturmian Library
#
# Authors:  James Avery (avery@diku.dk)
# Content:  Sturmian system: source code for main `configure' script
#
# Process this file with autoconf to produce a configure script.
#
# Created: February 16th 2006. Last modified: September 9th 2008 05:29

#define([svnversion], esyscmd([sh -c "svnversion|tr -d '\n'"]))

define([svnversion], esyscmd([sh -c "date '+%y-%j'|tr -d '\n'"]))

AC_INIT([The Space Discretization Library], [0.1.svnversion], [avery@diku.dk],[libspace])
AC_DEFINE(SVN_REVISION, "svnversion", [SVN Revision])

MAJORVERSION=0
MINORVERSION=1

LONGNAME="Space Discretization Library"
SHORTNAME="libspace"
WEBPAGE="http://libspace.kvante.org"
VERSION=$MAJORVERSION.$MINORVERSION.svnversion
RELEASESTATUS="release"
ISODATE=`date +%Y-%m-%d`
SOURCEDIR=`pwd`


# ----- HAND WAVING ------
AC_PREREQ(2.59)
AC_CONFIG_SRCDIR([./])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([foreign 1.7.3])



# ----- USER OPTIONS -----
AC_PROG_CXX()

AC_ARG_WITH(compiler-type,
	AC_HELP_STRING([--with-compiler-type={gnu|sun|intel|mipspro}],[Type of the C compiler $CC.]),
	compiler_type=$withval
)

AC_ARG_WITH(lapack-ldflags,
	AC_HELP_STRING(
	[--with-lapack-ldflags[=flags]],
	[Supply linker flags to the LAPACK library.]),
	LAPACK_LDFLAGS=$withval
)
AC_ARG_WITH(lapack-cflags,
	AC_HELP_STRING(
	[--with-lapack-cflags[=flags]],
	[Supply compiler flags to the LAPACK library.]),
	LAPACK_CFLAGS=$withval
)
AC_ARG_ENABLE(lapack,
	AC_HELP_STRING([--enable-lapack],[Build with LAPACK. Use --with-lapack-{ldflags,cflags} to specify how.]),
	enable_lapack=$enableval,
	enable_lapack=no
)

PKG_CHECK_MODULES(NETCDF,      [netcdf >= 3],
			       [AC_DEFINE(HAS_NETCDF) REQUIRED_PACKAGES="$REQUIRED_PACKAGES netcdf"])
PKG_CHECK_MODULES(NETCDF_CXX,  [netcdf-c++ >= 3],
			       [AC_DEFINE(HAS_NETCDF_CXX) REQUIRED_PACKAGES="$REQUIRED_PACKAGES netcdf-c++"])
PKG_CHECK_MODULES(LIBMESH,  libmesh,
			    [AC_DEFINE(HAS_LIBMESH) REQUIRED_PACKAGES="$REQUIRED_PACKAGES libmesh"])
PKG_CHECK_MODULES(DEALII,   deal.II,
			    [AC_DEFINE(HAS_DEALII) REQUIRED_PACKAGES="$REQUIRED_PACKAGES deal.II"])
PKG_CHECK_MODULES(DOLFIN,   dolfin,
			    [AC_DEFINE(HAS_DOLFIN) REQUIRED_PACKAGES="$REQUIRED_PACKAGES dolfin"])
PKG_CHECK_MODULES(TRILINOS, [trilinos >= 6],
			    [AC_DEFINE(HAS_TRILINOS) REQUIRED_PACKAGES="$REQUIRED_PACKAGES trilinos"])

AC_SUBST(NETCDF_CFLAGS NETCDF_LIBS LIBMESH_CFLAGS LIBMESH_LIBS DOLFIN_CFLAGS
         DOLFIN_LIBS TRILINOS_CFLAGS TRILINOS_LIBS DEALII_CFLAGS DEALII_LIBS REQUIRED_PACKAGES)


# ----- PROGRAMS -----
AC_PROG_CC
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_LIBTOOL
AC_C_INLINE
LT_INIT


# ----- LIBRARIES -----
# Check for math-library. Without it, we're not much use!
AC_CHECK_LIB([m], [sin])

# ----- HEADER FILES -----
AC_FUNC_ALLOCA
AC_HEADER_STDC
AC_CHECK_HEADERS([limits.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_SIZE_T
# What is the bit order?
AC_C_BIGENDIAN
# Is inline supported?
AC_C_INLINE

# ----- LIBRARY FUNCTIONS -----
AC_FUNC_STRTOD
AC_FUNC_VPRINTF
AC_CHECK_FUNCS([floor memset pow sqrt strchr strpbrk strspn strtol])


# ----- PUTTING IT ALL TOGETHER ----
if test x$compiler_type == x
then
	case $host_os
	in
	  *solaris*) compiler_type=sun 
	  ;;
	  *irix*)    compiler_type=mipspro 
	  ;;
	  *)         compiler_type=gnu 
	  ;;
	esac
	AC_MSG_NOTICE([Compiler-type not specified -- assuming default for $host_os : $compiler_type ])
fi

case $compiler_type
in
 gnu) OMP_CFLAGS="-fopenmp" 
      OMP_LDFLAGS="-lgomp"
      DEFAULT_LAPACK_LDFLAGS="-llapack -lgfortran"
      DEFAULT_LAPACK_CFLAGS="-DLAPACK_HAS_NO_CLAPACK"
      ;;
 sun) 
      OMP_CFLAGS="-xopenmp=parallel -xvpara"
      OMP_LDFLAGS="-lumem"
      DEFAULT_LAPACK_CFLAGS="-dalign -DLAPACK_IS_SUNPERF"
      DEFAULT_LAPACK_LDFLAGS="-xlibmil -xlibmopt -xlic_lib=sunperf"
      ;;
 intel) 
	OMP_CFLAGS="-openmp" 
	DEFAULT_LAPACK_CFLAGS="-DLAPACK_IS_MKL"
	DEFAULT_LAPACK_LDFLAGS="-L$LIBRARY_PATH -lmkl_intel_thread -lmkl_intel_lp64 -lmkl_lapack -lmkl_core"
	DEFAULT_CFLAGS="-fast -fno-alias"
	DEFAULT_LDFLAGS="-Bdynamic"
      ;;
esac

if test x"$CFLAGS" == x
then
	CFLAGS=$DEFAULT_CFLAGS
fi
if test x"$LDFLAGS" == x
then
	LDFLAGS=$DEFAULT_LDFLAGS
fi
if test x"$LAPACK_CFLAGS" == x
then
	LAPACK_CFLAGS=$DEFAULT_LAPACK_CFLAGS
fi
if test x"$LAPACK_LDFLAGS" == x
then
	LAPACK_LDFLAGS=$DEFAULT_LAPACK_LDFLAGS
fi

AC_MSG_NOTICE([LAPACK_CFLAGS=$LAPACK_CFLAGS, LAPACK_LDFLAGS=$LAPACK_LDFLAGS]) 


# Is LAPACK available?
# TODO: Offer choice between several installed LAPACK if they exist -- sequential, parallel, etc.
if test x$lapackpath != xno
then 
 if test x$USE_LAPACK != xno
  then
	AC_DEFINE(USE_LAPACK)
  else
	AC_MSG_WARN([[ 

	I was unable to find a LAPACK library on your system. Because of this,
	I will use slow, homerolled routines for calculating determinants and
	solving eigensystems.

	If you do have LAPACK installed, you may specify its location by
	calling 'configure' with the parameter
	
		--enable-lapack=(path to liblapack)
		--with-lapack-ldflags=(flags)
		--with-lapack-cflags=(flags)
  ]])
 fi
fi
AM_CONDITIONAL(USE_LAPACK,test x$USE_LAPACK = xyes)


AC_SUBST(SOURCEDIR LONGNAME WEBPAGE SHORTNAME RELEASESTATUS OPTIONAL_TARGETS ISODATE MAJORVERSION MINORVERSION VERSION EXTRA_SUBDIRS LAPACK_LDFLAGS LAPACK_CFLAGS USE_LAPACK)

AC_CONFIG_FILES([Makefile
	         README
		 libspace.pc])

AC_OUTPUT
