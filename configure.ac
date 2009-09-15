# This file is part of libspace
#
# Authors:  James Avery (avery@diku.dk)
#
# Process this file with autoconf to produce a configure script.

#define([svnversion], esyscmd([sh -c "svnversion|tr -d '\n'"]))

define([svnversion], esyscmd([sh -c "date '+%y-%j'|tr -d '\n'"]))

AC_INIT([The Space Discretization Library], [0.2.svnversion], [avery@diku.dk],[libspace])
AC_DEFINE(SVN_REVISION, "svnversion", [SVN Revision])

MAJORVERSION=0
MINORVERSION=2
LTAGE=`date +%j | sed 's/^0*//'`

LONGNAME="Space Discretization Library"
SHORTNAME="libspace"
WEBPAGE="http://libspace.kvante.org"
VERSION=$MAJORVERSION.$MINORVERSION.$LTAGE
RELEASESTATUS="release"
ISODATE=`date +%Y-%m-%d`
SOURCEDIR=`pwd`


# ----- HAND WAVING ------
AC_PREREQ(2.59)
AC_CONFIG_SRCDIR([./])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE(libspace,$VERSION)

AC_CONFIG_HEADERS([./include/space/config.h])
AX_PREFIX_CONFIG_H([./include/space/libspace-config.h])


# ----- USER OPTIONS -----
AC_PROG_CXX()

AC_ARG_WITH(dimensions,
	AC_HELP_STRING([--with-dimensions=(1|2|3)+],[List of space-dimensions for which to build. (default=1,2,3)]),
	build_dimensions=$withval,
	build_dimensions="123"
)

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
			       [AC_DEFINE(HAS_NETCDF,1,[System has netcdf]) REQUIRED_PACKAGES="$REQUIRED_PACKAGES netcdf"],[AC_MSG_WARN([No netcdf])])
PKG_CHECK_MODULES(NETCDF_CXX,  [netcdf-c++ >= 3],
			       [AC_DEFINE(HAS_NETCDF_CXX,1,[System has netcdf-c++]) REQUIRED_PACKAGES="$REQUIRED_PACKAGES netcdf-c++"],[AC_MSG_WARN([No netcdf-c++])])
PKG_CHECK_MODULES(LIBMESH,  libmesh,
			       [AC_DEFINE(HAS_LIBMESH,1,[System has the libmesh FEM library (http://libmesh.sf.net)])
			        REQUIRED_PACKAGES="$REQUIRED_PACKAGES libmesh"],[AC_MSG_WARN([No libmesh])])

PKG_CHECK_MODULES(DEALII,   deal.II,
			    [AC_DEFINE(HAS_DEALII,1,[System has the deal.II FEM library (http://dealii.org)]) 
			    REQUIRED_PACKAGES="$REQUIRED_PACKAGES deal.II"])
PKG_CHECK_MODULES(DOLFIN,   dolfin,
   			     [AC_DEFINE(HAS_DOLFIN,1,[System has the Dolfin FEM library (http://fenics.org)])
                              REQUIRED_PACKAGES="$REQUIRED_PACKAGES dolfin"],[AC_MSG_WARN([No dolfin])])

PKG_CHECK_MODULES(TRILINOS, [trilinos >= 6],
			    [AC_DEFINE(HAS_TRILINOS,1,[System has Trilinos (http://trilinos.sandia.gov)]) 
			    REQUIRED_PACKAGES="$REQUIRED_PACKAGES trilinos"])

AC_SUBST(NETCDF_CFLAGS)
AC_SUBST(NETCDF_LIBS)
AC_SUBST(LIBMESH_CFLAGS)
AC_SUBST(LIBMESH_LIBS)
AC_SUBST(DOLFIN_CFLAGS)
AC_SUBST(DOLFIN_LIBS)
AC_SUBST(TRILINOS_CFLAGS)
AC_SUBST(TRILINOS_LIBS)
AC_SUBST(DEALII_CFLAGS)
AC_SUBST(DEALII_LIBS)
AC_SUBST(REQUIRED_PACKAGES)


# ----- PROGRAMS -----
AC_PROG_CC
AC_PROG_MAKE_SET
AC_PROG_LIBTOOL
AC_C_INLINE
LT_INIT


# ----- LIBRARIES -----
# Check for math-library. Without it, we're not much use!
AC_CHECK_LIB([m], [sin])

# ----- HEADER FILES -----
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
AM_CONDITIONAL(BUILD_1D,[echo $build_dimensions | grep 1])
AM_CONDITIONAL(BUILD_2D,[echo $build_dimensions | grep 2])
AM_CONDITIONAL(BUILD_3D,[echo $build_dimensions | grep 3])

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
	AC_DEFINE(USE_LAPACK,1,[A LAPACK (http://netlib.org) implementation exists and is used.])
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


 AC_SUBST(SOURCEDIR)
 AC_SUBST(LONGNAME)
 AC_SUBST(WEBPAGE)
 AC_SUBST(SHORTNAME)
 AC_SUBST(RELEASESTATUS)
 AC_SUBST(OPTIONAL_TARGETS)
 AC_SUBST(ISODATE)
 AC_SUBST(MAJORVERSION)
 AC_SUBST(MINORVERSION)
 AC_SUBST(LTAGE)
 AC_SUBST(VERSION)
 AC_SUBST(EXTRA_SUBDIRS)
 AC_SUBST(LAPACK_LDFLAGS)
 AC_SUBST(LAPACK_CFLAGS)
 AC_SUBST(USE_LAPACK)

AC_CONFIG_FILES([Makefile
		 src/Makefile
		 test/Makefile
	         README
		 libspace.pc])

AC_OUTPUT
