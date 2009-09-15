#ifndef ___INCLUDE_SPACE_LIBSPACE_CONFIG_H
#define ___INCLUDE_SPACE_LIBSPACE_CONFIG_H 1
 
/* ./include/space/libspace-config.h. Generated automatically at end of configure. */
/* ./include/space/config.h.  Generated from config.h.in by configure.  */
/* ./include/space/config.h.in.  Generated from configure.ac by autoheader.  */

/* System has the deal.II FEM library (http://dealii.org) */
#ifndef LIBSPACE_HAS_DEALII 
#define LIBSPACE_HAS_DEALII  1 
#endif

/* System has the Dolfin FEM library (http://fenics.org) */
#ifndef LIBSPACE_HAS_DOLFIN 
#define LIBSPACE_HAS_DOLFIN  1 
#endif

/* System has the libmesh FEM library (http://libmesh.sf.net) */
#ifndef LIBSPACE_HAS_LIBMESH 
#define LIBSPACE_HAS_LIBMESH  1 
#endif

/* System has netcdf */
#ifndef LIBSPACE_HAS_NETCDF 
#define LIBSPACE_HAS_NETCDF  1 
#endif

/* System has netcdf-c++ */
#ifndef LIBSPACE_HAS_NETCDF_CXX 
#define LIBSPACE_HAS_NETCDF_CXX  1 
#endif

/* System has Trilinos (http://trilinos.sandia.gov) */
#ifndef LIBSPACE_HAS_TRILINOS 
#define LIBSPACE_HAS_TRILINOS  1 
#endif

/* Define to 1 if you have the <dlfcn.h> header file. */
#ifndef LIBSPACE_HAVE_DLFCN_H 
#define LIBSPACE_HAVE_DLFCN_H  1 
#endif

/* Define to 1 if you don't have `vprintf' but do have `_doprnt.' */
/* #undef HAVE_DOPRNT */

/* Define to 1 if you have the `floor' function. */
#ifndef LIBSPACE_HAVE_FLOOR 
#define LIBSPACE_HAVE_FLOOR  1 
#endif

/* Define to 1 if you have the <inttypes.h> header file. */
#ifndef LIBSPACE_HAVE_INTTYPES_H 
#define LIBSPACE_HAVE_INTTYPES_H  1 
#endif

/* Define to 1 if you have the `m' library (-lm). */
#ifndef LIBSPACE_HAVE_LIBM 
#define LIBSPACE_HAVE_LIBM  1 
#endif

/* Define to 1 if you have the <limits.h> header file. */
#ifndef LIBSPACE_HAVE_LIMITS_H 
#define LIBSPACE_HAVE_LIMITS_H  1 
#endif

/* Define to 1 if you have the <memory.h> header file. */
#ifndef LIBSPACE_HAVE_MEMORY_H 
#define LIBSPACE_HAVE_MEMORY_H  1 
#endif

/* Define to 1 if you have the `memset' function. */
#ifndef LIBSPACE_HAVE_MEMSET 
#define LIBSPACE_HAVE_MEMSET  1 
#endif

/* Define to 1 if you have the `pow' function. */
#ifndef LIBSPACE_HAVE_POW 
#define LIBSPACE_HAVE_POW  1 
#endif

/* Define to 1 if you have the `sqrt' function. */
#ifndef LIBSPACE_HAVE_SQRT 
#define LIBSPACE_HAVE_SQRT  1 
#endif

/* Define to 1 if you have the <stdint.h> header file. */
#ifndef LIBSPACE_HAVE_STDINT_H 
#define LIBSPACE_HAVE_STDINT_H  1 
#endif

/* Define to 1 if you have the <stdlib.h> header file. */
#ifndef LIBSPACE_HAVE_STDLIB_H 
#define LIBSPACE_HAVE_STDLIB_H  1 
#endif

/* Define to 1 if you have the `strchr' function. */
#ifndef LIBSPACE_HAVE_STRCHR 
#define LIBSPACE_HAVE_STRCHR  1 
#endif

/* Define to 1 if you have the <strings.h> header file. */
#ifndef LIBSPACE_HAVE_STRINGS_H 
#define LIBSPACE_HAVE_STRINGS_H  1 
#endif

/* Define to 1 if you have the <string.h> header file. */
#ifndef LIBSPACE_HAVE_STRING_H 
#define LIBSPACE_HAVE_STRING_H  1 
#endif

/* Define to 1 if you have the `strpbrk' function. */
#ifndef LIBSPACE_HAVE_STRPBRK 
#define LIBSPACE_HAVE_STRPBRK  1 
#endif

/* Define to 1 if you have the `strspn' function. */
#ifndef LIBSPACE_HAVE_STRSPN 
#define LIBSPACE_HAVE_STRSPN  1 
#endif

/* Define to 1 if you have the `strtol' function. */
#ifndef LIBSPACE_HAVE_STRTOL 
#define LIBSPACE_HAVE_STRTOL  1 
#endif

/* Define to 1 if you have the <sys/stat.h> header file. */
#ifndef LIBSPACE_HAVE_SYS_STAT_H 
#define LIBSPACE_HAVE_SYS_STAT_H  1 
#endif

/* Define to 1 if you have the <sys/types.h> header file. */
#ifndef LIBSPACE_HAVE_SYS_TYPES_H 
#define LIBSPACE_HAVE_SYS_TYPES_H  1 
#endif

/* Define to 1 if you have the <unistd.h> header file. */
#ifndef LIBSPACE_HAVE_UNISTD_H 
#define LIBSPACE_HAVE_UNISTD_H  1 
#endif

/* Define to 1 if you have the `vprintf' function. */
#ifndef LIBSPACE_HAVE_VPRINTF 
#define LIBSPACE_HAVE_VPRINTF  1 
#endif

/* Define to the sub-directory in which libtool stores uninstalled libraries.
   */
#ifndef LIBSPACE_LT_OBJDIR 
#define LIBSPACE_LT_OBJDIR  ".libs/" 
#endif

/* Name of package */
#ifndef LIBSPACE_PACKAGE 
#define LIBSPACE_PACKAGE  "libspace" 
#endif

/* Define to the address where bug reports for this package should be sent. */
#ifndef LIBSPACE_PACKAGE_BUGREPORT 
#define LIBSPACE_PACKAGE_BUGREPORT  "avery@diku.dk" 
#endif

/* Define to the full name of this package. */
#ifndef LIBSPACE_PACKAGE_NAME 
#define LIBSPACE_PACKAGE_NAME  "The Space Discretization Library" 
#endif

/* Define to the full name and version of this package. */
#ifndef LIBSPACE_PACKAGE_STRING 
#define LIBSPACE_PACKAGE_STRING  "The Space Discretization Library 0.2.09-258" 
#endif

/* Define to the one symbol short name of this package. */
#ifndef LIBSPACE_PACKAGE_TARNAME 
#define LIBSPACE_PACKAGE_TARNAME  "libspace" 
#endif

/* Define to the version of this package. */
#ifndef LIBSPACE_PACKAGE_VERSION 
#define LIBSPACE_PACKAGE_VERSION  "0.2.09-258" 
#endif

/* Define to 1 if you have the ANSI C header files. */
#ifndef LIBSPACE_STDC_HEADERS 
#define LIBSPACE_STDC_HEADERS  1 
#endif

/* SVN Revision */
#ifndef LIBSPACE_SVN_REVISION 
#define LIBSPACE_SVN_REVISION  "09-258" 
#endif

/* A LAPACK (http://netlib.org) implementation exists and is used. */
#ifndef LIBSPACE_USE_LAPACK 
#define LIBSPACE_USE_LAPACK  1 
#endif

/* Version number of package */
#ifndef LIBSPACE_VERSION 
#define LIBSPACE_VERSION  "0.2.258" 
#endif

/* Define to 1 if your processor stores words with the most significant byte
   first (like Motorola and SPARC, unlike Intel and VAX). */
/* #undef WORDS_BIGENDIAN */

/* Define to empty if `const' does not conform to ANSI C. */
/* #undef const */

/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name.  */
#ifndef __cplusplus
/* #undef inline */
#endif

/* Define to `unsigned int' if <sys/types.h> does not define. */
/* #undef size_t */
 
/* once: ___INCLUDE_SPACE_LIBSPACE_CONFIG_H */
#endif
