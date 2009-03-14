
AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

LIBSPACE_LT_CURRENT =  @MAJORVERSION@
LIBSPACE_LT_REVISION = @MINORVERSION@
LIBSPACE_LT_AGE =      0

lib_LTLIBRARIES = libspace.la

INCLUDE_FILES = \
$(top_srcdir)/include/shibuya-wulfman/discretization.h \
$(top_srcdir)/include/shibuya-wulfman/fem/deal.II/mesh.h \
$(top_srcdir)/include/shibuya-wulfman/fem/deal.II/feoperator.h \
$(top_srcdir)/include/shibuya-wulfman/fem/feshibuya-wulfman.h \
$(top_srcdir)/include/shibuya-wulfman/fem/libmesh/autoheader-undef.h\
$(top_srcdir)/include/shibuya-wulfman/fem/libmesh/feoperator.h\
$(top_srcdir)/include/shibuya-wulfman/fem/libmesh/mesh.h\
$(top_srcdir)/include/shibuya-wulfman/function.h\
$(top_srcdir)/include/shibuya-wulfman/grid/regulargrid.h\
$(top_srcdir)/include/shibuya-wulfman/storage/storage.h\
$(top_srcdir)/include/shibuya-wulfman/storage/simple.h\
$(top_srcdir)/include/shibuya-wulfman/storage/smallvector.h\
$(top_srcdir)/include/shibuya-wulfman/storage/smallmatrix.h

libspace_la_SOURCES = \
misc.cc grid/regulargrid.cc fem/fespace_pointwise.cc \
fem/deal.II/mesh.cc instances.cc                     \
$(INCLUDE_FILES)

libspace_la_LDFLAGS = $(NO_UNDEFINED_FLAG) -version-info			\
$(LIBSPACE_LT_CURRENT):$(LIBSPACE_LT_REVISION):$(LIBSPACE_LT_AGE)	\
$(LDFLAGS) $(AM_LDFLAGS) $(EXTRALIBS)
libspace_la_CXXFLAGS = $(CXXFLAGS) $(AM_CXXFLAGS) @DEALII_CFLAGS@ @NETCDF_CFLAGS@ @NETCDF_CXX_CFLAGS@ @LIBMESH_CFLAGS@ 

