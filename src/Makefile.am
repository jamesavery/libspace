
AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

LIBSPACE_LT_CURRENT =  @MAJORVERSION@
LIBSPACE_LT_REVISION = @MINORVERSION@
LIBSPACE_LT_AGE =      0

lib_LTLIBRARIES = libspace.la

INCLUDE = \
$(top_srcdir)/include/space/libspace-config.h \
$(top_srcdir)/include/space/discretization.h \
$(top_srcdir)/include/space/function.h

INCLUDE_FEM=\
$(top_srcdir)/include/space/fem/fespace.h 

INCLUDE_FEM_DEALII =\
$(top_srcdir)/include/space/fem/deal.II/mesh.h \
$(top_srcdir)/include/space/fem/deal.II/feoperator.h 

INCLUDE_FEM_LIBMESH =\
$(top_srcdir)/include/space/fem/libmesh/feoperator.h\
$(top_srcdir)/include/space/fem/libmesh/mesh.h

INCLUDE_GRID=\
$(top_srcdir)/include/space/grid/regulargrid.h

INCLUDE_STORAGE=\
$(top_srcdir)/include/space/storage/storage.h\
$(top_srcdir)/include/space/storage/simple.h\
$(top_srcdir)/include/space/storage/smallvector.h\
$(top_srcdir)/include/space/storage/smallmatrix.h

libspace_la_SOURCES = \
misc.cc grid/regulargrid.cc fem/fespace_pointwise.cc \
fem/deal.II/mesh.cc instances.cc                     \
$(INCLUDE) $(INCLUDE_GRID) $(INCLUDE_STORAGE) $(INCLUDE_FEM)\
$(INCLUDE_FEM_DEALII) $(INCLUDE_FEM_LIBMESH)

libspace_la_LDFLAGS = $(NO_UNDEFINED_FLAG) -version-info			\
$(LIBSPACE_LT_CURRENT):$(LIBSPACE_LT_REVISION):$(LIBSPACE_LT_AGE)	\
$(LDFLAGS) $(AM_LDFLAGS) $(EXTRALIBS)
libspace_la_CXXFLAGS = $(CXXFLAGS) $(AM_CXXFLAGS) @DEALII_CFLAGS@ @NETCDF_CFLAGS@ @NETCDF_CXX_CFLAGS@ @LIBMESH_CFLAGS@ 


library_includedir=$(includedir)/space
library_include_HEADERS=$(INCLUDE)

library_grid_includedir=$(includedir)/space/grid
library_grid_include_HEADERS=$(INCLUDE_GRID)

library_storage_includedir=$(includedir)/space/storage
library_storage_include_HEADERS=$(INCLUDE_STORAGE)

library_fem_includedir=$(includedir)/space/fem
library_fem_include_HEADERS=$(INCLUDE_FEM)

library_fem_dealii_includedir=$(includedir)/space/fem/deal.II
library_fem_dealii_include_HEADERS=$(INCLUDE_FEM_DEALII)

library_fem_libmesh_includedir=$(includedir)/space/fem/libmesh
library_fem_libmesh_include_HEADERS=$(INCLUDE_FEM_LIBMESH)

