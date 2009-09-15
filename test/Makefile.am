ACLOCAL_AMFLAGS=-I m4

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

EXTRALIBS=@DEALII_LIBS@ @NETCDF_LIBS@ @NETCDF_CXX_LIBS@ @LIBMESH_LIBS@ 

test_SOURCES=test.cc
test_CXXFLAGS=$(CXXFLAGS) $(AM_CXXFLAGS)
test_LDFLAGS=$(LDFLAGS) $(AM_LDFLAGS) $(EXTRALIBS)
test_LDADD=$(top_srcdir)/src/libspace.la

# test-dealii.cc uses deal.II directly, thus CFLAGS
test_dealii_SOURCES=test-dealii.cc
test_dealii_CXXFLAGS=$(CXXFLAGS) $(AM_CXXFLAGS) @DEALII_CFLAGS@ 
test_dealii_LDFLAGS=$(LDFLAGS) $(AM_LDFLAGS) $(EXTRALIBS)
test_dealii_LDADD=$(test_LDADD)

test_libmesh_SOURCES=test-libmesh.cc
test_libmesh_CXXFLAGS=$(CXXFLAGS) $(AM_CXXFLAGS) @LIBMESH_CFLAGS@
test_libmesh_LDFLAGS=$(LDFLAGS) $(AM_LDFLAGS) $(EXTRALIBS)
test_libmesh_LDADD=$(test_LDADD)


test_trilinos_SOURCES=test-trilinos.cc
test_trilinos_CXXFLAGS=$(CXXFLAGS) $(AM_CXXFLAGS) @DEALII_CFLAGS@ @TRILINOS_CFLAGS@
test_trilinos_LDFLAGS=$(LDFLAGS) $(AM_LDFLAGS) $(EXTRALIBS)
test_trilinos_LDADD=$(test_LDADD)

# test_pfinterval_SOURCES=test-pfinterval.cc
# test_pfinterval_CXXFLAGS=$(CXXFLAGS) $(AM_CXXFLAGS)
# test_pfinterval_LDFLAGS=$(LDFLAGS) $(AM_LDFLAGS) $(EXTRALIBS)
# test_pfinterval_LDADD=$(top_srcdir)/src/libspace.la 

test_msh_SOURCES=test-msh.cc
test_msh_CXXFLAGS=$(CXXFLAGS) $(AM_CXXFLAGS) @DEALII_CFLAGS@
test_msh_LDFLAGS=$(LDFLAGS) $(AM_LDFLAGS) $(EXTRALIBS)
test_msh_LDADD=$(test_LDADD)


bin_PROGRAMS=
#test test-dealii test-libmesh test-trilinos test-msh

