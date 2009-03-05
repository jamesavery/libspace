ACLOCAL_AMFLAGS=-I m4

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

EXTRALIBS=@DEALII_LIBS@ @NETCDF_LIBS@ @NETCDF_CXX_LIBS@ @LIBMESH_LIBS@ 

test_SOURCES=test/test.cc
test_CXXFLAGS=$(CXXFLAGS) $(AM_CXXFLAGS)
test_LDFLAGS=$(LDFLAGS) $(AM_LDFLAGS) $(EXTRALIBS)
test_LDADD=libspace.la 

# test-dealii.cc uses deal.II directly, this CFLAGS
test_dealii_SOURCES=test/test-dealii.cc
test_dealii_CXXFLAGS=$(CXXFLAGS) $(AM_CXXFLAGS) @DEALII_CFLAGS@ 
test_dealii_LDFLAGS=$(LDFLAGS) $(AM_LDFLAGS) $(EXTRALIBS)
test_dealii_LDADD=$(top_srcdir)/libspace.la

bin_PROGRAMS=test test-dealii test-libmesh 

