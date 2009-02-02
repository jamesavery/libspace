LIBMESH_CXXFLAGS=$(shell pkg-config --cflags libmesh)
LIBMESH_LDFLAGS=$(shell pkg-config --libs libmesh)
DEALII_CXXFLAGS=$(shell pkg-config --cflags deal.II)
DEALII_LDFLAGS=$(shell pkg-config --libs deal.II) -lgfortran
NETCDF_LDFLAGS=$(shell pkg-config --libs netcdf-c++)
NETCDF_CXXFLAGS=$(shell pkg-config --cflags netcdf-c++)


CXX=icpc
OPTFLAGS=${ICPC_OPTFLAGS}
#CXXFLAGS+=-I.. $(OPTFLAGS)
CXX=g++
CXXFLAGS+=-I.. -g3

#CXXFLAGS+=$(LIBMESH_CXXFLAGS) $(OPTFLAGS) -I..
#LDFLAGS+=$(LIBMESH_LDFLAGS)

#CXXFLAGS+=$(DEALII_CXXFLAGS) $(OPTFLAGS) -I..
#LDFLAGS+=$(DEALII_LDFLAGS)

CXXFLAGS+=$(DEALII_CXXFLAGS) 
LDFLAGS+=$(DEALII_LDFLAGS)

CXXFLAGS+=$(NETCDF_CXXFLAGS) 
LDFLAGS+=$(NETCDF_LDFLAGS) -lm

etags:
	find . -name "*.cc" -or -name "*.h" -or -name "*.c" -or -name "*.cpp" -or -name "*.C" | xargs etags -l c++

empty:  misc.o grid/regulargrid.o fem/fespace_pointwise.o fem/deal.II/mesh.o instances.o test-dealii.o
	$(CXX) -o $@ $(CXXFLAGS) $(LDFLAGS)  $^ -Wall


test: misc.o grid/regulargrid.o fem/fespace_pointwise.o fem/deal.II/mesh.o instances.o test.o
	$(CXX) -o $@ $(CXXFLAGS) $(LDFLAGS) $^ -Wall


test-dealii: misc.o grid/regulargrid.o fem/fespace_pointwise.o fem/deal.II/mesh.o instances.o test-dealii.o	
	$(CXX) -o $@ $(CXXFLAGS) $(LDFLAGS)  $^ -Wall


clean:
	find . -\( -name "*~" -or -name "\#*\#" -or -name "*.o" -\) -exec rm {} \;

clean-output:
	rm -f *.vtk *.gpl *.eps
	
