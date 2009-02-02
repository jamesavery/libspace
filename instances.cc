// Instantiation of the most commonly used numbers of dimensions
#include <libdiscretization/grid/regulargrid.h>
#include <libdiscretization/grid/regulargrid.cc>

template class RegularGrid<1>;
template class RegularGrid<2>;
template class RegularGrid<3>;

#include <libdiscretization/fem/fespace.h>
#include <libdiscretization/fem/fespace_pointwise.cc>

#include <libdiscretization/fem/deal.II/mesh.h>
#include <libdiscretization/fem/deal.II/mesh.cc>

template class dealII::FESpace<1>;
template class dealII::FESpace<2>;
template class dealII::FESpace<3>;
