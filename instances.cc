// Instantiation of the most commonly used numbers of dimensions
#include <libspace/grid/regulargrid.h>
#include <libspace/grid/regulargrid.cc>

template class RegularGrid<1>;
template class RegularGrid<2>;
template class RegularGrid<3>;

#include <libspace/fem/fespace.h>
#include <libspace/fem/fespace_pointwise.cc>

#include <libspace/fem/deal.II/mesh.h>
#include <libspace/fem/deal.II/mesh.cc>

template class dealII::FESpace<1>;
template class dealII::FESpace<2>;
template class dealII::FESpace<3>;
