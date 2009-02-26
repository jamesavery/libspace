// Instantiation of the most commonly used numbers of dimensions
#include <space/grid/regulargrid.h>
#include "grid/regulargrid.cc"

template class RegularGrid<1>;
template class RegularGrid<2>;
template class RegularGrid<3>;

#include <space/fem/fespace.h>
#include "fem/fespace_pointwise.cc"

#include <space/fem/deal.II/mesh.h>
#include "fem/deal.II/mesh.cc"

template class dealii::FESpace<1>;
template class dealii::FESpace<2>;
template class dealii::FESpace<3>;
