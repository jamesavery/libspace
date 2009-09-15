
#include <space/fem/fespace.h>
#include "fem/fespace_pointwise.cc"

#include <space/fem/deal.II/mesh.h>
#include "fem/deal.II/mesh.cc"

template class dealii::FESpace<1>;
template class dealii::FESpace<2>;
template class dealii::FESpace<3>;
