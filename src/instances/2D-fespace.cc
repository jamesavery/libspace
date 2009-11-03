
#include <space/fem/fespace.h>
#include "fem/fespace_pointwise.cc"

#include <space/fem/deal.II/mesh.h>
#include "fem/deal.II/mesh-init.cc"
#include "fem/deal.II/mesh-output.cc"
#include "fem/deal.II/mesh.cc"

template class dealii::FESpace<2>;
