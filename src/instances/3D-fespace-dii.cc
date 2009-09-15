// XXX: Kaempe hack. Hvad gaar galt?
#include <space/fem/fespace.h>
#include <space/fem/deal.II/mesh.h>
#include "fem/fespace_pointwise.cc"

#include <vector>
#include <lac/vector.h>
template class FESpace<3, dealii::PointFunction_, dealii::FEFunction_, dealii::FEOperator_<3> >;
