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


// XXX: Kaempe hack. Hvad gaar galt?
#include <vector>
#include <lac/vector.h>
template class FESpace<1, PointFunction_Simple<std::vector<double, std::allocator<double> >, double, double>, FEFunction_Simple<dealii::Vector<double>, double, double>, dealii::FEOperator_<1> >;
template class FESpace<2, PointFunction_Simple<std::vector<double, std::allocator<double> >, double, double>, FEFunction_Simple<dealii::Vector<double>, double, double>, dealii::FEOperator_<2> >;
template class FESpace<3, PointFunction_Simple<std::vector<double, std::allocator<double> >, double, double>, FEFunction_Simple<dealii::Vector<double>, double, double>, dealii::FEOperator_<3> >;
