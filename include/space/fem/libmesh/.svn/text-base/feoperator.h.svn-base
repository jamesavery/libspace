#ifndef LIBDISC_FEOPERATOR_H
# define LIBDISC_FEOPERATOR_H

#include <libmesh/numeric_vector.h>
#include <libmesh/sparse_matrix.h>
#include <space/fem/libmesh/mesh.h>

namespace libMesh {
  typedef PointFunction_Simple<std::vector<double>,double,double>  PointFunction_; 
  typedef FEFunction_Simple<NumericVector<double>,double,double>     FEFunction_; 

  template <int dim> class FESpace;

  template <int dim> class FEOperator : public SparseMatrix<double> {
  public:
    typedef FESpace<dim> FESpace;
    const   FESpace& space;

    FEOperator(const FESpace& G, const typename FESpace::FEFunction *f = NULL) 
      : SparseMatrix<double>(G.sparsity_pattern), space(G) {  if(f != NULL) BuildFromFEFunction(*f);  }

    void BuildFromFEFunction(const typename FESpace::FEFunction& f){
    }

  };
}

#endif
