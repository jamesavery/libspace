#ifndef LIBDISC_FEOPERATOR_H
# define LIBDISC_FEOPERATOR_H

#include <space/fem/deal.II/mesh.h>
/* deal.II */
#include <dofs/dof_handler.h>
#include <fe/fe_values.h>
#include <base/quadrature_lib.h>

#include <lac/full_matrix.h>
#include <lac/solver_gmres.h>
#include <lac/solver_cg.h>
#include <lac/trilinos_block_vector.h>
#include <lac/trilinos_sparse_matrix.h>
#include <lac/trilinos_block_sparse_matrix.h>
#include <lac/trilinos_precondition.h>

#include <lac/full_matrix.h>
#include <lac/solver_gmres.h>
#include <lac/solver_cg.h>
#include <lac/trilinos_block_vector.h>
#include <lac/trilinos_sparse_matrix.h>
#include <lac/trilinos_block_sparse_matrix.h>
#include <lac/trilinos_precondition.h>

/* /deal.II */
/* deal.II+Trilinos */
#include <lac/full_matrix.h>
#include <lac/solver_gmres.h>
#include <lac/solver_cg.h>
#include <lac/trilinos_block_vector.h>
#include <lac/trilinos_sparse_matrix.h>
#include <lac/trilinos_block_sparse_matrix.h>
#include <lac/trilinos_precondition.h>
/* /deal.II+Trilinos */
/* Trilinos */
#include <Epetra_Map.h>
#include <Epetra_FEVector.h>
/* /Trilinos */

namespace dealii {
  typedef PointFunction_Epetra  PointFunction_; 
  typedef FEFunction_Epetra     FEFunction_; 

  template <int dim> class FESpace;

  template <int dim> class FEOperator_ : public SparseMatrix<double> {
  public:
    typedef dealii::FESpace<dim> FESpace;
    const   FESpace& space;

    FEOperator_(const FESpace& G, const typename FESpace::FEFunction *f = NULL) 
      : SparseMatrix<double>(G.sparsity_pattern), space(G) {  if(f != NULL) BuildFromFEFunction(*f);  }

    void BuildFromFEFunction(const typename FESpace::FEFunction& f){
      using namespace dealii;
      // XXX: FEFieldFunction is *much* too slow! ...Wow, still running! 
      //    FEFieldFunction<dim> f_interpolated(dof_handler,f.coefficients);
      //    MatrixTools::create_mass_matrix(dof_handler,quadrature_formula,A,&f_interpolated);

      // NB: The following is more than 100 times faster for 27000 cells (30x30x30 fe_order=1).
      // TODO: This is still a slow way. Precalculate shape_value(i,q), perhaps.
      typename DoFHandler<dim>::active_cell_iterator
	cell = space.dof_handler.begin_active(),
	endc = space.dof_handler.end();
      
      const size_t N_dof   = space.n_cell_dof;
      const size_t N_q     = space.quadrature_formula.size();

      FullMatrix<double>   cell_matrix (N_dof,N_dof);
      FEValues<dim> fe_values(space.fe,space.quadrature_formula,update_quadrature_points|update_JxW_values|update_values);
      vector<double> f_values(N_q);
      vector<unsigned int> local_dof_indices (N_dof);

      for (size_t pt=0; cell!=endc; ++cell){
	fe_values.reinit (cell);
	fe_values.get_function_values(f.coefficients, f_values);

	cell_matrix = 0;
 
	for (size_t i=0; i<N_dof; ++i)
	  for (size_t j=0; j<N_dof; ++j)
	    for (size_t q=0; q<N_q; ++q)
	      cell_matrix(i,j) += (f_values[q] *
				   fe_values.shape_value(i, q) *
				   fe_values.shape_value(j, q) *
				   fe_values.JxW (q));
 
	cell->get_dof_indices (local_dof_indices);
 
	for (size_t i=0; i<N_dof; ++i)
	  for (size_t j=0; j<N_dof; ++j)
	    add (local_dof_indices[i],local_dof_indices[j], cell_matrix(i,j));
      }
    }

  };
}

#endif
