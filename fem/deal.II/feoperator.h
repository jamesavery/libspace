#ifndef LIBDISC_FEOPERATOR_H
# define LIBDISC_FEOPERATOR_H

#include <libdiscretization/fem/deal.II/mesh.h>
/* deal.II */
#include  <dofs/dof_handler.h>
#include <lac/sparse_matrix.h>	
#include <fe/fe_values.h>
#include <base/quadrature_lib.h>
/* /deal.II */

namespace dealII {
  typedef PointFunction_Simple<std::vector<double>,double,double>  PointFunction_; 
  typedef FEFunction_Simple<dealii::Vector<double>,double,double>     FEFunction_; 

  template <int dim> class FESpace;

  template <int dim> class FEOperator : public dealii::SparseMatrix<double> {
  public:
    typedef FESpace<dim> FESpace;
    const   FESpace& space;

    FEOperator(const FESpace& G, const typename FESpace::FEFunction *f = NULL) 
      : dealii::SparseMatrix<double>(G.sparsity_pattern), space(G) {  if(f != NULL) BuildFromFEFunction(*f);  }

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
      
      const size_t n_dof   = space.n_dof;
      const size_t n_q_pts = space.quadrature_formula.size();

      FullMatrix<double>   cell_matrix (n_dof,n_dof);
      FEValues<dim> fe_values(space.fe,space.quadrature_formula,update_quadrature_points|update_JxW_values|update_values);
      vector<double> f_values(n_q_pts);
      vector<unsigned int> local_dof_indices (n_dof);

      for (size_t pt=0; cell!=endc; ++cell){
	fe_values.reinit (cell);
	fe_values.get_function_values(f.coefficients, f_values);

	cell_matrix = 0;
 
	for (size_t i=0; i<n_dof; ++i)
	  for (size_t j=0; j<n_dof; ++j)
	    for (size_t q=0; q<n_q_pts; ++q)
	      cell_matrix(i,j) += (f_values[q] *
				   fe_values.shape_value(i, q) *
				   fe_values.shape_value(j, q) *
				   fe_values.JxW (q));
 
	cell->get_dof_indices (local_dof_indices);
 
	for (size_t i=0; i<n_dof; ++i)
	  for (size_t j=0; j<n_dof; ++j)
	    add (local_dof_indices[i],local_dof_indices[j], cell_matrix(i,j));
      }
    }

  };
}

#endif
