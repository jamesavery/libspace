#include <libspace/fem/deal.II/mesh.h>
#include <numerics/vectors.h>
#include <numerics/fe_field_function.h>
#include <base/function.h>

#include <stdio.h>

#  include <iostream>
#  include <fstream>
#  include <grid/grid_out.h>
#  include <numerics/data_out.h>
#  include <lac/sparsity_pattern.h>
#  include <grid/grid_refinement.h>

// Auxiliary stuff. Perhaps move to separate file.
int lookup_format(const string *supported_formats, const string& path);

namespace dealII {
  using namespace dealii;
  using namespace std;

#define fespace_member(returntype) template <int dim> returntype FESpace<dim>::

  fespace_member(void) update_nonuniform()
  {
    printf("Updating hanging nodes.\n");
    hanging_node_constraints.clear ();
    DoFTools::make_hanging_node_constraints (dof_handler, hanging_node_constraints);
    hanging_node_constraints.close ();
    hanging_node_constraints.condense (sparsity_pattern);
  }

  fespace_member(void) update()
  {
    // Distribute degrees of freedom
    printf("Distributing DoFs.\n");
    dof_handler.distribute_dofs(fe);  
    
    n_cells   = triangulation.n_cells();
    n_q_pts   = quadrature_formula.size();
    n_dof     = fe_values.dofs_per_cell;   // Rename! n_cell_dofs
    n_weights = dof_handler.n_dofs();	   // Rename! n_dofs

    printf("Allocating sparsity pattern.\n");
    sparsity_pattern.reinit (dof_handler.n_dofs(),
			     dof_handler.n_dofs(),
			     dof_handler.max_couplings_between_dofs());
    printf("Building and compressing sparsity pattern.\n");
    DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
    sparsity_pattern.compress(); 
    
    // Create overlap matrix ("mass matrix") \f$m_{ij} = \int_\Omega \phi_i(x) \phi_j(x) dx\f$ 
    // and Laplacian matrix  \f$m_{ij} = \int_\Omega \nabla\phi_i(x)\cdot \nabla\phi_j(x) dx\f$
    printf("Building overlap matrix.\n");
    dealii::ConstantFunction<dim> one(1.0);
    overlap_matrix.reinit(sparsity_pattern);
    fe_weights.resize(dof_handler.n_dofs());
    MatrixTools::create_mass_matrix(dof_handler,quadrature_formula,overlap_matrix,one,fe_weights.coefficients);
    printf("Building Laplace matrix.\n"); // Slow; Should only generate laplace_matrix when asked; 
				          // Perhaps have a PreparePoisson (Also builds preconditioner?)
    laplace_matrix.reinit(sparsity_pattern);
    MatrixTools::create_laplace_matrix(dof_handler,quadrature_formula,laplace_matrix);

    printf("n_dofs  = %d\n"
	   "n_cells = %d\n"
	   "dofs_per_cell = %d\n",
	   dof_handler.n_dofs(),
	   triangulation.n_cells(),
	   fe.dofs_per_cell);

    printf("Updating PointFunction auxiliary data.\n");
    get_positions();
  }
  

  // TODO: What does weight do exactly? Perhaps get rid of it!
  fespace_member(void)
  LoadFunctionToMesh(double weight, const ScalarFunction& f, 
				     const coordinate& center, FEFunction& density) const
  {
    printf("LoadFunctionToMesh()\n");
    ScalarFunctionWrap F(f,&center); 
    
    ConstraintMatrix cm;
    cm.close();
    QGauss<dim> quadrature(quadrature_order+2);
    density.resize(dof_handler.n_dofs());
    VectorTools::project(dof_handler,cm, quadrature, F, density.coefficients);
  }

  fespace_member(void)
  LoadGaussianToMesh(double weight, double exponent, 
		     const coordinate& center, 
		     double eps, FEFunction& density) const
  {
    GaussianFunction<dim> f2(exponent);
    ClampedGaussianFunction<dim> f(exponent, eps/weight);

    LoadFunctionToMesh(weight, f, center, density);
  }

  fespace_member(void)
  LoadScalar1DFunctionToMesh(FEFunction& density,  Scalar1DFunctionClass& f, 
			     double range,  coordinate& center,double weight) const
  {
    RadialFunction<dim> f3D(f);
    LoadFunctionToMesh(weight, f3D, center, density);  
  }


  fespace_member(double) Integrate(const FEFunction& f) const 
  {
    double sum = 0;
    typename FEFunction::const_iterator c = f.begin();
    typename FEFunction::const_iterator weight = fe_weights.begin();

    for(; c != f.end(); c++,weight++) sum += (*c)*(*weight);

    return sum;
  }

  fespace_member(double) Integrate(const FEFunction& f, const FEFunction& g) const 
  {
    return overlap_matrix.matrix_scalar_product(f.coefficients,g.coefficients);
  }

  fespace_member(double)
  Integrate(const FEFunction& f, const FEFunction& V, const FEFunction& g) const
  {
    return 0;			// TODO: FEFunction -> PointFunction; 
  }

  fespace_member(void)
  ConstructPointFunction(const FEFunction& f, PointFunction& fp) const
  {
    FEValues<dim> fe_values(fe, quadrature_formula, update_values | update_quadrature_points); // evt. update_values
    vector<double> f_values(n_q_pts);

    fp.resize(n_cells*n_q_pts);
    typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();

    for (size_t pt=0; cell!=endc; ++cell){
      fe_values.reinit (cell);
      fe_values.get_function_values(f.coefficients, f_values);
      
      for(size_t q=0;q<n_q_pts;++q,++pt) fp[pt] = f_values[q];
    }    
  }
  fespace_member(void) get_positions()
  {
    FEValues<dim> fe_values(fe, quadrature_formula, update_quadrature_points|update_values|update_JxW_values);

    point_positions.resize(n_cells*n_q_pts);
    point_weights.resize(n_cells*n_q_pts);
    typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();

    for (size_t pt=0; cell!=endc; ++cell){
      fe_values.reinit (cell);
      const vector< Point< dim > > &q_points = fe_values.get_quadrature_points();

      for(size_t q=0;q<n_q_pts;++q,++pt){
	double w = 0;
      	double *x = point_positions[pt].x;

	for(size_t i=0;i<dim;i++) x[i] = q_points[q][i];
	for(size_t i=0;i<n_dof;i++) w += fe_values.shape_value(i,q)*fe_values.JxW(q);
	point_weights[pt] = w;
      }
    }    
  }

  fespace_member(double)
  Integrate(const FEFunction& f, const FEOperator& V, const FEFunction& g) const
  {
//     SparsityPattern::const_iterator ij = sparsity_pattern.begin(), end = sparsity_pattern.end();
//     double sum = 0;
//     for(;ij != end;ij++){
//       size_t i = ij->row(), j = ij->column();
//       sum += f[i]*g[j]*V(i,j);
//     }
//     return sum;
    return V.matrix_scalar_product(f.coefficients,g.coefficients);
  }


  fespace_member(double)
  LaplaceElement(const size_t i, const size_t j) const 
  {
    return laplace_matrix(i,j);
  }

  fespace_member(void)
  SolvePoisson(const FEFunction& density, FEFunction& result) 
  {
    FEValues<dim> fe_values(fe, quadrature_formula, update_quadrature_points | update_values | update_JxW_values);

    typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();

    const size_t N_q   = quadrature_formula.size();
    const size_t N_dof = fe.dofs_per_cell;

    result.resize(dof_handler.n_dofs());
    Vector<double>       rhs(dof_handler.n_dofs());
    FullMatrix<double>   cell_matrix (N_dof, N_dof);
    vector<double>       cell_rho (N_q);
    double               dof_rho[N_dof];
    vector<size_t>       local_dof_indices (N_dof);


    for (; cell!=endc; ++cell){
      fe_values.reinit(cell);
      fe_values.get_function_values(density.coefficients,cell_rho);

      for (size_t i=0; i<N_dof; ++i){
	double integral = 0;
	for (size_t q=0; q<N_q; ++q) integral += fe_values.shape_value (i, q)*cell_rho[q]*fe_values.JxW(q);
	dof_rho[i] = integral;
      }

      cell->get_dof_indices(local_dof_indices);
      for (size_t i=0; i<N_dof; ++i) 
	rhs(local_dof_indices[i]) += cell_rho[i];
    }    

    // TODO: Implement real boundary values, given on initialization of FE-space.
    map<unsigned int,double> boundary_values;
    VectorTools::interpolate_boundary_values (dof_handler,
					      0,
					      ZeroFunction<dim>(),
					      boundary_values);

  
    MatrixTools::apply_boundary_values (boundary_values,
					laplace_matrix,
					result.coefficients,
					rhs);

    // SOLVE
    SolverControl solver_control (2000/*max iterations - should depend on problem size?*/, 1e-10/*tolerance*/);
    SolverCG<>    cg (solver_control);

    PreconditionSSOR<> prec;
    prec.initialize(laplace_matrix, 1.2); // magic parameter -- investigate
    cg.solve (laplace_matrix, result.coefficients, rhs, prec);
  
    hanging_node_constraints.distribute (result.coefficients); // distribute solution
  }

  fespace_member(void)
  SolvePoisson(const PointFunction& density, PointFunction& point_result) 
  {
  }

  // Grid stuff
  fespace_member(void) 
  absolute_error_estimate(const FEFunction& fe_function, const ScalarFunction& function, 
			  Vector<double>& error_vector/*[n_active_cells()]*/) const
  {
    QGauss<dim> fine_quadrature(quadrature_formula.size()+1);
    const size_t N_q = fine_quadrature.size();

    error_vector.resize(triangulation.n_active_cells());

    vector<double>       approx_values(N_q), real_values(N_q);
    vector< Point<dim> > qpoints(N_q);

    ScalarFunctionWrap function_dii(function);

    typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
  
    FEValues<dim> fe_values(fe, fine_quadrature, update_quadrature_points | update_values | update_JxW_values);

    for (size_t cellnum=0; cell!=endc; ++cell,++cellnum){
      fe_values.reinit (cell);
      qpoints = fe_values.get_quadrature_points();
      fe_values.get_function_values(fe_function.coefficients, approx_values);
      function_dii.value_list(qpoints,real_values);

      double error = 0;
      for(size_t q=0;q<N_q;q++)
	error += ::fabs((approx_values[q]-real_values[q])*fe_values.JxW(q));

      error_vector[cellnum] = error;
    }
  }

  fespace_member(void) refine_grid(const Vector<double>& estimated_error_per_cell)
  {
    GridRefinement::refine_and_coarsen_fixed_number (triangulation,estimated_error_per_cell, 0.30, 0.10);
    triangulation.execute_coarsening_and_refinement ();
    update();    
  }
  fespace_member(void) refine_grid(const Vector<float>& estimated_error_per_cell)
  {
    GridRefinement::refine_and_coarsen_fixed_number (triangulation,estimated_error_per_cell, 0.30, 0.10);
    triangulation.execute_coarsening_and_refinement ();
    update();    
  }

  // Ouput stuff.
  fespace_member(void) write_dof_sparsity(const string& path) const 
  {
    const string supported_formats_str[] = {"txt","gpl",""};
    enum {TXT,GNUPLOT} supported_formats;

    ofstream file(path.c_str());
    
    switch(lookup_format(supported_formats_str,path)) {
    case TXT:
      sparsity_pattern.print(file);
      break;
    default:
      sparsity_pattern.print_gnuplot(file);
    };

    file.close();
  }



  fespace_member(void) write_mesh(const string& path) const 
  {
    const string supported_formats_str[] = {"dx","msh","ucd","eps","xfig","gpl",""};
    enum {DX,MSH,UCD,EPS,XFIG,GNUPLOT} supported_formats;

    ofstream file(path.c_str());
    GridOut out;

    switch(lookup_format(supported_formats_str,path)) {
    case DX:
      out.write_dx(triangulation,file);
      break;
    case MSH:
      out.write_msh(triangulation,file);
      break;
    case UCD:
      out.write_ucd(triangulation,file);
      break;
    case EPS:
      out.write_eps(triangulation,file);
      break;
    case XFIG:
      out.write_xfig(triangulation,file);
      break;
    case GNUPLOT:
    default:
      out.write_gnuplot(triangulation,file);
    };

    file.close();
  }




  fespace_member(void) write_function(const string& path, const FEFunction& f) const
  {
    const string supported_formats_str[] = {"dx","eps","gmv","pov","plt","pltx","ucd","vtk","dealII","gpl",""};
    enum {DX,EPS,GMV,POVRAY,TECPLOT,TECPLOT_BINARY,UCD,VTK,DEALII,GNUPLOT} supported_formats;

    ofstream file(path.c_str());
    DataOut<dim> out;

    out.attach_dof_handler (dof_handler);
    out.add_data_vector (f.coefficients, "f");
    out.build_patches ();
   
    switch(lookup_format(supported_formats_str,path)) {
    case DX:               out.write_dx(file);       break;
    case EPS:              {
      DataOutBase::EpsFlags eps_flags;
      eps_flags.z_scaling = 4;

      out.set_flags (eps_flags);
      out.write_eps(file);  
    }
      break;
    case GMV:              out.write_gmv(file);      break;
    case POVRAY:           out.write_povray(file);   break;
    case TECPLOT:          out.write_tecplot(file); break;
    case TECPLOT_BINARY:   out.write_tecplot_binary(file); break;
    case UCD:              out.write_ucd(file);      break;
    case VTK:              out.write_vtk(file);      break;
    case DEALII:           out.write_deal_II_intermediate(file);  break;
    case GNUPLOT:
    default:
      out.write_gnuplot(file);
    };

    file.close();  
  }

}
#undef fespace_member


