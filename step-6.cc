#include <math.h>
#include <base/quadrature_lib.h>
#include <base/function.h>
#include <base/logstream.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>
#include <grid/tria.h>
#include <dofs/dof_handler.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria_boundary_lib.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <fe/fe_values.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <numerics/data_out.h>

#include <fstream>
#include <iostream>

#include <fe/fe_q.h>
#include <grid/grid_out.h>
#include <dofs/dof_constraints.h>
#include <grid/grid_refinement.h>
#include <numerics/error_estimator.h>

using namespace dealii;
template <int dim>
class LaplaceProblem 
{
public:
  LaplaceProblem (const Function<dim>& rhsfun, const Function<dim> *coefficient = NULL);
  ~LaplaceProblem ();

  void run ();

  const Function<dim>& rhsfun;
  const Function<dim> *coefficient;
  
  double Integrate(const Vector<double>& fe_fun) const;
private:
  void setup_system ();
  void assemble_system ();
  void solve ();
  void refine_grid ();
  void output_results (const unsigned int cycle) const;

  Triangulation<dim>   triangulation;

  DoFHandler<dim>      dof_handler;
  FE_Q<dim>            fe;

  ConstraintMatrix     hanging_node_constraints;

  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;

  Vector<double>       solution;
  Vector<double>       system_rhs;
};

template <int dim>
class Coefficient : public Function<dim> 
{
  public:
    Coefficient () : Function<dim>() {}
    
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const
  {
    return 1;
  }
};


template <int dim> class GaussianCharge : public Function<dim>
{
public: 
  const double Q, sigma;
  GaussianCharge(double Q, double sigma) : Q(Q), sigma(sigma) {}

  virtual double value (const Point<dim>   &p,
			const unsigned int  component = 0) const {

    double norm = sigma * sqrt(2.0 * M_PI);
    norm        = norm*norm*norm;

    return (Q/norm)*expl(-p.square()/(2.0*sigma*sigma));
  }
};


template <int dim> class GaussianPotential : public Function<dim>
{
public: 
  const double Q, sigma;
  GaussianPotential(const GaussianCharge<dim>& q) : Q(q.Q), sigma(q.sigma) {}
  GaussianPotential(double Q, double sigma) : Q(Q), sigma(sigma) {}

  virtual double value (const Point<dim>   &p,
			const unsigned int  component = 0) const {
    double r = sqrt(p.square());

    return (Q/(4. * M_PI * r)) * erf(r /( sqrt(2) * sigma ));
  }
};



template <int dim> LaplaceProblem<dim>::LaplaceProblem (const Function<dim>& rhsfun, const Function<dim> *coefficient) : dof_handler (triangulation), fe (2), rhsfun(rhsfun), coefficient(coefficient) {}
template <int dim> LaplaceProblem<dim>::~LaplaceProblem (){  dof_handler.clear (); }


template <int dim>
void LaplaceProblem<dim>::setup_system ()
{
  dof_handler.distribute_dofs (fe);

  sparsity_pattern.reinit (dof_handler.n_dofs(),
			   dof_handler.n_dofs(),
			   dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);

  hanging_node_constraints.clear ();
  DoFTools::make_hanging_node_constraints (dof_handler,
					   hanging_node_constraints);
  hanging_node_constraints.close ();
  hanging_node_constraints.condense (sparsity_pattern);

  sparsity_pattern.compress();

  system_matrix.reinit (sparsity_pattern);
}

template <int dim> double LaplaceProblem<dim>::Integrate(const Vector<double>& fe_fun) const
{
  const QGauss<dim>  quadrature_formula(3);
  FEValues<dim> fe_values(fe, quadrature_formula, update_values | update_JxW_values |update_quadrature_points); 
  const int N_q   = quadrature_formula.size();
  const int N_dof = fe_values.dofs_per_cell;
  std::vector<double> f_values(N_q);

  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

  double w = 0;
  for (size_t pt=0; cell!=endc; ++cell){
    fe_values.reinit (cell);
    fe_values.get_function_values(fe_fun, f_values);
    
    for(size_t i=0;i<N_dof;i++)
      for(size_t q=0;q<N_q;++q,++pt) w += fe_values.shape_value(i,q)*fe_values.JxW(q)*f_values[q];
  }      
  return w;
}

template <int dim> void LaplaceProblem<dim>::assemble_system()
{
  const QGauss<dim>  quadrature_formula(3);
  solution.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());
  

  system_matrix.reinit(sparsity_pattern);
  MatrixTools::create_laplace_matrix(dof_handler,quadrature_formula,system_matrix,rhsfun,system_rhs);

  hanging_node_constraints.condense (system_matrix);
  hanging_node_constraints.condense (system_rhs);

  int n = sqrt(system_matrix.n_nonzero_elements());
  printf("||T||^2 = %g ~ %g (%d)\n",system_matrix.frobenius_norm()/(double)n,
	   system_matrix.frobenius_norm(),n);

  std::map<unsigned int,double> boundary_values;
  VectorTools::interpolate_boundary_values (dof_handler,
					    0,
					    ZeroFunction<dim>(),
					    boundary_values);

  MatrixTools::apply_boundary_values (boundary_values,
				      system_matrix,
				      solution,
				      system_rhs);
}

template <int dim>
void LaplaceProblem<dim>::solve () 
{
  SolverControl           solver_control (1000, 1e-12);
  SolverCG<>              cg (solver_control);

  PreconditionSSOR<> preconditioner;
  preconditioner.initialize(system_matrix, 1.2);

  cg.solve (system_matrix, solution, system_rhs,
	    preconditioner);

  hanging_node_constraints.distribute (solution);

  Vector<double> difference(dof_handler.n_dofs());
  GaussianPotential<dim> exact(1.0,1.0/5.0);
  VectorTools::integrate_difference(dof_handler,solution,exact,difference,QGauss<dim>(3),VectorTools::L2_norm);
  
  Point<dim> x;
  for(int i=0;i<dim;i++) x(i) = 0.1;
  printf("phi(0) = %g\n",exact.value(x));
  printf("||solution-exact|| = %g (\\int f = %g)\n", difference.l2_norm(),Integrate(difference));

}

template <int dim>
void LaplaceProblem<dim>::refine_grid ()
{
  Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

  KellyErrorEstimator<dim>::estimate (dof_handler,
				      QGauss<dim-1>(3),
				      typename FunctionMap<dim>::type(),
				      solution,
				      estimated_error_per_cell);


  GridRefinement::refine_and_coarsen_fixed_number (triangulation,
						   estimated_error_per_cell,
						   0.3, 0.03);

  triangulation.execute_coarsening_and_refinement ();
}


template <int dim>
void LaplaceProblem<dim>::output_results (const unsigned int cycle) const
{
  Assert (cycle < 10, ExcNotImplemented());

  std::string filename = "grid-";
  filename += ('0' + cycle);
  filename += ".eps";
  
  std::ofstream output (filename.c_str());

  GridOut grid_out;
  grid_out.write_eps (triangulation, output);


  {
    DataOut<dim> data_out;

    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (solution, "solution");
    data_out.build_patches ();
  
    char path[0x100];

    sprintf(path,"solution%d.gpl",cycle);
    std::ofstream output (path);
    data_out.write_gnuplot (output);
  }
  {
    DataOut<dim> data_out;

    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (system_rhs, "rho");
    data_out.build_patches ();
  
    char path[0x100];

    sprintf(path,"rho%d.gpl",cycle);
    std::ofstream output (path);
    data_out.write_gnuplot (output);
  }
}


template <int dim>
void LaplaceProblem<dim>::run () 
{
  for (unsigned int cycle=0; cycle<7; ++cycle)
    {
      std::cout << "Cycle " << cycle << ':' << std::endl;

      if (cycle == 0)
	{
	  //	  GridGenerator::hyper_ball (triangulation);
	  GridGenerator::hyper_cube (triangulation,-.5,.5);

	  //	  static const HyperBallBoundary<dim> boundary;
	  //	  triangulation.set_boundary (0, boundary);

	  triangulation.refine_global (4);
	}
      else
	refine_grid ();
      

      std::cout << "   Number of active cells:       "
		<< triangulation.n_active_cells()
		<< std::endl;

      setup_system ();

      std::cout << "   Number of degrees of freedom: "
		<< dof_handler.n_dofs()
		<< std::endl;
      
      assemble_system ();
      solve ();
      output_results (cycle);
    }


  DataOut<dim> data_out;

  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (solution, "solution");
  data_out.build_patches ();
  
  std::ofstream output ("final-solution.gpl");
  data_out.write_gnuplot (output);
}
int main () 
{

  try
    {
      deallog.depth_console (0);
      ConstantFunction<2> one(1.0);
      GaussianCharge<2> charge(1.0,1.0/5.0);
      Point<2> x;
      x(0) = 0;
      x(1) = 0;

      printf("rho(0) = %g\n",charge.value(x));

      LaplaceProblem<2> laplace_problem_2d(charge,&one);
      laplace_problem_2d.run ();
    }

  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      std::cerr << "Exception on processing: " << std::endl
		<< exc.what() << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;

      return 1;
    }

  catch (...) 
    {
      std::cerr << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      std::cerr << "Unknown exception!" << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      return 1;
    }

  return 0;
}
