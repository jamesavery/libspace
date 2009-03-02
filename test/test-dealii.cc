#include <space/fem/deal.II/mesh.h>
#include <space/function.h>
#include <numerics/error_estimator.h>
#include <fe/fe_values.h>
#include <fe/fe_q.h>

#include <math.h>
#include <sstream>
#  include <iostream>
#  include <fstream>
#  include <grid/grid_out.h>
#  include <numerics/data_out.h>
#  include <grid/grid_refinement.h>

using namespace std;

typedef struct {
  double error_cutoff;
  size_t fe_order, gauss_order, dimensions;
  size_t initial_refinement, max_refinement_steps;
  bool write_mesh;
} options_t;




template <int dim> class GaussianCharge : public Function<dim,double>
{
public: 
  typedef typename Function<dim,double>::coordinate coordinate;
  const double Q, sigma;
  GaussianCharge(double Q, double sigma) : Q(Q), sigma(sigma) {}

  double operator() (const coordinate& x,
		     const off_t  component = 0) const {

    double norm = sigma * sqrt(2.0 * M_PI);
    norm        = norm*norm*norm;

    return (Q/norm)*expl(-x.dot(x)/(2.0*sigma*sigma));
  }
};


template <int dim> class GaussianPotential : public Function<dim,double>
{
public: 
  typedef typename Function<dim,double>::coordinate coordinate;

  const double Q, sigma;
  GaussianPotential(const GaussianCharge<dim>& q) : Q(q.Q), sigma(q.sigma) {}
  GaussianPotential(double Q, double sigma) : Q(Q), sigma(sigma) {}

  double operator() (const coordinate& x, const off_t  component = 0) const {
    double r = sqrt(x.dot(x));

    return (Q/(4. * M_PI * r)) * erf(r /( sqrt(2) * sigma ));
  }
};

template <int dim> void test_refinement(const options_t& options)
{
  using namespace dealii;
  const double exact[3] = {0.546291971785147991730352057218L,
			   0.298434918436904929474588820641L,
			   0.163032600042436609935236443504L};
  const double upperleft_[3]  = {-.5,-.5,-.5};
  const double dimensions_[3] = {1,1,1};
  const size_t npts[3] = {options.initial_refinement,options.initial_refinement,options.initial_refinement};

  typedef dealii::FESpace<dim> grid;

  typename grid::coordinate zero;
  typename grid::coordinate upperleft(upperleft_);
  typename grid::coordinate dimensions(dimensions_);
  typename grid::FEFunction density;

  grid G(npts,upperleft,dimensions,options.fe_order,options.gauss_order);


  GaussianCharge<dim>    gauss(1.0,1.0/5.0);
  double relative_error = INFINITY;
  typename grid::FEFunction vHartree;
  
  for(size_t refinement_step = 0; refinement_step < options.max_refinement_steps 
	                      && fabs(relative_error)>options.error_cutoff; refinement_step++){
    ostringstream path;
    path << "step" << refinement_step << ".eps";

    Vector<float> error_estimate;

    if(options.write_mesh) G.write_mesh(path.str());
    G.LoadFunctionToMesh(1.0,gauss,zero,density);
    double normsqr = G.Integrate(density,density);
    relative_error = (exact[dim-1]-normsqr)/exact[dim-1];

    printf("%d: <f,f> =~ %g (relative error = %g; %d dofs)\n",
	   refinement_step, normsqr, relative_error, G.dof_handler.n_dofs());

    if(refinement_step < options.max_refinement_steps && fabs(relative_error)>options.error_cutoff){
      printf("Refining grid...\n");
      G.absolute_error_estimate(density,gauss,error_estimate);
      G.refine_grid(error_estimate);
      printf("Done.\n");
    } 
  }
}

template <int dim> void test_poisson(const options_t& options)
{
  using namespace dealii;

  const double upperleft_[3]  = {-.5,-.5,-.5};
  const double dimensions_[3] = {1,1,1};
  const size_t npts[3] = {options.initial_refinement,options.initial_refinement,options.initial_refinement};

  typedef dealii::FESpace<dim> grid;

  typename grid::coordinate zero;
  typename grid::coordinate upperleft(upperleft_);
  typename grid::coordinate dimensions(dimensions_);

  grid G(npts,upperleft,dimensions,options.fe_order,options.gauss_order);

  GaussianCharge<dim> gauss(1.0,1.0/5.0);
  GaussianPotential<dim> exact_solution(1.0,1.0/5.0);
  typename grid::ScalarFunctionWrap exact(exact_solution);
  
  double relative_error = INFINITY;
  typename grid::FEFunction fe_gauss, solution;


  for(size_t refinement_step = 0; refinement_step < options.max_refinement_steps 
	                      && fabs(relative_error)>options.error_cutoff; refinement_step++){
    typename grid::FEFunction difference;
    difference.resize(G.n_dofs);

    G.LoadFunctionToMesh(1.0,gauss,zero,fe_gauss);

    G.write_function("rho.gpl",fe_gauss);

    printf("||rho||^2 = %g\n",G.Integrate(fe_gauss,fe_gauss)); 
    G.SolvePoisson(fe_gauss,solution);

    VectorTools::integrate_difference(G.dof_handler,solution.coefficients,exact,difference.coefficients,
				      G.quadrature_formula, VectorTools::L2_norm);    

    printf("||solution-exact||^2 = %g\n",G.Integrate(difference,difference)); 

    if(options.write_mesh){ 	// XXX: write_mesh -> separat fra write_solution
      ostringstream path;
      path << "poisson" << refinement_step << ".gpl";
      G.write_function(path.str(), solution);
    }

    if(refinement_step+1 < options.max_refinement_steps && fabs(relative_error)>options.error_cutoff){
      Vector<float> error_estimate(G.triangulation.n_active_cells());

      KellyErrorEstimator<dim>::estimate (G.dof_handler,
					  QGauss<dim-1>(options.gauss_order+1),
					  typename FunctionMap<dim>::type(),
					  solution.coefficients,
					  error_estimate);

      relative_error = error_estimate.l2_norm();
      printf("%d: Total error estimate = %g\n",refinement_step,relative_error);
      G.refine_grid(error_estimate);
      //G.refine_grid(1);
    }
  }
    if(options.write_mesh){ 	// XXX: write_mesh -> separat fra write_solution
      ostringstream path;
      path << "poisson-last.gpl";
      G.write_function(path.str(), solution);
    }
}

void test_3d(const options_t& options)
{
  using dealii::FESpace;
  const long double exact[4] = {0.346574791727846916020668859333L,
				0.163032600042436609935236443504L,
				0.0940858601000585939092204780071L,
				0.0619638909335590820510280935045L
  };

  double upperleft_[3]  = {-.5,-.5,-.5};
  double dimensions_[3] = {1,1,1};
  const size_t npts[3] = {10,10,10};

  typedef FESpace<3> grid3d;
  typedef FESpace<2> grid2d;
  typedef FESpace<1> grid1d;

  grid3d::coordinate upperleft(upperleft_);
  grid3d::coordinate dimensions(dimensions_);

  grid3d G(npts,upperleft,dimensions,2,3);
  grid3d::FEFunction density;
  G.LoadGaussianToMesh(1.0,5.0,grid3d::coordinate(),0, density);

  double l1 = G.Integrate(density);
  double l2 = G.Integrate(density,density);
  printf("||f||_1 = %g (err = %Lg)\n",l1,(l1 - exact[0])/exact[0]);
  printf("<f|f>   = %g (err = %Lg)\n",l2,(l2 - exact[1])/exact[1]);

  // Pre-nicification hack
  {
    using namespace dealii;
    printf("Building operator from DoFFunction...\n");
    grid3d::FEOperator V(G, &density);
    printf("Done.\n");
    double l3 = G.Integrate(density,V,density);
    printf("<f|f|f>_3 = %g (err = %Lg)\n",l3,(l3 - exact[2])/exact[2]);

//     printf("Doing 1000 weighted inner products: <f|V|g>\n");
//     for(size_t i=0;i<1000;i++) l3 = G.Integrate(density,V,density);
//     printf("Done.\n");
  }

  {
    grid3d::PointFunction f;
    printf("Constructing approximate point-function (n_cells*n_q_pts = %d points).\n",G.point_weights.size());
    G.ConstructPointFunction(density,f);
    
    double l1 = G.Integrate(f);
    double l2 = G.Integrate(f,f);
    double l3 = G.Integrate(f,f,f);
    double l4 = G.Integrate(f*f*f*f);
    printf("||f||_1 = %g (err = %Lg)\n",l1,(l1 - exact[0])/exact[0]);
    printf("||f||_2 = %g (err = %Lg)\n",l2,(l2 - exact[1])/exact[1]);
    printf("||f||_3 = %g (err = %Lg)\n",l3,(l3 - exact[2])/exact[2]);
    printf("||f||_4 = %g (err = %Lg)\n",l4,(l4 - exact[3])/exact[3]);

//     printf("Doing 10.000 weighted inner products: <f|V|g>\n");
//     for(size_t i=0;i<10000;i++) l3 = G.Integrate(f,f,f);
//     printf("Done.\n");

//     printf("Doing 100.000 simple integrals: \\int_\\Omega f\n");
//     for(size_t i=0;i<100000;i++) l3 = G.Integrate(f);
//     printf("Done.\n");
  }
}


#include <stdio.h>
#include <stdlib.h>
#define MATCH(x)  if(!strncmp(av[i],#x"=",sizeof(#x))){ options.x = strtol(av[i]+sizeof(#x),NULL,10); }
#define MATCHD(x) if(!strncmp(av[i],#x"=",sizeof(#x))){ options.x = strtod(av[i]+sizeof(#x),NULL);    }

int main(int ac, char **av)
{
  options_t options;
  options.error_cutoff = 5e-10;;
  options.fe_order = 1, options.gauss_order = 2, options.dimensions = 1;
  options.initial_refinement = 10, options.max_refinement_steps = 5;
  options.write_mesh = false;

  for(int i=1;i<ac;i++){
    MATCH(dimensions);
    MATCH(fe_order);
    MATCH(gauss_order);
    MATCH(initial_refinement);
    MATCH(max_refinement_steps);
    MATCHD(error_cutoff);
    MATCH(write_mesh);
  }

  printf("dimensions  = %d\n"
	 "fe_order    = %d\n"
	 "gauss_order = %d\n"
	 "initial_refinement   = %d\n"
	 "max_refinement_steps = %d\n"
	 "error_cutoff            = %g\n"
	 "write_mesh  = %s\n",
	 options.dimensions,options.fe_order,options.gauss_order,
	 options.initial_refinement,options.max_refinement_steps,
	 options.error_cutoff,options.write_mesh? "true":"false"
	 );


  switch(options.dimensions){
  case 1: test_poisson<1>(options); break;
  case 2: test_poisson<2>(options); break;
  case 3: test_poisson<3>(options); break;
  default:
    break;
  }
  return 0;
}
