#include <space/fem/deal.II/mesh.h>
#include <space/function.h>
#include <numerics/error_estimator.h>
#include <fe/fe_values.h>
#include <fe/fe_q.h>

#include <math.h>
#include <sstream>
#include <iostream>
#include <fstream>

using namespace std;

typedef struct {
  double error_cutoff;
  size_t fe_order, gauss_order, dimensions;
  size_t initial_refinement, max_refinement_steps;
  bool write_mesh;
} options_t;


template <int dim> class CoulombAttraction : public Function<dim,double>
{
public: 
  typedef typename Function<dim,double>::coordinate coordinate;
  const double Z;
  CoulombAttraction(double Z) : Z(Z) {}

  double operator() (const coordinate& x,
		     const off_t  component = 0) const {

    const double r = sqrt(x.dot(x));
    
    return r==0? -INFINITY : -Z/r;
  }
};

// template <int dim> class Repulsion : public Function<dim,double>
// {
// public: 
//   typedef typename Function<dim,double>::coordinate coordinate;
//   const Function<dim,double> *density;
//   CoulombAttraction(const Function<dim,double> *density) : density(density) {}

//   double operator() (const coordinate& x,
// 		     const off_t  component = 0) const {

//     const double r = sqrt(x.dot(x));
    
//     return r==0? -INFINITY : -Z/r;
//   }
// };

template <int dim> void test_hydrogen(const options_t& options)
{
  using namespace dealii;

  const double upperleft_[3]  = {-10,-10,-10};
  const double dimensions_[3] = {20,20,20};
  const size_t npts[3] = {options.initial_refinement,options.initial_refinement,options.initial_refinement};

  typedef dealii::FESpace<dim> grid;
 
  typename grid::coordinate zero;
  typename grid::coordinate upperleft(upperleft_);
  typename grid::coordinate dimensions(dimensions_);

  grid G(npts,upperleft,dimensions,options.fe_order,options.gauss_order);
  
  
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
  case 1: test_hydrogen<1>(options); break;
  case 2: test_hydrogen<2>(options); break;
  case 3: test_hydrogen<3>(options); break;
  default:
    break;
  }
  return 0;
}
