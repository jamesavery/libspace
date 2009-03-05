#include <space/grid/regulargrid.h>
#include <math.h>
#include <stdlib.h>

#define grid_member(returntype) template <int dim,class PointFunction> returntype RegularGrid<dim,PointFunction>::

grid_member(void) LoadFunctionToMesh(double weight, const ScalarFunction& f, 
			  const coordinate& position, PointFunction& density) const
{
  density.resize(total_npts);
  typename PointFunction::iterator rho = density.coefficients.begin(), end = density.coefficients.end();
  for(size_t i=0; rho != end; rho++,i++){ // XXX: TODO: i should be global index, not local.
    const coordinate &x = IndexToPosition(i);
    *rho = weight*f(x-position) + (1.0-weight)*(*rho);
  }
}

grid_member(void) LoadGaussianToMesh(double weight, double exponent, 
		   const coordinate& center, 
		   double eps, PointFunction& density) const
{
  GaussianFunction<dim> f2(exponent);
  ClampedGaussianFunction<dim> f(exponent, eps/weight);

  LoadFunctionToMesh(weight, f, center, density);
}

grid_member(void) LoadScalar1DFunctionToMesh(PointFunction& density,  Scalar1DFunctionClass& f, 
					     double range,  coordinate& center,double weight) const
{
  RadialFunction<dim> f3D(f);
  LoadFunctionToMesh(weight, f3D, center, density);  
}


grid_member(double) Integrate(const PointFunction& f) const 
{
  double sum = 0;
  typename PointFunction::const_iterator c = f.begin(), end = f.end();

  for(; c != end; c++) sum += *c;
  return sum*dV;
}

grid_member(double)
Integrate(const PointFunction& f, const PointFunction& g) const 
{
  double sum = 0;
  typename PointFunction::const_iterator fc = f.begin(), end = f.end(), gc = g.begin();

  for(;fc != end; fc++) sum += (*fc)*(*gc);  
  return sum*dV;
}

grid_member(double) Integrate(const PointFunction& f, const PointFunction& V, const PointFunction& g) const
{
  double sum = 0;
  typename PointFunction::const_iterator fc = f.begin(), end = f.end(),
                                         Vc = V.begin(), gc = g.begin();

  for(;fc != end; fc++) sum += (*fc)*(*Vc)*(*gc);
  return sum*dV;
}


grid_member(double) LaplaceElement(const size_t i, const size_t j) const {
  double w = 0;
  if(i==j) w = (1<<dim)/dV;
  else {
    int shared_cell = single_shared_cell(i,j);
    if(shared_cell >= 0) w = -1/dx[shared_cell];
  }

  return w; // TODO: Normalization may be suspect. What is the correct normalization?
}

// ******************************************
// Functionality specific to regular grids 
// ******************************************
grid_member(void) PositionToNode(const coordinate& x, size_t *idx) const {
  const coordinate offset((x-leftcorner)/coordinate(dx));
  for(size_t i=0;i<dim;i++) idx[i] = int(offset.x[i]);
}

grid_member(void) IndexToNode(const size_t global_idx, size_t *idx) const {
  size_t remainder = global_idx;
  for(size_t i=0;i<dim;i++){
    idx[i] = remainder % npts[i];
    remainder /= npts[i];
  }
}

grid_member(size_t) NodeToIndex(const size_t *idx) const {
  size_t global_idx = 0;
  for(size_t i=0;i<dim;i++){
    global_idx += idx[dim-i-1];
    if(i+1<dim) global_idx *= npts[i];
  }
  return global_idx;
}

grid_member(int) single_shared_cell(const size_t global_i, const size_t global_j) const {
  size_t idx_i[dim], idx_j[dim], shared = 0;
  IndexToNode(global_i,&idx_i[0]);
  IndexToNode(global_j,&idx_j[0]);

  for(size_t i=0;i<dim;i++)
    switch(::abs(idx_i[i] - idx_j[i])){
    case 0: break;
    case 1: ++shared; break;
    default: return -1; 
    }
    
  if(shared == 1)
    for(size_t i=0;i<dim;i++) if(::abs(idx_i[i] - idx_j[i]) == 1) return i;

  return -1;
}

grid_member(void) SolvePoisson(const PointFunction& density, PointFunction& vHartree)
{
}

