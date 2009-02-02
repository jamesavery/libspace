#include <libspace/fem/fespace.h>

#define fespace_member(returntype) \
template <int dim, class PointFunction, class FEF, class FEO> returntype FESpace<dim,PointFunction,FEF,FEO>::


fespace_member(double)
Integrate(const PointFunction& f) const 
{
  double sum = 0;
  typename PointFunction::const_iterator c = f.begin(), end = f.end();
  typename PointFunction::const_iterator weight = point_weights.begin();
  
  for(; c != end; c++,weight++) sum += (*c)*(*weight);
  
  return sum;
}

fespace_member(double)
Integrate(const PointFunction& f,const PointFunction& g) const 
{
  double sum = 0;
  typename PointFunction::const_iterator cf = f.begin(), cg = g.begin(), end = f.end();
  typename PointFunction::const_iterator weight = point_weights.begin();
  
  for(; cf != end; cf++,cg++,weight++) sum += (*cf)*(*cg)*(*weight);
  
  return sum;
}

fespace_member(double)
Integrate(const PointFunction& f,const PointFunction& V, const PointFunction& g) const 
{
  double sum = 0;
  typename PointFunction::const_iterator cf = f.begin(), cg = g.begin(), 
                                         cV = V.begin(), end = f.end();
  typename PointFunction::const_iterator weight = point_weights.begin();
  
  for(; cf != end; cf++,cg++,cV++,weight++) sum += (*cf)*(*cV)*(*cg)*(*weight);
  
  return sum;
}


fespace_member(void)
LoadFunctionToMesh(double weight, const ScalarFunction& f, 
			const coordinate& center, PointFunction& density) const 
{
  density.resize(point_weights.size());
  typename std::vector<coordinate>::const_iterator position = point_positions.begin();
  typename PointFunction::iterator rho = density.begin(), end = density.end();

  for(;rho!=end;rho++,position++){
    const coordinate &x = *position;
    (*rho) = weight*f(x-center) + (1.0-weight)*(*rho);
  }
}

fespace_member(void)
LoadGaussianToMesh(double weight, double exponent, const coordinate& center, 
			double eps, PointFunction& density) const 
{
  ClampedGaussianFunction<dim> f(exponent, eps/weight);
  LoadFunctionToMesh(weight, f, center, density);  
}

fespace_member(void)
LoadScalar1DFunctionToMesh(PointFunction& density,  Scalar1DFunctionClass& f, 
				double range,  coordinate& center,double weight) const 
{
  RadialFunction<dim> f3D(f);
  LoadFunctionToMesh(weight, f3D, center, density);    
}

#undef fespace_member
