#ifndef LIBDISC_FUNCTION_H
# define LIBDISC_FUNCTION_H

#include <math.h>
#include <space/storage/storage.h>
#include <space/volumes.h>

template <int dim, typename Q> class Function {
 public:
  typedef SmallVector<dim,double> coordinate;

  virtual Q operator()(const coordinate& x, const off_t component=0) const = 0;
  virtual ~Function() {}
};

template <typename Q> class Function1D {
 public:
  virtual Q operator()(const double x, const off_t component=0) const = 0;
  virtual ~Function1D() {}
};

typedef Function1D<double> Scalar1DFunctionClass;
typedef Function<2,double> Scalar2DFunctionClass;
typedef Function<3,double> Scalar3DFunctionClass;


template <int dim, typename Q> class ConstantFunction : public Function<dim,Q> {
 public:
  typedef typename Function<dim,Q>::coordinate coordinate;

  const Q c;
  ConstantFunction(const Q& c) : c(c) {}
  
  Q operator()(const coordinate& x, const off_t component=0) const {
    return c;
  }
};

template <int dim, typename Q> class SuperPositionFunction: public Function<dim,Q> {
 public:
  typedef typename Function<dim,Q>::coordinate coordinate;

  std::vector<double> coefficients;
  std::vector<coordinate> centers;
  std::vector< const Function<dim,Q> const * > funs;

  SuperPositionFunction(){}

  void add_function(const Function<dim,Q> const *fun){
    coefficients.push_back(1.0);
    centers.push_back(coordinate());
    funs.push_back(fun);
  }

  void add_function(const double coefficient, const coordinate& center, const Function<dim,Q> const *fun){
    coefficients.push_back(coefficient);
    centers.push_back(center);
    funs.push_back(fun);
  }

  Q operator()(const coordinate& x, const off_t component=0) const {
    const size_t m = funs.size();
    double value = 0;
    for(size_t i=0;i<m;i++){
      const coordinate &center = centers[i];
      value += coefficients[i]*(*funs[i])(x - center);
    }
    return value;
  }
  
  ~SuperPositionFunction(){
    for(size_t i=0;i<funs.size();i++)
      delete funs[i];
    funs.resize(0);
  }
};

/* MISSING: Normalization constant */
template <int dim> class GaussianFunction : public Function<dim,double> {
public:
  typedef typename Function<dim,double>::coordinate coordinate;

  double     exponent;
  GaussianFunction(double a) : exponent(a) { }
 
  double operator()(const coordinate& x, const off_t component=0) const
  {
    return exp(-exponent*x.dot(x));
  }
};

/* MISSING: Normalization constant */
template <int dim> class ClampedGaussianFunction : public Function<dim,double> {
public:
  typedef typename Function<dim,double>::coordinate coordinate;

  double exponent,rCutSqr;
  ClampedGaussianFunction(double a, double eps) : exponent(a)
  {
    double rCut=sqrt(-log(eps)/exponent);
    rCutSqr = rCut*rCut;
  }
 
  double operator()(const coordinate& x, const off_t component=0) const
  {
    double dsqr = x.dot(x);
    return dsqr < rCutSqr? exp(-exponent*dsqr) : 0;
  }
};

template <int dim> class RadialFunction : public Function<dim,double> {
public:
  typedef typename Function<dim,double>::coordinate coordinate;

  const Scalar1DFunctionClass& f;
  double range;

  RadialFunction(const Scalar1DFunctionClass& f, double range=INFINITY) 
    : f(f), range(range) {}
 
  double operator()(const coordinate& x, const off_t component=0) const
  {
    double r = sqrt(x.dot(x)); // Can be made a lot faster, if f is function of r^2 instead.
    return r<=range? f(r) : 0;	
  }
};

template <int dim> class VolumeFunction : public Function<dim,double> {
 public:
  double dfault;
  typedef SmallVector<dim> coordinate;
  std::vector<ConstantVolume<dim> > volumes;

 VolumeFunction(){  }
 VolumeFunction(const std::vector<ConstantVolume<dim> >& volumes, const double dfault=1.0) : volumes(volumes), dfault(dfault) {  }

  double operator()(const coordinate& x, const off_t component=0) const {
    int i = volume_number(x);
    if(i>=0 /* && !volumes[i].fixed_potential*/){
      return volumes[i].value;
    } else 			/* Default is vacuum */
      return dfault;
  }

  int volume_number(const coordinate& x) const {
    for(size_t i=0;i<volumes.size();i++)
      if(volumes[i].shape != NULL && volumes[i].shape->point_inside_volume(x)) return i;
    return -1;
  }

  size_t closest_volume_number(const coordinate& x, double *distance_return = NULL) const {
    int minnum = volume_number(x);
    double mindist = INFINITY;
    if(minnum>=0) mindist = 0.0;
    else 
      for(size_t i=0;i<volumes.size();i++){
	const double dist = volumes[i].shape->approximate_distance(x);
	if(dist < mindist){
	  mindist = dist;
	  minnum = i;
	}
      }
    if(distance_return != NULL) *distance_return = mindist;
    return minnum;
  }

  friend std::ostream& operator<<(std::ostream& S, const VolumeFunction& v){
    S << "VolumeFunction with default value " << v.dfault << "." << endl;
    for(size_t i=0;i<v.volumes.size();i++) S << "\t" << v.volumes[i]  << endl;
    return S;
  } 
};


#endif
