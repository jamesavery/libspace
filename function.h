#ifndef LIBDISC_FUNCTION_H
# define LIBDISC_FUNCTION_H

#include <math.h>
#include <libdiscretization/storage/storage.h>

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
  
  double operator()(const coordinate& x, const off_t component=0) const {
    return c;
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



#endif
