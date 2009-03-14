/* \file discretization.h
   \brief Abstract space discretization class. 
 */
#ifndef LIBDISC_DISCRETIZATION_H
# define LIBDISC_DISCRETIZATION_H

#include <space/libspace-config.h>
#include <space/storage/storage.h>
#include <space/function.h>

template <int dim, class PF = PointFunction_Simple<>, class PO = PF >
  class Discretization {	/* Specialized to double -- perhaps templatize field instead. */
  
 public:
  typedef SmallVector<dim,double> coordinate;
  typedef Function<dim,double> ScalarFunction;
  typedef PF DiscreteFunction;

  /* Adding scalar functions to a mesh density */
  virtual void LoadFunctionToMesh(double weight, const ScalarFunction& f, 
				  const coordinate& center, PF& density) const = 0;
  virtual void LoadGaussianToMesh(double weight, double exponent, const coordinate& center, 
				  double eps, PF& density) const = 0;
  virtual void LoadScalar1DFunctionToMesh(PF& density,  Scalar1DFunctionClass& f, 
				  double range,  coordinate& center,double weight=1.0) const = 0;

  /* Integration and inner products */
  virtual double Integrate(const PF& f) const = 0;
  virtual double Integrate(const PF& f, const PF& g) const = 0;
  virtual double Integrate(const PF& f, const PO& V, const PF& g) const = 0;

  virtual double LaplaceElement(const size_t i, const size_t j) const = 0;
  /* Solving PDE's */
  virtual void SolvePoisson(const PF& density, PF& vHartree) = 0;

 private:
};

#endif
