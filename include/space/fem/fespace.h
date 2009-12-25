#ifndef LIBDISC_FESPACE_H
# define LIBDISC_FESPACE_H

#include <space/discretization.h>

template <int dim_, class PointFunction_, class FEFunction_, class FEOperator_> 
  class FESpace : public Discretization<dim_,PointFunction_> /* FEM Mesh type 1 */, 
                  public Discretization<dim_,FEFunction_,FEOperator_> /* FEM Mesh type 2 */ { 
 public:
  typedef typename Discretization<dim_,FEFunction_,FEOperator_>::coordinate     coordinate;
  typedef typename Discretization<dim_,FEFunction_,FEOperator_>::ScalarFunction ScalarFunction;
  typedef PointFunction_ PointFunction;
  typedef FEFunction_ FEFunction;
  typedef FEOperator_ FEOperator;




  static const int dim = dim_;

  /* PointFunction-interface. TODO: Det er lidt rodet med PointFunction og FEFunction sammen. Split? */
  /* Adding scalar functions to a mesh density */
  void LoadFunctionToMesh(const ScalarFunction& f, PointFunction& density) const;
  void LoadFunctionToMesh(double weight, const ScalarFunction& f, 
			  const coordinate& center, PointFunction& density) const;
  void LoadGaussianToMesh(double weight, double exponent, const coordinate& center, 
			  double eps, PointFunction& density) const;
  void LoadScalar1DFunctionToMesh(PointFunction& density,  const Scalar1DFunctionClass& f, 
				  double range,  const coordinate& center,double weight=1.0) const;


  /* Integration and inner products */
  double Integrate(const PointFunction& f) const;
  double Integrate(const PointFunction& f, const PointFunction& g) const;
  double Integrate(const PointFunction& f, const PointFunction& V, const PointFunction& g) const;
  /* Solving PDE's. Library-dependent. */
  virtual void SolvePoisson(const PointFunction& density, PointFunction& vHartree) const = 0;


  /* FEFunction-interface. Library-dependent. */
  void LoadFunctionToMesh(const ScalarFunction& f, FEFunction& density) const;
  virtual void LoadFunctionToMesh(double weight, const ScalarFunction& f, 
			  const coordinate& position, FEFunction& density) const = 0;
  virtual void LoadGaussianToMesh(double weight, double exponent, const coordinate& position, 
			  double eps, FEFunction& density) const = 0;
  virtual void LoadScalar1DFunctionToMesh(FEFunction& density,  const Scalar1DFunctionClass& f, 
				  double range,  const coordinate& cartesianPosition,double weight=1.0) const = 0;

  virtual double Integrate(const FEFunction& f) const = 0;
  virtual double Integrate(const FEFunction& f, const FEFunction& g) const = 0;
  virtual double Integrate(const FEFunction& f, const FEOperator& V, const FEFunction& g) const = 0;

  virtual double LaplaceElement(const size_t i, const size_t j) const = 0;
  virtual void SolvePoisson(const FEFunction& density,    FEFunction& vHartree) const = 0;


  /* Point-wise auxiliary data. */
  std::vector<coordinate> point_positions;      /* Compact */
  PointFunction           point_weights;        /* Compact   point-function */

  /* FE-wise auxiliary data. */
  FEFunction              fe_weights; 

  /* Output -- perhaps move to a separate class.  */
  virtual void write_mesh(const std::string& path) const = 0;
  virtual void write_function(const std::string& path, const FEFunction& f) const = 0;
  virtual void write_dof_sparsity(const string& path) const = 0;

  /* Useful stuff for refinement */
  class PointFunctional {
  public:
    virtual PointFunction& fill(const std::vector<coordinate>& positions) = 0;
  };
};

#endif
