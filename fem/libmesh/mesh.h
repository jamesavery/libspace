b#ifndef LIBDISC_LIBMESH_H
# define LIBDISC_LIBMESH_H

#include <libspace/function.h>
#include <libspace/fem/fespace.h>
#include <libspace/fem/libmesh/feoperator.h>

#include <libmesh/libmesh.h>
#include <libmesh/mesh.h>

namespace libMesh {

  template <class FESpace> class PointWrap;
  template <class FESpace, typename Q=double> class ScalarFunctionWrap; 

  template <int dim> class FESpace : public ::FESpace<dim,PointFunction_,FEFunction_,FEOperator<dim> > {
  public:
    /* Test. */
    //    typedef TrilinosWrappers::SparseMatrix    SparseMatrix;
    //    typedef TrilinosWrappers::SparsityPattern SparsityPattern;
    typedef SparseMatrix<double>    SparseMatrix;
    typedef SparsityPattern         SparsityPattern;
    typedef ConstraintMatrix        ConstraintMatrix;

    /* Imported types from base classes. */
    typedef ::FESpace<dim,PointFunction_,FEFunction_,FEOperator<dim> > BaseType;
    typedef FESpace<dim>                                             SelfType;

    typedef typename BaseType::coordinate     coordinate;
    typedef typename BaseType::ScalarFunction ScalarFunction;
    typedef ScalarFunctionWrap<FESpace> ScalarFunctionWrap;
    typedef typename BaseType::PointFunction  PointFunction;
    typedef typename BaseType::FEFunction     FEFunction;
    typedef typename BaseType::FEOperator     FEOperator;
    typedef FEFunction                        DiscreteFunction;
    using BaseType::fe_weights;

    /* PointFunction interface is inherited from FESpace. */
    using BaseType::point_weights;
    using BaseType::point_positions;

    double Integrate(const PointFunction& f) const { return BaseType::Integrate(f); }
    double Integrate(const PointFunction& f, const PointFunction& g) const { 
      return BaseType::Integrate(f,g); 
    }
    double Integrate(const PointFunction& f, const PointFunction& V, const PointFunction& g) const { 
      return BaseType::Integrate(f,V,g); 
    }

    void LoadFunctionToMesh(double weight, const ScalarFunction& f, 
			    const coordinate& center, PointFunction& density) const {
      return BaseType::LoadFunctionToMesh(weight,f,center,density);
    }
    void LoadGaussianToMesh(double weight, double exponent, const coordinate& center, 
			    double eps, PointFunction& density) const {
      return BaseType::LoadGaussianToMesh(weight,exponent,center,eps,density);
    }

    /* FEFunction interface */
    void LoadFunctionToMesh(double weight, const ScalarFunction& f, 
			    const coordinate& position, FEFunction& density) const;
    void LoadGaussianToMesh(double weight, double exponent, const coordinate& position, 
			    double eps, FEFunction& density) const;
    void LoadScalar1DFunctionToMesh(FEFunction& density,  Scalar1DFunctionClass& f, 
				    double range,  coordinate& cartesianPosition,double weight=1.0) const;

    /* Integration and inner products */
    double Integrate(const FEFunction& fvalues) const;
    double Integrate(const FEFunction& fvalues, const FEFunction& gvalues) const;
    double Integrate(const FEFunction& fvalues, const FEFunction& Vvalues, const FEFunction& gvalues) const;

    double Integrate(const FEFunction& f, const FEOperator& V, const FEFunction& g) const;

    double LaplaceElement(const size_t i, const size_t j) const;
    /* Solving PDE's */
    void SolvePoisson(const FEFunction& density, FEFunction& result);
    void SolvePoisson(const PointFunction& density, PointFunction& result);

    /* <messy> */
    void   ConstructPointFunction(const FEFunction& f, PointFunction& fp) const;

    size_t       n_q_pts;	/**< Number of quadrature points per cell. */
    size_t       n_cells;	/**< Number of quadrature points per cell. */
    size_t       n_dof;		/**< Number of DOF per cell. */
    size_t       n_weights;	/**< Total number of weights. */
    size_t       quadrature_order;

    //  VectorValues nodepositions;	/**< Spatial coordinates of each node. Is this necessary? */
    /* </messy> */

    /* Functionality specific to Deal.II-meshes */
    FESpace(const size_t npts_[dim], const coordinate& leftcorner, const coordinate& dimensions, 
	    size_t fe_order = 1, size_t gauss_order=2) :
    {
    }

    void absolute_error_estimate(const FEFunction& fe_function, const ScalarFunction& function, 
				 NumericVector<double>& error/*[n_active_cells()]*/) const;

    void refine_grid(const NumericVector<double>& estimated_error_per_cell);
    void refine_grid(const NumericVector<float>& estimated_error_per_cell);

    /* Output -- perhaps move to a separate class.  */
    void write_mesh(const std::string& path) const;
    void write_function(const std::string& path, const FEFunction& f) const;
    void write_dof_sparsity(const string& path) const;

    /* Internal stuff. */
    void get_positions();
    void update();

    QBase                  quadrature_formula;
    DofConstraints         dof_constraints;
    
    SparsityPattern        sparsity_pattern;
    Triangulation<dim>     triangulation;
    FE_Q<dim>              fe;
    DoFHandler<dim>        dof_handler;

    SparseMatrix system_matrix;
    SparseMatrix overlap_matrix;
    SparseMatrix laplace_matrix;
  };

  template <class FESpace> class PointWrap : public Point {
  public:
    PointWrap(const typename FESpace::coordinate& x){
      for(size_t i=0;i<FESpace::dim;i++) (*this)(i) = x.x[i];
    }
  };

  template <class FESpace> class ScalarFunctionWrap : public BaseFunction {
  public:
    static const int dim = FESpace::dim;
    const Function<FESpace::dim,Q>& f;
    const typename FESpace::coordinate *center;
    
  ScalarFunctionWrap(const Function<FESpace::dim,double>& f, const typename FESpace::coordinate *center = NULL) 
    : f(f), center(center) { }

    Number value(const Point &x) const {
      double x_[dim];
      for(size_t i=0;i<dim;i++) x_[i] = x(i);
	
      typename FESpace::coordinate xp(x_);
      if(center != NULL) xp -= *center;
      return f(xp);
    }
  };
}



#endif
