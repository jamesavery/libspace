#ifndef LIBDISC_LIBMESH_H
# define LIBDISC_LIBMESH_H

#include <libdiscretization/function.h>
#include <libdiscretization/fem/fespace.h>
#include <libdiscretization/fem/deal.II/feoperator.h>


#include <libmesh/libmesh.h>
#include <libmesh/mesh.h>
  


namespace libMesh {

  template <class FESpace> class PointWrap;
  template <class FESpace, typename Q=double> class ScalarFunctionWrap; 

  template <int dim> class FESpace : public ::FESpace<dim,PointFunction_,FEFunction_,FEOperator<dim> > {
  public:
    /* Test. */
    //    typedef dealii::TrilinosWrappers::SparseMatrix    SparseMatrix;
    //    typedef dealii::TrilinosWrappers::SparsityPattern SparsityPattern;
    typedef dealii::SparseMatrix<double>    SparseMatrix;
    typedef dealii::SparsityPattern         SparsityPattern;
    typedef dealii::ConstraintMatrix        ConstraintMatrix;

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
/*     void LoadScalar1DFunctionToMesh(PointFunction& density,  Scalar1DFunctionClass& f,  */
/* 				    double range,  coordinate& cartesianPosition,double weight=1.0) const */
/*     { */
/*       return BaseType::LoadScalar1DFunctionToMesh(density,f,range,weight); */
/*     } */

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
    quadrature_order(gauss_order), fe(fe_order), dof_handler(triangulation), 
      quadrature_formula(gauss_order), 
      fe_values(fe, quadrature_formula, 
		dealii::update_values|dealii::update_JxW_values|dealii::update_quadrature_points|dealii::update_gradients)
    {
      PointWrap<FESpace> p1(leftcorner), p2(leftcorner+dimensions);
      std::vector<size_t> npts(dim); for(size_t i=0;i<dim;i++) npts[i] = npts_[i];

      printf("fe_order    = %d\n"
	     "gauss_order = %d\n", fe_order, gauss_order);

      dealii::GridGenerator::subdivided_hyper_rectangle(triangulation, npts,p1,p2);
      update();  
    }

    void absolute_error_estimate(const FEFunction& fe_function, const ScalarFunction& function, 
				 dealii::Vector<double>& error/*[n_active_cells()]*/) const;

    void refine_grid(const dealii::Vector<double>& estimated_error_per_cell);
    void refine_grid(const dealii::Vector<float>& estimated_error_per_cell);

    /* Output -- perhaps move to a separate class.  */
    void write_mesh(const std::string& path) const;
    void write_function(const std::string& path, const FEFunction& f) const;
    void write_dof_sparsity(const string& path) const;

    /* Internal stuff. */
    void get_positions();
    void update();
    ConstraintMatrix               hanging_node_constraints;
    SparsityPattern                sparsity_pattern;
    dealii::Triangulation<dim>     triangulation;
    dealii::FE_Q<dim>              fe;
    dealii::DoFHandler<dim>        dof_handler;

    dealii::QGauss<dim>  quadrature_formula;
    dealii::FEValues<dim> fe_values;

    SparseMatrix system_matrix;
    SparseMatrix overlap_matrix;
    SparseMatrix laplace_matrix;
  };

  template <class FESpace> class PointWrap : public dealii::Point<FESpace::dim> {
  public:
    PointWrap(const typename FESpace::coordinate& x){
      for(size_t i=0;i<FESpace::dim;i++) (*this)[i] = x.x[i];
    }
  };

  template <class FESpace, typename Q> class ScalarFunctionWrap : public dealii::Function<FESpace::dim> {
  public:
    static const int dim = FESpace::dim;
    const Function<FESpace::dim,Q>& f;
    const typename FESpace::coordinate *center;
    
  ScalarFunctionWrap(const Function<FESpace::dim,Q>& f, const typename FESpace::coordinate *center = NULL) 
    : f(f), center(center) { }

    double value(const dealii::Point<dim> &x,const unsigned int component = 0) const {
      double x_[dim];
      for(size_t i=0;i<dim;i++) x_[i] = x[i];
	
      typename FESpace::coordinate xp(x_);
      if(center != NULL) xp -= *center;
      return f(xp);
    }
  };
}



#endif
