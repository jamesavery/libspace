#ifndef LIBDISC_DEALIIMESH_H
# define LIBDISC_DEALIIMESH_H

#include <space/function.h>
#include <space/volumes.h>
#include <space/fem/fespace.h>
#include <space/fem/deal.II/feoperator.h>
#include <math.h>
#include <string>
/* deal.II */
#include <grid/tria.h>
#include <dofs/dof_handler.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/persistent_tria.h>
#include <dofs/dof_accessor.h>
#include <fe/fe_q.h>
#include <dofs/dof_tools.h>
#include <fe/fe_values.h>
#include <fe/mapping_q.h>
#include <base/quadrature_lib.h>
#include <base/function.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>
/* /deal.II */


namespace dealii {

  template <class FESpace> class PointWrap;
  template <class FESpace, typename Q=double> class ScalarFunctionWrap; 


  template <int dim> class FESpace : public ::FESpace<dim,PointFunction_,FEFunction_,FEOperator_<dim> > {
    typedef unsigned int   uint_t;
    typedef unsigned short ushort_t;
    typedef unsigned char  uchar_t;
  public:
    /* Test. */
    typedef dealii::SparseMatrix<double>    SparseMatrix;
    typedef dealii::SparsityPattern         SparsityPattern;
    typedef dealii::Vector<double>          cellVector;
    typedef dealii::Vector<double>          dofVector;
    typedef dealii::ConstraintMatrix        ConstraintMatrix;

    /* Imported types from base classes. */
    typedef ::FESpace<dim,PointFunction_,FEFunction_,FEOperator_<dim> > BaseType;
    typedef FESpace<dim>                                             SelfType;

    typedef typename BaseType::coordinate     coordinate;
    typedef typename BaseType::ScalarFunction ScalarFunction;
    typedef dealii::ScalarFunctionWrap<FESpace> ScalarFunctionWrap;
    typedef typename BaseType::PointFunction  PointFunction;
    typedef typename BaseType::FEFunction     FEFunction;
    typedef typename BaseType::FEOperator     FEOperator;
    typedef FEFunction                        DiscreteFunction;
    using BaseType::fe_weights;

    /* PointFunction interface is inherited from FESpace. */
    using BaseType::point_weights;
    using BaseType::point_positions;

    typedef typename BaseType::PointFunctional PointFunctional;
    typedef enum {DIRICHLET, NEUMANN,PERIODIC,MULTIPOLE} BoundaryType;


    double Integrate(const PointFunction& f) const { return BaseType::Integrate(f); }
    double Integrate(const PointFunction& f, const PointFunction& g) const { 
      return BaseType::Integrate(f,g); 
    }
    double Integrate(const PointFunction& f, const PointFunction& V, const PointFunction& g) const { 
      return BaseType::Integrate(f,V,g); 
    }

    void LoadFunctionToMesh(const ScalarFunction& f, PointFunction& density) const {
      return BaseType::LoadFunctionToMesh(f,density);
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
    void LoadFunctionToMesh(const ScalarFunction& f, FEFunction& density) const;
    void LoadGaussianToMesh(double weight, double exponent, const coordinate& position, 
			    double eps, FEFunction& density) const;
    void LoadScalar1DFunctionToMesh(FEFunction& density,  const Scalar1DFunctionClass& f, 
				    double range,  const coordinate& cartesianPosition,double weight=1.0) const;

    /* Integration and inner products */
    double Integrate(const FEFunction& fvalues) const;
    double Integrate(const FEFunction& fvalues, const FEFunction& gvalues) const;
    double Integrate(const FEFunction& fvalues, const FEFunction& Vvalues, const FEFunction& gvalues) const;

    double Integrate(const FEFunction& f, const FEOperator& V, const FEFunction& g) const;

    double LaplaceElement(const size_t i, const size_t j) const;
    /* Solving PDE's */
    void SolvePoisson(const FEFunction& density, FEFunction& result) const;
    void SolvePoisson(const PointFunction& density, FEFunction& result) const;
    void SolvePoisson(const PointFunction& density, PointFunction& result) const; 

    double Value (const FEFunction& f, const coordinate& x) const;

    void set_boundary(unsigned char b_id, const BoundaryType b_type, const double b_value);
    void set_fixed_region(unsigned char material_id, const double potential);
    void set_dielectric_region(unsigned char material_id, const double value);

    void set_fixed_regions(const std::vector<ConstantVolume<dim> >& regions);
    void set_dielectric_regions(const std::vector<ConstantVolume<dim> >& regions);
    
    /* <messy> */
    void ConstructPointFunction(const FEFunction& f, PointFunction& fp) const;
    void ConstructFEFunction(const PointFunction& fp, FEFunction& f) const;
    /* Point-wise point-management -- move to fespace.h/fespace_pointwise.cc */

    uint_t       n_q_pts;	/**< Number of quadrature points per cell. */
    uint_t       n_cells;	/**< Number of cells. */
    uint_t       n_cell_dof;	/**< Number of DOF per cell. */
    uint_t       n_dofs;	/**< Total number of dofs. */
    uint_t       quadrature_order;

    std::map<uint_t,double> dielectric_material, fixed_material; /* Fixed- and dielectric regions defined by material id  */
    VolumeFunction<dim> dielectric_regions; /* Dielectric regions defined by spatial volumes*/
    VolumeFunction<dim> fixed_regions;      /* Fixed regions defined by spatial volumes */

    std::map<uint_t,double> boundary_values; /* Interpolated Dirichlet boundary values */
    std::map<uint_t,double> fixed_dof;	     /* Fixed degrees of freedom -- "internal boundary condition"  */
    std::set<unsigned char> homogeneous_neumann_boundaries;
    typename FunctionMap<dim>::type dirichlet_boundaries, neumann_boundaries;
    Vector<double> neumann_rhs;
    //  VectorValues nodepositions;	/**< Spatial coordinates of each node. Is this necessary? */
    /* </messy> */

    /* Functionality specific to Deal.II-meshes */
    FESpace(const std::string meshfile,uint_t fe_order = 1, uint_t gauss_order=2, const string restartfile = "");
    FESpace(const uint_t npts[dim], const coordinate& leftcorner, const coordinate& dimensions, 
	    uint_t fe_order = 1, uint_t gauss_order=2,bool colorize=false, const string restartfile = ""); 

    FESpace(const uint_t npts[dim], const double cell[dim*dim],
	    uint_t fe_order = 1, uint_t gauss_order=2,bool colorize=false, const string restartfile = "");

    void absolute_error_estimate(const FEFunction& fe_function, const ScalarFunction& function, 
				 cellVector& error/*[n_active_cells()]*/) const;

    void increase_quadrature_order(const int n);
    void refine_grid(const cellVector& estimated_error_per_cell,bool update_at_end=true);
    void refine_grid(const Vector<float>& estimated_error_per_cell,bool update_at_end=true);
    void refine_to_density(const ScalarFunction& density, const double dE,bool update_at_end=true);
    void refine_to_density(PointFunctional& density, const double dE,bool update_at_end=true);
    void refine_to_density_curvature(PointFunctional& density_fe, 
				     const double dC, bool update_at_end=true);
    void refine_max_cell_size(const double max_diameter, const bool update_at_end=true);
    void refine_around_points(const std::vector<coordinate> xs, const double max_diameter=INFINITY, 
			      const double near_diameter=0,const bool update_at_end=true);
    void refine_around_regions(const VolumeFunction<dim>& V, const double max_diameter, 
			       const bool update_at_end=true);
    void refine_grid(const size_t n,bool update_at_end=true);
    void refine_to_regularize(double max_aspect_ratio=2.0);

    /* Output -- perhaps move to a separate class.  */
    void write_mesh(const std::string& path) const;
    void write_function(const std::string& path, const FEFunction& f, 
			const std::string name = std::string("f")) const;
    void write_dof_sparsity(const string& path) const;

    /* Internal stuff. */
    void get_positions();
    void update();
    void update_mass_matrix();
    void update_laplace_matrix();
    void update_boundary_conditions();
    void update_hanging_nodes();
    ConstraintMatrix               hanging_node_constraints;
    SparsityPattern                sparsity_pattern;
    Triangulation<dim>             coarse_grid;
    PersistentTriangulation<dim>   triangulation;
    FE_Q<dim>              fe;
    DoFHandler<dim>        dof_handler;

    QGauss<dim>  quadrature_formula;
    FEValues<dim> fe_values;

    SparseMatrix system_matrix;
    SparseMatrix overlap_matrix;
    SparseMatrix laplace_matrix;

    PreconditionSSOR<> *prec;
  private:
    void update_neumann_boundary(); /* Needs separate function because of special dim=1 case */
  };

  template <class FESpace> class PointWrap : public Point<FESpace::dim> {
  public:
    PointWrap(const typename FESpace::coordinate& x){
      for(size_t i=0;i<FESpace::dim;i++) (*this)[i] = x.x[i];
    }
  };

  template <class FESpace> class WrapPoint : public FESpace::coordinate {
  public:
    using FESpace::coordinate::x;

    WrapPoint(const typename dealii::Point<FESpace::dim>& X){
      for(size_t i=0;i<FESpace::dim;i++) x[i] = X[i];
    }
  };

  template <class FESpace, typename Q> class ScalarFunctionWrap : public Function<FESpace::dim> {
  public:
    static const int dim = FESpace::dim;
    const ::Function<FESpace::dim,Q>& f;
    const typename FESpace::coordinate *center;
    
  ScalarFunctionWrap(const ::Function<FESpace::dim,Q>& f, const typename FESpace::coordinate *center = NULL) 
    : f(f), center(center) { }

    double value(const Point<dim> &x,const unsigned int component = 0) const {
      double x_[dim];
      for(size_t i=0;i<dim;i++) x_[i] = x[i];
	
      typename FESpace::coordinate xp(x_);
      if(center != NULL) xp -= *center;
      return f(xp);
    }
  };

  template <int dim> class FEInterpolationFunction : public ::Function<dim,double> {
  public:
    typedef typename FESpace<dim>::FEFunction FEFunction;

    FEInterpolationFunction(const FEFunction& f, const FESpace<dim>& space){
      /* For each cell c
       *  fe_values.reinit(c) 
       *  fe_values.get_function_values(f,values)
       *  fe_values.get_quadrature_points(points)
       *  fip = interpolationfunction(points,values);
       *  ...noget med f.eks. BSP-trae. 
       */
      /* Alternativ: 
       * allvalues = space.ConstructPointFunction(f);
       * allpoints = space.point_positions
       * fip = interpolationfunction(allpoints,allvalues) -- brug eksternt lib. 
       */
    }
  };

}



#endif
