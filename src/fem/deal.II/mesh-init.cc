// \file mesh-init.cc This file contains the implementation of FESpace's constructors.

#include <space/fem/deal.II/mesh.h>
#include <space/volumes.h>

#include <stdio.h>
#include <iostream>
#include <fstream>

#  include <grid/tria.h>
#  include <grid/grid_generator.h>

#  include <grid/grid_out.h>
#  include <numerics/data_out.h>
#  include <grid/grid_refinement.h>
#  include <grid/grid_tools.h>
#  include <grid/grid_in.h>
#  include <dofs/dof_renumbering.h>

#  include <fe/fe_tools.h>
#  include <lac/solver_selector.h>



#define fespace_member(returntype) template <int dim> returntype FESpace<dim>::
namespace dealii {

  template <int dim> class lineartransform {
  public:
    double A [dim*dim];
    double lowerleft[dim];

    lineartransform(const double A_[dim*dim]) { 
      memcpy(A,A_,sizeof(A)); 
    }
    Point<dim> operator () (const Point<dim>& x) const {
      Point<dim> y;
      for(size_t i=0;i<dim;i++){
	double yi = 0;
	for(size_t j=0;j<dim;j++) yi += A[j*dim+i]*x(j);
	y(i) = yi;
      }
      return y;
    }
  };


/// Construct a FESpace from a GMSH2 .msh-file
  fespace_member() FESpace(const string meshfile,uint_t fe_order , uint_t gauss_order, 
			   const string restartfile):
    quadrature_order(gauss_order), fe(fe_order), dof_handler(triangulation), 
      quadrature_formula(gauss_order), 
  fe_values(fe, quadrature_formula, update_values|update_JxW_values|update_quadrature_points|update_gradients),
  triangulation(coarse_grid)
 {
    GridIn<dim> gridreader;
    gridreader.attach_triangulation(coarse_grid);
    ifstream input_file(meshfile.c_str());
    try {
      gridreader.read_msh(input_file);
    } catch(dealii::ExceptionBase e){
      cerr << "Exception data:\n";
      e.print_exc_data(cerr);
      cerr << "Exception info:\n";
      e.print_info(cerr);
      cerr << "Stack trace:\n";
      e.print_stack_trace(cerr);
      abort();
    }

    cout << "Mesh diameter of coarse grid: " << GridTools::diameter(coarse_grid) << endl;
    // TODO: Reorder DoFs to improve preconditioner performance.
    prec = new PreconditionSSOR<>();

    if(!restartfile.empty()){
      ifstream f(restartfile.c_str());
      if(f.is_open()){
	triangulation.read_flags(f);
	f.close();
      } else {
	fprintf(stderr,"Can't open %s for reading.\n",restartfile.c_str());
      }
    }
    triangulation.restore();
    update();
    cout << "Mesh diameter of refined grid: " << GridTools::diameter(triangulation) << endl;
    
    cerr << "Original cells have material ids: [";
    for(typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin(); cell != dof_handler.end(); cell++)
      fprintf(stderr," %d",cell->material_id());
    cerr << " ]" << endl;
  }

/// Construct a rectangular FESpace mesh 
  fespace_member() FESpace(const uint_t npts_[dim], const coordinate& leftcorner, 
			   const coordinate& dimensions, uint_t fe_order, uint_t gauss_order, 
			   bool colorize, const string restartfile) :
  triangulation(coarse_grid),  
    quadrature_order(gauss_order), fe(fe_order), dof_handler(triangulation), 
      quadrature_formula(gauss_order), fe_values(fe, quadrature_formula, 
      update_values|update_JxW_values|update_quadrature_points|update_gradients)
    {
      PointWrap<FESpace> p1(leftcorner), p2(leftcorner+dimensions);
      std::vector<uint_t> npts(dim); for(size_t i=0;i<dim;i++) npts[i] = npts_[i];

      GridGenerator::subdivided_hyper_rectangle(coarse_grid, npts,p1,p2,colorize);

      if(!restartfile.empty()){
	ifstream f(restartfile.c_str());
	triangulation.read_flags(f);
	f.close();
      }
      triangulation.restore();

      update();  
    }

/// Construct a rhombic FESpace mesh with lower left corner (0,0,0)
/// and cartesian coordinate system with "unit" vectors given by cell[0],...,cell[2].
///  Good for cells in a periodic grid structure.
  fespace_member() FESpace(const uint_t npts[dim], const double cell[dim*dim],
			   uint_t fe_order, uint_t gauss_order, bool colorize, const string restartfile) 
  : triangulation(coarse_grid),
    quadrature_order(gauss_order), fe(fe_order), dof_handler(triangulation), quadrature_formula(gauss_order), 
    fe_values(fe, quadrature_formula, update_values|update_JxW_values|update_quadrature_points|update_gradients)
    {
      vector<uint_t> npts_(dim);
      Point<dim> p1, p2;
      for(size_t i=0;i<dim;i++){
	p1(i) =   0;
	p2(i) = 1.0;
	npts_[i] = npts[i];
      }

      GridGenerator::subdivided_hyper_rectangle(coarse_grid, npts_, p1, p2,colorize);

      cout << "Linear transform: " << endl;
      for(size_t i=0;i<3;i++){
	cout << '\t';
	for(size_t j=0;j<3;j++)
	  cout << cell[i*3+j] << " ";
	cout << '\n';
      }

      cout << "Mesh diameter before linear transform: " << GridTools::diameter(coarse_grid) << endl;
      GridTools::transform(lineartransform<dim>(cell),coarse_grid);
      cout << "Mesh diameter after linear transform: " << GridTools::diameter(coarse_grid) << endl;

      if(!restartfile.empty()){
	ifstream f(restartfile.c_str());
	triangulation.read_flags(f);
	f.close();
      }
      triangulation.restore();
      update();
    }

}
#undef fespace_member
