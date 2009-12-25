#include <space/fem/deal.II/mesh.h>
#include <space/volumes.h>

#include <algorithm>
#include <numerics/vectors.h>
#include <numerics/fe_field_function.h>
#include <base/function.h>

#include <stdio.h>
#include <iostream>
#include <fstream>

#  include <grid/tria.h>
#  include <grid/grid_generator.h>
#  include <grid/tria_accessor.h>
#  include <grid/tria_iterator.h>
#  include <grid/tria_boundary_lib.h>

#  include <grid/grid_out.h>
#  include <numerics/data_out.h>
#  include <lac/sparsity_pattern.h>
#  include <grid/grid_refinement.h>
#  include <grid/grid_tools.h>
#  include <grid/grid_in.h>
#  include <dofs/dof_renumbering.h>

#  include <fe/fe_tools.h>
#  include <lac/solver_selector.h>

// Auxiliary stuff. Perhaps move to separate file.
#define fespace_member(returntype) template <int dim> returntype FESpace<dim>::

namespace dealii {
  using namespace std;  

  fespace_member(void) update_hanging_nodes()
  {
    fprintf(stderr,"Updating hanging nodes.\n");

    hanging_node_constraints.clear ();
    DoFTools::make_hanging_node_constraints (dof_handler, hanging_node_constraints);
    hanging_node_constraints.close ();
    hanging_node_constraints.condense (sparsity_pattern);
  }

  fespace_member(void) update_mass_matrix()
  {
    printf("Building overlap matrix.\n");
    dealii::ConstantFunction<dim> one(1.0);
    overlap_matrix.clear();
    overlap_matrix.reinit(sparsity_pattern);
    fe_weights.resize(dof_handler.n_dofs());
   
    MatrixTools::create_mass_matrix(dof_handler,quadrature_formula,overlap_matrix,one,fe_weights.coefficients);    
    
  }
  
  fespace_member(void) set_boundary(unsigned char b_id, const BoundaryType b_type, const double b_value)
  {
    switch(b_type){
    case DIRICHLET: 
      cerr << "Dirichlet boundary " << static_cast<int>(b_id) << " with value " << b_value << endl;
      dirichlet_boundaries[b_id] = new dealii::ConstantFunction<dim>(b_value); // XXX: Memory leak
      break;
    case NEUMANN: 
      if(b_value == 0.0)
	homogeneous_neumann_boundaries.push_back(b_id);
      else 
	neumann_boundaries[b_id] = new dealii::ConstantFunction<dim>(b_value); // XXX: Memory leak

      break;
    case PERIODIC:
      cerr << "Periodic boundary conditions are not implemented for FEM yet." << endl;
      break;
    case MULTIPOLE:
      cerr << "Multipole boundary conditions are not implemented for FEM yet." << endl;
      break;
    }
  }


  fespace_member(void) set_dielectric_region(unsigned char material_id, const double value = 0)
  {
    dielectric_material[material_id] = value;
  }

  fespace_member(void) set_fixed_region(unsigned char material_id, const double potential = 0)
  {
    cerr << "Fixing material " << int(material_id) << " to value " << potential << endl;
    fixed_material[material_id] = potential;
  }

  fespace_member(void) set_fixed_regions(const vector<ConstantVolume<dim> >& regions)
  {
    fixed_regions = VolumeFunction<dim>(regions,-1.0);
    cerr << "Fixed regions: " << fixed_regions << endl;
  }

  fespace_member(void) set_dielectric_regions(const vector<ConstantVolume<dim> >& regions){
    dielectric_regions = VolumeFunction<dim>(regions);
    cerr << "Dielectric regions: " << dielectric_regions << endl;
  }

  fespace_member(void) update_laplace_matrix()
  {

    vector<double> material_values(n_q_pts);
    // and Laplacian matrix  \f$m_{ij} = \int_\Omega \nabla\phi_i(x)\cdot \nabla\phi_j(x) dx\f$
    cerr << "Building Laplace matrix.\n"; // Slow; Should only generate laplace_matrix when asked; 
				          // Perhaps have a PreparePoisson (Also builds preconditioner?)
    laplace_matrix.reinit(sparsity_pattern);

    // TODO: 
    // create_laplace_matrix needs material constant to be a function of space; 
    // What to do when we have it as a function of cell->material_id()?
    {
      ScalarFunctionWrap material_fun(dielectric_regions);
      MatrixTools::create_laplace_matrix(dof_handler,quadrature_formula,laplace_matrix,&material_fun);
    }

    hanging_node_constraints.condense (laplace_matrix);
  }

  fespace_member(void) update_boundary_conditions()
  {
    boundary_values.clear();
    fixed_dof.clear();		// XXX: If material_id is reintroduced, this must be dealt with.

    // Set fixed regions
    {

      typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();   
    
      vector<uint_t> global_dof_indices (n_cell_dof);

      for (size_t cellnum=0; cell!=endc; ++cell,++cellnum){
	const WrapPoint<SelfType> x(cell->center());
	double vol_dist;
	int vol_idx = fixed_regions.closest_volume_number(x,&vol_dist);

	if(vol_dist < cell->diameter()){  
	  // Cell is likely to intersect with a fixed region -- fix the whole bloody cell, man!
	  cerr << "HIGHLY DEPRECATED: Fixed region by volume. "
	           "Don't do this unless you know what you're doing!\n";
	  const double value = fixed_regions.volumes[vol_idx].value; 

	  cell->get_dof_indices (global_dof_indices);
	  
	  for(size_t i=0;i<n_cell_dof;i++)
	    fixed_dof[global_dof_indices[i]] = value;

	} else if(fixed_material.find(cell->material_id()) != fixed_material.end()){
	  const double value = fixed_material[cell->material_id()];

	  cerr << "Material " << int(cell->material_id()) << " fixed to " << value << endl;

	  cell->get_dof_indices (global_dof_indices);
	  
	  for(size_t i=0;i<n_cell_dof;i++)
	    fixed_dof[global_dof_indices[i]] = value; 	  
	}
      }
    }

    // Homogeneous von Neumann boundary condition is a noop -- hom. von Neumann is default behaviour
    cerr << "Homogeneous von Neumann boundaries: "<<endl;
    for(size_t i=0;i<homogeneous_neumann_boundaries.size();i++)
      cerr << "\t" << static_cast<int>(homogeneous_neumann_boundaries[i]) << endl;

    // Apply all Dirichlet boundary conditions
    if(!dirichlet_boundaries.empty()){
      cerr << "Dirichlet boundaries: "<<endl;
      for(typename FunctionMap<dim>::type::const_iterator i=dirichlet_boundaries.begin(); i!=dirichlet_boundaries.end();i++)
	cerr << "\t" << static_cast<int>(i->first) << endl;
      VectorTools::interpolate_boundary_values (MappingQ1<dim>(),dof_handler,dirichlet_boundaries,boundary_values);
    }
    // Merge Dirichlet BCs with fixed regions
    cerr << boundary_values.size() << " fixed dofs from Dirichlet BVC, " << fixed_dof.size()
	 << " fixed dofs from set_fixed_regions()\n";
    for(map<uint_t,double>::const_iterator m=fixed_dof.begin();m!=fixed_dof.end();m++)
      boundary_values[m->first] = m->second;

    // Apply all nonhomogeneous Von Neumann boundary conditions
    if(!neumann_boundaries.empty()){
      cerr << "Houston: We have a non-homogeneous Neumann boundary condition!" << endl;
      update_neumann_boundary();
    }
  }

  template<> void FESpace<1>::update_neumann_boundary(){}
  fespace_member(void) update_neumann_boundary(){
    // TODO: (1) User provides either desired gradient Dphi (such that neumann_value(x) = dot(Dphi(x),normal(x)))
    //       or, as now, simply neumann_value(x).
    // TODO: (2) Cache boundary faces instead of having to iterate through *all* cells whenever one wants to do surface stuff
    //       -- as there are many, many more interior cells than there are boundary faces. 
    // TODO: (3) This is unnecessarily complicated, and the code is repeated many places. Factor this out somehow.

    typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
    typename FunctionMap<dim>::type::const_iterator fun;
    Vector<double>  cell_rhs (n_cell_dof);
    std::vector<double> neumann_values;
    std::vector<unsigned int> global_dof_indices (n_cell_dof);
    QGauss<dim-1> face_quadrature(quadrature_order);

    size_t n_face_q_points = face_quadrature.size();
    neumann_values.resize(n_face_q_points);
    neumann_rhs.resize(n_dofs);

    FEFaceValues<dim> fe_face_values (fe, face_quadrature, 
				      update_values         |  update_quadrature_points  |
				      update_normal_vectors |  update_JxW_values);
    neumann_rhs = 0;
    for(;cell!=endc;cell++){
      for (size_t face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
	if (cell->face(face)->at_boundary() && 
	    (fun = neumann_boundaries.find(cell->face(face)->boundary_indicator())) != neumann_boundaries.end())
	  {

	    cell_rhs = 0;
	    fe_face_values.reinit(cell,face);
	    fun->second->value_list(fe_face_values.get_quadrature_points(),neumann_values);

	    for (size_t q=0; q<n_face_q_points; ++q)
	      {
		// For TODO(1), here neumann_values[q] should be set to dot(Dpi(x_q),fe_face_values.normal_vector(q)) 

		for (unsigned int i=0; i<n_cell_dof; ++i)
		  cell_rhs(i) += (neumann_values[q] *
				  fe_face_values.shape_value(i,q) *
				  fe_face_values.JxW(q));
	      }
	  }
	
      cell->get_dof_indices (global_dof_indices);
      for (size_t i=0; i<n_cell_dof; ++i)
	neumann_rhs(global_dof_indices[i]) += cell_rhs(i);
    }
  }

  fespace_member(void) increase_quadrature_order(const int n)
  {
    quadrature_formula = QGauss<dim>(quadrature_order+n);
    quadrature_order+=n;
    n_q_pts = quadrature_formula.size();

    get_positions();
  }


  fespace_member(void) update()
  {
    //    DoFRenumbering::Cuthill_McKee(dof_handler);

    dof_handler.distribute_dofs(fe);  
    
    n_cells        = triangulation.n_active_cells();
    n_q_pts        = quadrature_formula.size();
    n_cell_dof     = fe_values.dofs_per_cell; 
    n_dofs         = dof_handler.n_dofs();	   

    sparsity_pattern.reinit (dof_handler.n_dofs(),
			     dof_handler.n_dofs(),
			     dof_handler.max_couplings_between_dofs());

    DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);

    update_hanging_nodes();
    sparsity_pattern.compress();
    update_mass_matrix();
#if 0    
    for(map<uint_t,double>::const_iterator i = material.begin();i!=material.end();i++){
      const size_t material_id       = i->first;
      const double material_constant = i->second;
      if(material_constant <= 0) set_fixed_region(material_id,material_constant);
    }
#endif

//     cerr << "fixed_dof = \n";
//     for(map<size_t,double>::const_iterator i=fixed_dof.begin();i!=fixed_dof.end();i++){
//       cerr << "\t" << i->first  << " -> " << i->second << endl;
//     }
    update_laplace_matrix();

    update_boundary_conditions();
    get_positions();
  }
  // density(x) += f(x)
  fespace_member(void)
  LoadFunctionToMesh(const ScalarFunction& f, FEFunction& density) const
  {
    coordinate zero;
    ScalarFunctionWrap F(f); 
    FEFunction fe_f;
    QGauss<dim> quadrature(quadrature_order+2);
    density.resize(dof_handler.n_dofs());
    fe_f.resize(dof_handler.n_dofs());
    VectorTools::project(dof_handler,hanging_node_constraints, quadrature, F, fe_f.coefficients);

    density += fe_f;
  }

  // density(x) += weight*f(x-center)
  fespace_member(void)
  LoadFunctionToMesh(double weight, const ScalarFunction& f, 
				     const coordinate& center, FEFunction& density) const
  {
    ScalarFunctionWrap F(f,&center); 
    FEFunction fe_f;
    QGauss<dim> quadrature(quadrature_order+2);
    density.resize(dof_handler.n_dofs());
    fe_f.resize(dof_handler.n_dofs());
    VectorTools::project(dof_handler,hanging_node_constraints, quadrature, F, fe_f.coefficients);

    fe_f *= weight;
    density += fe_f;

    //hanging_node_constraints.condense(density.coefficients);
  }

  // f(x) = weight*exp(-exponent*|x-center|^2)
  // density += f(x)>=eps? f(x) : 0
  fespace_member(void)
  LoadGaussianToMesh(double weight, double exponent, 
		     const coordinate& center, 
		     double eps, FEFunction& density) const
  {
    GaussianFunction<dim> f2(exponent);
    ClampedGaussianFunction<dim> f(exponent, eps/weight);

    LoadFunctionToMesh(weight, f, center, density);
  }

  // density(x) += weight*f(|x-center|)*(|x-center|<range? 1:0)
  fespace_member(void)
  LoadScalar1DFunctionToMesh(FEFunction& density,  const Scalar1DFunctionClass& f, 
			     double range,  const coordinate& center,double weight) const
  {
    // Add support range
    RadialFunction<dim> f3D(f);
    LoadFunctionToMesh(weight, f3D, center, density);  
  }


  fespace_member(double) Integrate(const FEFunction& f) const 
  {
    double sum = 0;
    typename FEFunction::const_iterator c = f.begin();
    typename FEFunction::const_iterator weight = fe_weights.begin();

    for(; c != f.end(); c++,weight++) sum += (*c)*(*weight);

    return sum;
  }

  fespace_member(double) Integrate(const FEFunction& f, const FEFunction& g) const 
  {
    return overlap_matrix.matrix_scalar_product(f.coefficients,g.coefficients);
  }

  fespace_member(double)
  Integrate(const FEFunction& f, const FEFunction& V, const FEFunction& g) const
  {
    return 0;			// TODO: FEFunction -> PointFunction; 
  }

  fespace_member(void)
  ConstructPointFunction(const FEFunction& f, PointFunction& fp) const
  {
    FEValues<dim> fe_values(fe, quadrature_formula, update_values | update_quadrature_points); // evt. update_values
    vector<double> f_values(n_q_pts);

    fp.resize(n_cells*n_q_pts);
    typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();

    for (size_t pt=0; cell!=endc; ++cell){
      fe_values.reinit (cell);
      fe_values.get_function_values(f.coefficients, f_values);
      
      for(size_t q=0;q<n_q_pts;++q,++pt) fp[pt] = f_values[q];
    }
  }

  fespace_member(void)
  ConstructFEFunction(const PointFunction& fp, FEFunction& f) const 
  {
    f.resize(dof_handler.n_dofs());
    // Move computation of X to update()?
    FullMatrix<double> X(n_cell_dof,n_q_pts); // TODO: n_q_pts -> n_cell_q; Homogenisering af navne
    FETools::compute_projection_from_quadrature_points_matrix(fe,quadrature_formula,quadrature_formula,X);
    {
      typename DoFHandler<dim>::active_cell_iterator
	cell = dof_handler.begin_active(),
	endc = dof_handler.end();

      Vector<double> U(n_q_pts), V(n_cell_dof);
      std::vector<unsigned int> global_dof_indices (n_cell_dof);

      for (size_t pt=0; cell!=endc; ++cell, pt += n_q_pts){
	cell->get_dof_indices (global_dof_indices);

	//	copy(&fp[pt],&fp[pt+n_q_pts],U.begin());
	for(size_t i=0;i<n_q_pts;i++) U[i] = fp[pt+i];

	X.vmult(V,U);
	for(size_t i=0;i<n_cell_dof;i++)
	  f[global_dof_indices[i]] = V[i]; // XXX: Hmm... skal de virkelig laegges sammen? Eller maaske averages?
      }
    }
  }

  fespace_member(void) get_positions()
  {
    FEValues<dim> fe_values(fe, quadrature_formula, update_quadrature_points|update_values|update_JxW_values);

    point_positions.resize(n_cells*n_q_pts);
    point_weights.resize(n_cells*n_q_pts);
    typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();

    map<coordinate,bool> point_seen;

    for (size_t pt=0; cell!=endc; ++cell){
      fe_values.reinit (cell);
      const vector< Point< dim > > &q_points = fe_values.get_quadrature_points();

      for(size_t q=0;q<n_q_pts;++q,++pt){
	double w = 0;
      	point_positions[pt] = WrapPoint<SelfType>(q_points[q]);
	point_seen[point_positions[pt]] = true;

	for(size_t i=0;i<n_cell_dof;i++) w += fe_values.shape_value(i,q)*fe_values.JxW(q);
	point_weights[pt] = w;
      }
    }
    cerr << point_positions.size() << " quadrature points; " << point_seen.size() << " unique points." << endl;
  }

  fespace_member(double)
  Integrate(const FEFunction& f, const FEOperator& V, const FEFunction& g) const
  {
//     SparsityPattern::const_iterator ij = sparsity_pattern.begin(), end = sparsity_pattern.end();
//     double sum = 0;
//     for(;ij != end;ij++){
//       size_t i = ij->row(), j = ij->column();
//       sum += f[i]*g[j]*V(i,j);
//     }
//     return sum;
    return V.matrix_scalar_product(f.coefficients,g.coefficients);
  }


  fespace_member(double)
  LaplaceElement(const size_t i, const size_t j) const 
  {
    return laplace_matrix(i,j);
  }

  fespace_member(void)
  SolvePoisson(const FEFunction& density, FEFunction& result) const
  {
    SparseMatrix lhsmatrix;	// TODO: Find ud af, hvad der er galt, i stedet for det her heis.
    lhsmatrix.reinit(laplace_matrix.get_sparsity_pattern());
    lhsmatrix.copy_from(laplace_matrix);

    result.resize(dof_handler.n_dofs());

    dofVector rhs(n_dofs);
    overlap_matrix.vmult(rhs,density.coefficients);

    if(!neumann_boundaries.empty()){
      cerr << "Adding nonhomogeneous von Neumann boundary conditions."<<endl;
      rhs += neumann_rhs; 	// Add nonhomogeneous von Neumann boundary conditions
    }

    hanging_node_constraints.condense(rhs);

    MatrixTools::apply_boundary_values (boundary_values,
					lhsmatrix,
					result.coefficients,
					rhs);

    // SOLVE
    SolverControl solver_control (1000/*max iterations - should depend on problem size?*/, 1e-12/*tolerance*/);
    SolverCG<>    cg (solver_control);
    
    PreconditionSSOR<> prec;
    prec.initialize(lhsmatrix); // magic parameter -- investigate
    cg.solve (lhsmatrix, result.coefficients, rhs, prec);


    //    printf("\\int u = %g\n",Integrate(result)); 
    hanging_node_constraints.distribute (result.coefficients); // distribute solution
    //    printf("\\int u = %g\n",Integrate(result)); 

  }

  fespace_member(void)
  SolvePoisson(const PointFunction& density, FEFunction& result) const
  {
    SparseMatrix lhsmatrix;
    lhsmatrix.reinit(laplace_matrix.get_sparsity_pattern());
    lhsmatrix.copy_from(laplace_matrix);

    FEValues<dim> fe_values(fe, quadrature_formula, update_values | update_JxW_values);

    result.resize(dof_handler.n_dofs());
    std::vector<unsigned int> global_dof_indices (n_cell_dof);
    dofVector  rhs(n_dofs);
    cellVector cell_rhs (n_cell_dof);

    // Calculate rhs[i] = \int_{supp(\phi_i)} \phi_i(x) density(x) for each 1<=i<=n_dofs.
    typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();

    size_t pt = 0;
    for (size_t cellnum=0; cell!=endc; ++cell,++cellnum){
      fe_values.reinit (cell);
      // For each local i with supp(\phi_{GLOBAL[i]}) supported by cell:
      //   1. Calculate restriction of rhs[GLOBAL[i]] to cell
      //   2. Find global index.
      //   3. Add local contribution to rhs[i].

      cell_rhs = 0;
      for(size_t q=0;q<n_q_pts;++q,++pt){
	double dV = fe_values.JxW(q);
	for(size_t i=0;i<n_cell_dof;i++)
	  cell_rhs[i] += density[pt]*fe_values.shape_value(i,q)*dV;
      }
      cell->get_dof_indices (global_dof_indices);
      for(size_t i=0;i<n_cell_dof;i++) rhs[global_dof_indices[i]] += cell_rhs[i];
    }
    // Finished constructing right hand side.
    if(!neumann_boundaries.empty()) rhs += neumann_rhs; 	// Add nonhomogeneous von Neumann boundary conditions

    hanging_node_constraints.condense(rhs);

//     cerr << "Fixed dofs:\n";
//     for(map<uint_t,double>::const_iterator m=boundary_values.begin();
// 	m!=boundary_values.end();m++){
//       cerr << m->first << " -> " << m->second << endl;
//     }
    MatrixTools::apply_boundary_values (boundary_values,
					lhsmatrix,
					result.coefficients,
					rhs);

    // SOLVE
    SolverControl solver_control (2000/*max iterations - should depend on problem size?*/, 1e-10/*tolerance*/);
    SolverCG<>    solver (solver_control);
      
    PreconditionSSOR<> prec;
    prec.initialize(lhsmatrix);
    solver.solve (lhsmatrix, result.coefficients, rhs, prec);

    hanging_node_constraints.distribute (result.coefficients); // distribute solution
  }

  fespace_member(void)
  SolvePoisson(const PointFunction& density, PointFunction& result) const {
    FEFunction result_fe;
    SolvePoisson(density,result_fe);
    ConstructPointFunction(result_fe,result);
    
  }

  // Grid stuff
  fespace_member(void) 
  absolute_error_estimate(const FEFunction& fe_function, const ScalarFunction& function, 
			  cellVector& error_vector/*[n_active_cells()]*/) const
  {
    QGauss<dim> fine_quadrature(quadrature_formula.size()+1);
    const size_t N_q = fine_quadrature.size();

    error_vector.reinit(triangulation.n_active_cells());

    vector<double>       approx_values(N_q), real_values(N_q);
    vector< Point<dim> > qpoints(N_q);

    ScalarFunctionWrap function_dii(function);

    typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
  
    FEValues<dim> fe_values(fe, fine_quadrature, update_quadrature_points | update_values | update_JxW_values);
    
    for (size_t cellnum=0; cell!=endc; ++cell,++cellnum){
      fe_values.reinit (cell);
      qpoints = fe_values.get_quadrature_points();
      fe_values.get_function_values(fe_function.coefficients, approx_values);
      function_dii.value_list(qpoints,real_values);

      double error = 0;
      for(size_t q=0;q<N_q;q++)
	error += ::fabs((approx_values[q]-real_values[q])*fe_values.JxW(q));

      error_vector(cellnum) = error;
    }
  }
  // TODO: Refactor refinement routines -- most code is common to all schemes.
  fespace_member(void) refine_around_regions(const VolumeFunction<dim>& V, const double max_diameter, 
					     const bool update_at_end){
    typedef typename DoFHandler<dim>::active_cell_iterator cell_iterator;
    list<cell_iterator> parents, candidates;
    cell_iterator cell;

    // Initially, all cells are candidates for refinement
    for(cell = dof_handler.begin_active(); cell != dof_handler.end(); cell++)
      candidates.push_back(cell);
    
    bool done = false;
    while(!done){
      done = true;
      parents.clear();
      size_t n_refined = 0;

      for(typename list<cell_iterator>::const_iterator c = candidates.begin(); c != candidates.end(); c++){
	const cell_iterator& cell(*c);

	set<size_t> intersecting_volumes;
	for (size_t v=0;v < GeometryInfo<dim>::vertices_per_cell; v++)
	  intersecting_volumes.insert(V.volume_number(WrapPoint<SelfType>(cell->vertex(v))));

	// Cell spans at least two regions and is too large
	if(intersecting_volumes.size() >= 2 && cell->diameter() > max_diameter){
	  cerr << "Volumes intersected by cell: "; 
	  for(set<size_t>::const_iterator i=intersecting_volumes.begin();i!=intersecting_volumes.end();i++)
	    cerr << *i << " ";
	  cerr << endl;
	  
	  cell->set_refine_flag();
	  parents.push_back(cell);
	  done = false;
	  n_refined++;
	}
      }
      cerr << "Regions: Refining "<<n_refined<<" cells." << endl;
      triangulation.execute_coarsening_and_refinement();
      
      candidates.clear();
      for (typename list<cell_iterator>::iterator p(parents.begin());p!=parents.end();p++){
	for(size_t i=0;i<(1<<dim);i++)
	  candidates.push_back((*p)->child(i));
      }
    }
    if(update_at_end) update();
  }

  double approximate_cell_distance(const CellAccessor<3>& c, const Point<3>& X);


  // This refines to a mesh with the following two conditions:
  //   1. Each cell contains at most one point in xs (e.g. the nuclei in a molecule).
  //   2. Each cell containing one point in xs has volume no more than maxV.
  fespace_member(void) refine_around_points(const vector<coordinate> xs, const double max_diameter, 
					    const double near_diameter,
					    const bool update_at_end)
  {
    typedef typename DoFHandler<dim>::active_cell_iterator cell_iterator;
    list<cell_iterator> parents, candidates;
    cell_iterator cell;

    cerr << "Refining around points. Maximum diameter of cells in vicinity " 
	 << near_diameter << " of centers is " << max_diameter << endl;

    // Initially, all cells are candidates for refinement
    for(cell = dof_handler.begin_active(); cell != dof_handler.end(); cell++)
      candidates.push_back(cell);
    
    bool done = false;
    while(!done){
      done = true;
      parents.clear();
      size_t n_refined = 0;

      for(typename list<cell_iterator>::const_iterator c = candidates.begin(); c != candidates.end(); c++){
	size_t num_points_in_cell = 0, num_points_near_cell = 0;
	const cell_iterator& cell(*c);
	
	for(typename vector<coordinate>::const_iterator x=xs.begin();x!=xs.end();x++){
	  const PointWrap<SelfType> X_a(*x);
	  
	  if(cell->point_inside(X_a)) num_points_in_cell++;
	  if(approximate_cell_distance(*cell,X_a) <= near_diameter){
	    //  cerr << "Cell is close!" << endl;
	    num_points_near_cell++;
	  }
	}

	if(num_points_in_cell>=2 || (num_points_near_cell>=1 && cell->diameter()>max_diameter)){
	  //	  cerr << "Refining cell." << endl;
	  cell->set_refine_flag();
	  parents.push_back(cell);
	  done = false;
	  n_refined++;
	}
      }
      cerr << "Refining "<<n_refined<<" cells." << endl;
      triangulation.execute_coarsening_and_refinement();
      
      candidates.clear();
      for (typename list<cell_iterator>::iterator p(parents.begin());p!=parents.end();p++){
	for(size_t i=0;i<(1<<dim);i++)
	  candidates.push_back((*p)->child(i));
      }
    }
    if(update_at_end) update();
  }

  // This is terribly close in functionality to 
  //  1. absolute_error_estimate(0,density,error_vector)
  //  2. for each i: if(error_vector[i]>2^dim dE) mark cell[i]
  //  3. if no cells are marked, stop.
  //  4. refine
  //  5. repeat from 1.
  // TODO: Coarsen after refinement has ended.
  fespace_member(void) refine_to_density(const ScalarFunction& density, const double dE, bool update_at_end)
  {
    // Refines mesh with the goal of making each cell contain close to $dE$ of $density$.
    const size_t N_q = n_q_pts;

    fprintf(stderr,"Refine to density %g.\n",dE);
    vector<double>       approx_values(N_q), real_values(N_q);
    vector< Point<dim> > qpoints(N_q);

    ScalarFunctionWrap density_dii(density);
    FEValues<dim> fe_values(fe, quadrature_formula, update_quadrature_points | update_JxW_values);

    // For each active cell:
    //  + Integrate density over cell volume (values at quadrature points times weights).
    //  + If integral is $2^dim dE$ or greater, mark for subdivision.
    bool done_refining = false;

    typename DoFHandler<dim>::active_cell_iterator
	cell = dof_handler.begin_active(),
	endc = dof_handler.end();      

    typedef typename DoFHandler<dim>::cell_iterator cell_iterator;
    list<cell_iterator> parents, candidates;
    for(;cell!=endc;cell++) candidates.push_back(cell);

    while(!done_refining){
      done_refining = true;
      size_t n_refined = 0;
      cerr << "Integrating densities and tagging cells for refinement..." << endl;
      cerr << candidates.size() << " candidates for refinement." << endl;
      parents.clear();

      double max_density = 0;

      for (typename list<cell_iterator>::iterator c(candidates.begin());c!=candidates.end();c++){
	DoFHandler<dim>::cell_iterator& cell(*c);
	fe_values.reinit (cell);
	qpoints = fe_values.get_quadrature_points();
	density_dii.value_list(qpoints,real_values);

	double cell_density = 0;
	for(size_t q=0;q<N_q;q++)
	  cell_density += ::fabs(real_values[q])*fe_values.JxW(q);

	if(cell_density >= pow(1.9,dim)*dE){
	  cell->set_refine_flag();
	  parents.push_back(cell);
	  done_refining = false;
	  n_refined++;

	  if(cell_density > max_density) max_density = cell_density;
	} 
      }
      cerr << "Refining "<<n_refined<<" cells." << endl;
      triangulation.execute_coarsening_and_refinement();
      
      candidates.clear();
      fprintf(stderr,"(max_density,8*dE) = (%g,%g)\n",max_density,8.*dE);
      if(max_density < pow(2.0,dim)*dE) done_refining = true;
      else 
	for (typename list<cell_iterator>::iterator p(parents.begin());p!=parents.end();p++){
	  for(size_t i=0;i<(1<<dim);i++)
	    candidates.push_back((*p)->child(i));
	}


    }
    if(update_at_end) update();
  }

  fespace_member(void) refine_to_density(PointFunctional& density, const double dE, bool update_at_end){
    // Refines mesh with the goal of making each cell contain close to $dE$ of $density$.
    const size_t N_q = n_q_pts;

    fprintf(stderr,"Refine to density %g.\n",dE);
    vector<double>       approx_values(N_q), real_values(N_q);
    vector< Point<dim> > qpoints(N_q);

    FEValues<dim> fe_points(fe, quadrature_formula, update_quadrature_points);
    FEValues<dim> fe_values(fe, quadrature_formula, update_JxW_values);

    // For each active cell:
    //  + Integrate density over cell volume (values at quadrature points times weights).
    //  + If integral is $2^dim dE$ or greater, mark for subdivision.
    bool done_refining = false;

    typename DoFHandler<dim>::active_cell_iterator
	cell = dof_handler.begin_active(),
	endc = dof_handler.end();      

    typedef typename DoFHandler<dim>::cell_iterator cell_iterator;
    list<cell_iterator> parents, candidates;
    vector<coordinate> positions;

    for(cell = dof_handler.begin_active();cell!=endc;cell++) candidates.push_back(cell);
    // Now recursively coarsen dense cells
    while(!done_refining){
      done_refining = true;
      size_t n_refined = 0;
      cerr << candidates.size() << " candidates for refinement." << endl;
      parents.clear();

      double max_density = 0;

      // Factor the next 8 lines out
      positions.resize(candidates.size()*N_q);
      typename list<cell_iterator>::iterator c;
      size_t i;

      for (c = candidates.begin(), i = 0;c!=candidates.end();c++,i++){
	fe_points.reinit(*c);
	qpoints = fe_points.get_quadrature_points();
	for(size_t q=0;q<N_q;q++)
	  for(size_t j=0;j<dim;j++) positions[i*N_q+q][j] = qpoints[q][j];
      }
      const PointFunction& density_pt(density.fill(positions));

      for (c = candidates.begin(), i=0;c!=candidates.end();c++){
	DoFHandler<dim>::cell_iterator& cell(*c);
	fe_values.reinit (cell);

	double cell_density = 0;
	for(size_t q=0;q<N_q;q++,i++)
	  cell_density += ::fabs(density_pt[i])*fe_values.JxW(q);

	if(cell_density >= pow(1.9,dim)*dE){
	  cell->set_refine_flag();
	  parents.push_back(cell);
	  done_refining = false;
	  n_refined++;

	  if(cell_density > max_density) max_density = cell_density;
	} 
      }
      cerr << "Refining "<<n_refined<<" cells." << endl;
      triangulation.execute_coarsening_and_refinement();
      
      candidates.clear();
      fprintf(stderr,"(max_density,8*dE) = (%g,%g)\n",max_density,8.*dE);
      if(max_density < pow(2.0,dim)*dE) done_refining = true;
      else 
	for (typename list<cell_iterator>::iterator p(parents.begin());p!=parents.end();p++){
	  for(size_t i=0;i<(1<<dim);i++)
	    candidates.push_back((*p)->child(i));
	}


    }
    if(update_at_end) update();    
  }

  
  fespace_member(void) refine_grid(const cellVector& estimated_error_per_cell, bool update_at_end)
  {
    GridRefinement::refine_and_coarsen_fixed_number (triangulation,estimated_error_per_cell, 0.30, 0.10);
    triangulation.execute_coarsening_and_refinement ();
    if(update_at_end) update();    
  }

  fespace_member(void) refine_grid(const Vector<float>& estimated_error_per_cell, bool update_at_end)
  {
    GridRefinement::refine_and_coarsen_fixed_number (triangulation,estimated_error_per_cell, 0.30, 0.10);
    triangulation.execute_coarsening_and_refinement ();
    if(update_at_end) update();    
  }

  fespace_member(void) refine_grid(const size_t n, bool update_at_end)
  {
    triangulation.refine_global(n);
    if(update_at_end) update();
  }

  fespace_member(double) Value(const FEFunction& f, const coordinate& x) const 
  {
    return VectorTools::point_value(dof_handler, f.coefficients, PointWrap<FESpace>(x));
  }


}
#undef fespace_member


