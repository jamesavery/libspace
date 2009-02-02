
  template <int dim, int spacedim> void myproject(const Mapping<dim, spacedim>       &mapping,
				    const DoFHandler<dim,spacedim>    &dof,
				    const ConstraintMatrix   &constraints,
				    const SparseMatrix<double>   &mass_matrix,
				    const Vector<double>         &mass_matrix_rhs,
				    const SparsityPattern        &sparsity_pattern,
				    const Quadrature<dim>    &quadrature,
				    const dealii::Function<spacedim>      &function,
				    Vector<double>                   &vec_result)
  {
    printf("myproject()\n");
    Assert (dof.get_fe().n_components() == function.n_components,
	    ExcDimensionMismatch(dof.get_fe().n_components(),
				 function.n_components));

    Assert (vec_result.size() == dof.n_dofs(),
	    ExcDimensionMismatch (vec_result.size(), dof.n_dofs()));
  
    Vector<double> vec(dof.n_dofs());
    Vector<double> tmp(mass_matrix_rhs);
    SparseMatrix<double> mm(mass_matrix);

    for (unsigned int i=0; i<vec.size(); ++i)
      vec_result(i) = vec(i);
    // Allow for a maximum of 5*n
    // steps to reduce the residual by
    // 10^-12. n steps may not be
    // sufficient, since roundoff
    // errors may accumulate for badly
    // conditioned matrices
    ReductionControl        control(5*tmp.size(), 0., 1e-12, false, false);
    GrowingVectorMemory<> memory;
    SolverCG<>              cg(control,memory);

    PreconditionSSOR<> prec;
    prec.initialize(mm, 1.2);

    // solve
    cg.solve (mm, vec, tmp, prec);
  
    // distribute solution
    constraints.distribute (vec);

    // copy vec into vec_result. we
    // can't use ve_result itself
    // above, since it may be of
    // another type than Vector<double>
    // and that wouldn't necessarily go
    // together with the matrix and
    // other functions
    for (unsigned int i=0; i<vec.size(); ++i)
      vec_result(i) = vec(i);

    printf("myproject() end\n");
  }
