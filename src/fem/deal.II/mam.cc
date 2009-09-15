    {

      std::map<unsigned int, double> boundary_values;
      // points to constrain a dof

      Point<dim> fix_point_A(-3.0, -7.7, 3.0); // point to fix

      DoFHandler<dim>::active_cell_iterator
        cell = dof_handler.begin_active(), endc = dof_handler.end();


      std::vector<bool> vertex_touched(triangulation.n_vertices(), false);

      for (unsigned int cellCount = 0; cell != endc; ++cell, ++cellCount) { 
	for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; ++v) {

	  if (vertex_touched[cell->vertex_index(v)] == false) {
	    vertex_touched[cell->vertex_index(v)] = true;

	    Point<dim> vertex_position = cell->vertex(v); 
	    if (vertex_position(0) == fix_point_A(0) && vertex_position(1) == fix_point_A(1)) { 
	      boundary_values[cell->vertex_dof_index(v, 0)] = 0.0;
	    }
	  }
	}

      }


      MatrixTools::apply_boundary_values(boundary_values, K, delta_d, f_residual);

    }
