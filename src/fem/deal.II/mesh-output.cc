#include <space/fem/deal.II/mesh.h>
#include <space/volumes.h>

#include <algorithm>
#include <numerics/vectors.h>
#include <numerics/fe_field_function.h>
#include <base/function.h>

#include <stdio.h>
#include <iostream>
#include <fstream>

#  include <grid/grid_out.h>
#  include <numerics/data_out.h>
#  include <lac/sparsity_pattern.h>
#  include <fe/fe_tools.h>

#define fespace_member(returntype) template <int dim> returntype FESpace<dim>::
int lookup_format(const string *supported_formats, const string& path);

namespace dealii{

  using namespace std;
  // Ouput stuff.
  fespace_member(void) write_dof_sparsity(const string& path) const 
  {
    const string supported_formats_str[] = {"txt","gpl",""};
    enum {TXT,GNUPLOT} supported_formats;

    ofstream file(path.c_str());
    
    switch(lookup_format(supported_formats_str,path)) {
    case TXT:
      sparsity_pattern.print(file);
      break;
    default:
      sparsity_pattern.print_gnuplot(file);
    };

    file.close();
  }

  fespace_member(void) write_mesh(const string& path) const 
  {
    const string supported_formats_str[] = {"dx","msh","ucd","eps","xfig","gpl",""};
    enum {DX,MSH,UCD,EPS,XFIG,GNUPLOT} supported_formats;

    ofstream file(path.c_str());
    GridOut out;

    switch(lookup_format(supported_formats_str,path)) {
    case DX:
      out.write_dx(triangulation,file);
      break;
    case MSH:
      // ... + indicators
      out.write_msh(triangulation,file);
      break;
    case UCD:
      out.write_ucd(triangulation,file);
      break;
    case EPS:
      out.write_eps(triangulation,file);
      break;
    case XFIG:
      out.write_xfig(triangulation,file);
      break;
    case GNUPLOT:
    default:
      out.write_gnuplot(triangulation,file);
    };

    file.close();
  }




  fespace_member(void) write_function(const string& path, const FEFunction& f) const
  {
    const string supported_formats_str[] = {"dx","eps","gmv","pov","plt","pltx","ucd","vtk","dealII","gpl",""};
    enum {DX,EPS,GMV,POVRAY,TECPLOT,TECPLOT_BINARY,UCD,VTK,DEALII,GNUPLOT} supported_formats;

    ofstream file(path.c_str());
    DataOut<dim> out;

    out.attach_dof_handler (dof_handler);
    out.add_data_vector (f.coefficients, "f");
    out.build_patches ();
   
    switch(lookup_format(supported_formats_str,path)) {
    case DX:               out.write_dx(file);       break;
    case EPS:              {
      DataOutBase::EpsFlags eps_flags;
      eps_flags.z_scaling = 4;

      out.set_flags (eps_flags);
      out.write_eps(file);  
    }
      break;
    case GMV:              out.write_gmv(file);      break;
    case POVRAY:           out.write_povray(file);   break;
    case TECPLOT:          out.write_tecplot(file); break;
    case TECPLOT_BINARY:   out.write_tecplot_binary(file); break;
    case UCD:              out.write_ucd(file);      break;
    case VTK:              out.write_vtk(file);      break;
    case DEALII:           out.write_deal_II_intermediate(file);  break;
    case GNUPLOT:
    default:
      out.write_gnuplot(file);
    };

    file.close();  
  }
}
