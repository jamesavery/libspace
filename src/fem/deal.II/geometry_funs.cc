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

#  include <iostream>
#  include <fstream>
#  include <grid/grid_out.h>
#  include <numerics/data_out.h>
#  include <lac/sparsity_pattern.h>
#  include <grid/grid_refinement.h>
#  include <grid/grid_tools.h>
#  include <grid/grid_in.h>
#  include <dofs/dof_renumbering.h>

#  include <fe/fe_tools.h>
#  include <lac/solver_selector.h>

namespace dealii { 
  using namespace std;

class Line {
public:
  Point<3> a,b;

  Line(const Point<3>& a, const Point<3>& b) : a(a), b(b) { }

  Point<3> operator()(const double t) const {
    return a + (b-a)*t;
  }
    
};


class Plane { 
public:
  Point<3> n;	        // Normal
  double d;		// Distance from origo

  Plane(const Point<3>& p0, const Point<3>& p1, const Point<3>& p2){
    const Point<3> U(p1-p0), V(p2-p0);
      
    n[0] = U[1]*V[2]-U[2]*V[1];
    n[1] = U[2]*V[0]-U[0]*V[2];
    n[2] = U[0]*V[1]-U[1]*V[0];

    d = n*p0;
  }

  double intersection(const Line& L) const {
    const Point<3> &x(L.a), &y(L.b);
    const double t = 		// XXX: Check om L er parallel med n
      (d - n[0]*x[0]-n[1]*x[1]-n[2]*x[2])/
      (n[0]*(y[0]-x[0]) + n[1]*(y[1]-x[1]) + n[2]*(y[2]-x[2]));
      
    return t;
  }
};



double approximate_cell_distance(const CellAccessor<3>& c, const Point<3>& X)
{
  double min_distance_sqr = INFINITY;

  if(c.point_inside(X)){ return 0; }

  const double delta = 1e-4;
  const Line l(X,c.center());
  // Distances to vertices
  for(size_t v=0;v<GeometryInfo<3>::vertices_per_cell;v++){
    const double d = X.distance(c.vertex(v));
    if(d*d < min_distance_sqr) min_distance_sqr = d*d;
  }
  // Distances to lines. Scheme:
  
  // Distances to faces. Scheme:
  //   1. Define line l(t) = X+t(c.center()-X)
  //   For each face:
  //       2.1. find intersection with l and plane defined by face.
  //       2.2. if t\in[0;1] and c.point_inside(l(t+delta)) then:
  //            2.2.1 p_sqr = |l(t)|^2 
  //            2.2.2 if p_sqr<min_distance_sqr then set min_distance_sqr = p_sqr
  for(size_t f=0;f<GeometryInfo<3>::faces_per_cell;f++){
    const TriaAccessor<2,3,3>& face(*c.face(f));
    const Plane faceplane(face.vertex(0),face.vertex(1),face.vertex(2));
    const double t = faceplane.intersection(l);
     
    if(t>=0 && t<=1 && (c.point_inside(l(t+delta)) || c.point_inside(l(t-delta)))){
      const Point<3> p(l(t));
      const double   p_sqr = p*p;
      if(p_sqr<min_distance_sqr) min_distance_sqr = p_sqr;
    }
  }

  return sqrt(min_distance_sqr);
}
};
