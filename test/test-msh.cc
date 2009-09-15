#include <space/fem/deal.II/mesh.h>
#include <space/function.h>
#include <ios>
using namespace std;

typedef dealii::FESpace<3> Space3D;
typedef Space3D::coordinate coordinate;
typedef dealii::Point<2> point_dii;

int main()
{
  Space3D space("test.msh",1,2);
  
  cerr << "   Number of active cells: "
       << space.triangulation.n_active_cells()
       << std::endl
       << "   Total number of cells: "
       << space.triangulation.n_cells()
       << endl;

  space.write_mesh("test2.msh");
  
  return 0;
}
