#include <stdio.h>
#include <sys/types.h>

#include <libdiscretization/discretization.h>
#include <libdiscretization/grid/regulargrid.h>

int main()
{
  double upperleft_[3]  = {-.5,-.5,-.5};
  double dimensions_[3] = {1,1,1};
  const size_t npts[3] = {20,20,20};

  typedef RegularGrid<3> grid3d;
  typedef RegularGrid<2> grid2d;
  typedef RegularGrid<1> grid1d;

  grid3d::coordinate upperleft(upperleft_);
  grid3d::coordinate dimensions(dimensions_);

  grid3d G(npts,upperleft,dimensions);
  grid3d::DiscreteFunction density;

  G.LoadGaussianToMesh(1.0,5.0,grid3d::coordinate(),0, density);

  printf("%g\n",G.Integrate(density));
  printf("%g\n",G.Integrate(density*density));

  // Test af NodeToIdx <-> IdxToNode
  size_t idx[3], global_idx = 313373;
  G.IndexToNode(global_idx,&idx[0]);
  printf("IndexToNode(%d) = [%d,%d,%d]\n",global_idx,idx[0],idx[1],idx[2]);
  printf("NodeToIndex([%d,%d,%d]) = %d\n",idx[0],idx[1],idx[2], G.NodeToIndex(&idx[0]));
  
  grid3d::coordinate x(G.IndexToPosition(global_idx));
  printf("IndexToPosition(%d) = [%g,%g,%g]\n",global_idx,x.x[0],x.x[1],x.x[2]);
  G.PositionToNode(x,&idx[0]);
  printf("PositionToNode([%g,%g,%g]) = [%d,%d,%d]\n",x.x[0],x.x[1],x.x[2],idx[0],idx[1],idx[2]);

  // Test af Laplace-element
  {
    printf("1D Laplacian\n");
    size_t npts[1] = {8};
    double left[1] = {0};
    double dims[1] = {8};
    grid1d::coordinate upperleft(left);
    grid1d::coordinate dimensions(dims);
    grid1d line(npts,upperleft,dimensions); // 8 point representation of [0;8]

    
    for(size_t i=0;i<npts[0];i++){
      printf("[ ");
      for(size_t j=0;j<npts[0];j++) printf("%2.2g ",line.LaplaceElement(i,j));
      printf("]\n");
    }
  }

  {
    printf("2D Laplacian\n");
    size_t npts[2] = {4,4};
    double left[2] = {0,0};
    double dims[2] = {4,4};
    grid2d::coordinate upperleft(left);
    grid2d::coordinate dimensions(dims);
    grid2d square(npts,upperleft,dimensions); // 16 point representation of [0;4]x[0;4]

    
    for(size_t i=0;i<npts[0]*npts[1];i++){
      printf("[ ");
      for(size_t j=0;j<npts[0]*npts[1];j++) printf("%2.2g ",square.LaplaceElement(i,j));
      printf("]\n");
    }
  }

  {
    printf("3D Laplacian\n");
    size_t npts[3] = {3,3,3};
    double left[3] = {0,0,0};
    double dims[3] = {3,3,3};
    grid3d::coordinate upperleft(left);
    grid3d::coordinate dimensions(dims);
    grid3d square(npts,upperleft,dimensions); // 27 point representation of [0;3]x[0;3]x[0;3]

    
    for(size_t i=0;i<npts[0]*npts[1]*npts[2];i++){
      printf("[ ");
      for(size_t j=0;j<npts[0]*npts[1]*npts[2];j++) printf("%2.2g ",square.LaplaceElement(i,j));
      printf("]\n");
    }
  }



  return 0;
}
