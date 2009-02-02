/* \file regulargrid.h
   \brief Regular grid class.
 */
#ifndef LIBDISC_REGULARGRID_H
# define LIBDISC_REGULARGRID_H

#include <sys/types.h>
#include <libdiscretization/discretization.h>

#include <math.h>

/* TODO: Make distributed. */
template <int dim, class PointFunction_ = PointFunction_Simple<> > 
  class RegularGrid : public Discretization<dim,PointFunction_>  {	/* Specialized to double -- perhaps templatize field instead. */
 public:
 typedef typename Discretization<dim,PointFunction_>::coordinate       coordinate;
 typedef typename Discretization<dim,PointFunction_>::ScalarFunction   ScalarFunction;
 typedef typename Discretization<dim,PointFunction_>::DiscreteFunction PointFunction; 

 RegularGrid(const size_t npts_[dim], const coordinate& leftcorner, const coordinate& dimensions) :
  leftcorner(leftcorner), dimensions(dimensions)
  {
    total_npts = 1;
    dV = 1;
    for(size_t i=0;i<dim;i++){
      npts[i] = npts_[i];
      total_npts *= npts[i];
      dx[i] = dimensions.x[i]/((double)npts[i]);
      dV   *= dx[i];
    }

  }
  
  /* Adding scalar functions to a mesh density */
  void LoadFunctionToMesh(double weight, const ScalarFunction& f, 
			  const coordinate& position, PointFunction& density) const;
  void LoadGaussianToMesh(double weight, double exponent, const coordinate& position, 
			  double eps, PointFunction& density) const;
  void LoadScalar1DFunctionToMesh(PointFunction& density,  Scalar1DFunctionClass& f, 
				  double range,  coordinate& cartesianPosition,double weight=1.0) const;

  /* Integration and inner products */
  double Integrate(const PointFunction& fvalues) const;
  double Integrate(const PointFunction& fvalues, const PointFunction& gvalues) const;
  double Integrate(const PointFunction& fvalues, const PointFunction& Vvalues, const PointFunction& gvalues) const;

  double LaplaceElement(const size_t i, const size_t j) const;
  /* Solving PDE's */
  void SolvePoisson(const PointFunction& density, PointFunction& vHartree);

  /* Functionality specific to regular grids */
  size_t npts[dim];			/**< Number of points on each axis. */
  coordinate leftcorner;		/**< Spatial position of [0,0,...,0]-node. */
  coordinate dimensions;		/**< Length of grid along each axis. */
  size_t total_npts;			/**< Total number of nodes: total_npts = npts[0]*...*npts[dim-1]. */
  double dx[dim];                       /**< Dimensions of a grid cell: dx[i] = dimensions.x[i]/npts[i]. */
  double dV;			        /**< Volume of a grid cell: dV = dx[0]*...*dx[dim-1]. */


  coordinate NodeToPosition(const size_t idx[dim]) const { /* NB: Trouble when moved to .cc */
    coordinate position(leftcorner);
    for(size_t i=0;i<dim;i++) position.x[i] += dimensions.x[i]*idx[i]/((double) npts[i]);

    return position;
  }
  coordinate IndexToPosition(const size_t global_idx) const { /* NB: Trouble when moved to .cc */
    size_t idx[dim];
    IndexToNode(global_idx,&idx[0]);
    return NodeToPosition(idx);
  }
  
  void   PositionToNode(const coordinate& x, size_t *idx) const;
  void   IndexToNode(const size_t global_idx, size_t *idx) const;
  size_t NodeToIndex(const size_t *idx) const;

  int    single_shared_cell(const size_t global_i, const size_t global_j) const;
};

#endif
