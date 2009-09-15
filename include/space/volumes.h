#ifndef LIBSPACE_VOLUMES_H
# define LIBSPACE_VOLUMES_H

#include <space/storage/smallvector.h>


template <int dim> class Volume {
 public:
  typedef SmallVector<dim> coordinate;

  coordinate bounding_box[2];
  
  bool point_inside_bounding_box(const coordinate& x) const {
    const coordinate& upperleft(bounding_box[0]);
    const coordinate& lengths(bounding_box[1]);
    
    for(size_t i=0;i<dim;i++)
      if(x[i] < upperleft[i] || x[i] > upperleft[i]+lengths[i]) return false;
    
    return true;
  }
  virtual bool point_inside_volume(const coordinate& x)  const = 0;

  /* Default implementation returns 0 if inside or distance to bounding box if outside */
  /* TODO: Check validity of this hack. */
  virtual double approximate_distance(const coordinate& x) const {
    double D = 0;
    if(point_inside_volume(x)) return 0.0;

    const coordinate &upperleft(bounding_box[0]), &lengths(bounding_box[1]);
    for(size_t i=0;i<dim;i++){
      const double d1 = fabs(upperleft[i]-x[i]);
      const double d2 = fabs(upperleft[i]+lengths[i]-x[i]);
      D += d1<d2? d1*d1:d2*d2;
    }
    return sqrt(D);
  }

  friend std::ostream& operator<<(std::ostream& S, const Volume& v){
    S << "volume with bounding box(corner:" << v.bounding_box[0] << "; lengths: " << v.bounding_box[1] << ")";
    return S;
  } 
};

template <int dim> class BoxVolume : public Volume<dim> {
public:
  typedef typename Volume<dim>::coordinate coordinate;

  BoxVolume(const coordinate& upperleft, const coordinate& lengths){
    Volume<dim>::bounding_box[0] = upperleft;
    Volume<dim>::bounding_box[1] = lengths;
  }
  bool point_inside_volume(const coordinate& x) const { 
    return Volume<dim>::point_inside_bounding_box(x); 
  }
};

template <int dim> class EllipsoidVolume :  public Volume<dim> {
public:
  typedef typename Volume<dim>::coordinate coordinate;

  coordinate center, R;

  // Volume defined by \f$\sum_{i=1}^{dim} (x_i-center_i)^2/R_i^2 \le 1\f$
  EllipsoidVolume(const coordinate& center, const coordinate& R) : center(center), R(R) {
    // Bounding box has faces at
    //    (x_i/R_i)^2 + 0 + 0 == 1  <=> |x_i| = |R_i|
    Volume<dim>::bounding_box[0] = center - R; // upper left corner
    Volume<dim>::bounding_box[1] = R*2.0;	  // length of BB sides
  }

  bool point_inside_volume(const coordinate& x) const {
    const coordinate r(x-center);
    double l = 0;
    for(size_t i=0;i<dim;i++) l += r[i]*r[i]/(R[i]*R[i]);
    return (l<=1.0);
  }

};

template <int dim> class ConstantVolume {
 public:
  double value;
  Volume<dim> *shape;

 ConstantVolume(double value, Volume<dim> *shape) : value(value), shape(shape) {}
  
  friend std::ostream& operator<<(std::ostream& S, const ConstantVolume& v){
    S <<  *v.shape << " has value " << v.value;
    return S; 
  }
};

  

#endif

