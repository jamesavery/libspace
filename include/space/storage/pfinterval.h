#ifndef LIBSPACE_PFINTERVAL_H
# define LIBSPACE_PFINTERVAL_H

#include <space/storage/storage.h>
#include <space/storage/interval.h>

/* Perhaps make template <class ContainerClass, typename Y=double, typename S=double */
class PFInterval : public DiscreteFunction<std::vector<double>,PFInterval> {
  typedef std::vector<double> ContainerClass;
  typedef DiscreteFunction<std::vector<double>,PFInterval> BaseType;

  static bool zero_streak_start(const ContainerClass& cs, const size_t i, const double eps=0){
    return((i+2<cs.size()) && ::fabs(cs[i]) <= eps && ::fabs(cs[i+1]) <= eps && ::fabs(cs[i+2]) <= eps);
  }
public:
  using BaseType::coefficients;

  typedef Interval<unsigned int>  interval_t;
  typedef std::list<interval_t> intervallist_t;

  Intervals<unsigned int> intervals;
  double eps;

  double  operator[](const size_t i) const { return coefficients[i]; }
  double& operator[](const size_t i) { return coefficients[i]; }

  PFInterval(const double eps=0): eps(eps){}
  PFInterval(const ContainerClass& cs, const double eps = 0) : BaseType(cs), eps(eps) 
  {
    update(eps);
  }

  void update(const double epsilon=0);

  void set_coefficients(const ContainerClass& cs, double epsilon = 0);
  void resize(const size_t n){ coefficients.resize(n); update(eps); }
  double Integrate(const ContainerClass& weights) const;

  PFInterval& operator +=(const PFInterval& y);
  PFInterval& operator -=(const PFInterval& y);
  PFInterval& operator *=(const PFInterval& y);
  PFInterval& operator *=(const double& y);

  PFInterval operator +(const PFInterval& y) const {
    PFInterval z(*this);
    return (z += y);
  }
  PFInterval operator -(const PFInterval& y) const {
    PFInterval z(*this);
    return (z -= y);
  }

  PFInterval operator *(const PFInterval& y) const {
    PFInterval z(*this);
    return (z *= y);
  }

  PFInterval operator *(const double& y) const {
    PFInterval z(*this);
    return (z *= y);
  }

  PFInterval& operator =(const double& rhs);

  friend std::ostream& operator<<(std::ostream& s, const PFInterval& X){ 
    s << "{\n\tIntervals: "<<X.intervals<<",\n"<<
         "\tValues: {";
    for(size_t i=0;i<X.coefficients.size();i++) s << X.coefficients[i]<<(i+1<X.coefficients.size()?", ":"}\n");
    s << "}";
    return s;
  }
  
  size_t number_nonzero() const;
  size_t number_intervals() const;
};

#endif
