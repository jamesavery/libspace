#include <space/storage/pfinterval.h>

#define member(T) T PFInterval::

member(void) update(const double epsilon){
  eps = epsilon;
  bool zero_streak = true;
  int a=-1;
  intervals.clear();
  for(size_t i=0;i<coefficients.size();i++){
    if(zero_streak){
      if(fabs(coefficients[i]) > eps){
	zero_streak = false;
	a = i;
      }
    } else if(zero_streak_start(coefficients,i,eps)){
      zero_streak = true;
      intervals.push_back(interval_t(a,i));
      a = -1;
    }
  }
  // Commit last half finished interval, if any
  if(zero_streak==false) 
    intervals.push_back(interval_t(a,coefficients.size()));
}

member(void) set_coefficients(const ContainerClass& cs, double epsilon){
  eps = epsilon;
  coefficients = cs;
  update(epsilon);
}

member(double) Integrate(const ContainerClass& weights) const {
  double sum = 0;
  for(intervallist_t::const_iterator i=intervals.begin();i!=intervals.end();i++)
    for(size_t j=i->a;j<i->b;j++) sum += coefficients[j] * weights[j];

  return sum;
}

member(PFInterval&) operator +=(const PFInterval& y){
  intervals += y.intervals;
  // y.interval <= this->interval
  for(intervallist_t::const_iterator i=y.intervals.begin();i!=y.intervals.end();i++)
    for(size_t j=i->a;j<i->b;j++) coefficients[j] += y.coefficients[j];

  return *this;
}

member(PFInterval&) operator -=(const PFInterval& y){
  intervals += y.intervals;
  // y.interval <= this->interval
  for(intervallist_t::const_iterator i=y.intervals.begin();i!=y.intervals.end();i++)
    for(size_t j=i->a;j<i->b;j++) coefficients[j] -= y.coefficients[j];

  return *this;
}

member(PFInterval&) operator *=(const PFInterval& y){
  intervallist_t iprod(intervals * y.intervals);
  //  cerr << "Complement("<<intervals <<"; "<<iprod<<") = " << (intervals-iprod) <<endl;
  intervallist_t zero_out(intervals-iprod);

  intervals = iprod;
  // Set A\B to zero
  //TODO:  zero_out -= y.intervals;
  for(intervallist_t::const_iterator i=zero_out.begin();i!=zero_out.end();i++)
    for(size_t j=i->a;j<i->b;j++) coefficients[j] = 0; // Evt. memset

  // Support is B\cup A
  for(intervallist_t::const_iterator i=intervals.begin();i!=intervals.end();i++)
    for(size_t j=i->a;j<i->b;j++) coefficients[j] *= y.coefficients[j];

  return *this;
}

member(PFInterval&) operator *=(const double& y){
  for(intervallist_t::const_iterator i=intervals.begin();i!=intervals.end();i++)
    for(size_t j=i->a;j<i->b;j++) coefficients[j] *= y;

  if(fabs(y)<=eps) intervals.clear();
  return *this;
}


member(size_t) number_nonzero() const {
  size_t nonzero = 0;
  for(intervallist_t::const_iterator i=intervals.begin();i!=intervals.end();i++)
    nonzero += (i->b - i->a);

  return nonzero;
}
member(size_t) number_intervals() const { return intervals.size(); }

member(PFInterval&) operator =(const double& rhs){ 
  for(ContainerClass::iterator c = coefficients.begin();c != coefficients.end(); c++)
    *c = rhs; 

  update(eps);
  return *this; 
}
