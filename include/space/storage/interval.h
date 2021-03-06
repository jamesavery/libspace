#ifndef LIBSPACE_INTERVAL_H
# define LIBSPACE_INTERVAL_H

#include <list>
#include <iostream>

/// Models half-open interval [a;b[ 
/// When Q is int, it is equivalent to the closed interval [a;b-1]
template <typename Q> class Interval { // Models half-open interval [a;b[ -- for integers equivalent to closed [a;b-1]
 public:
  Q a,b;
  
  Interval(Q a=0,Q b=0) : a(a), b(b){}

  /* Intersection */
  Interval& operator *=(const Interval& B){
    if(B.a>a) a = B.a;
    if(B.b<b) b = B.b;
    
    return *this;
  }
  Interval operator *(const Interval& B) const {
    Interval<Q> C(a,b);
    return (C *= B);
  }

  friend std::ostream& operator<<(std::ostream& s, const Interval<Q>& I){ 
    s << "[" << I.a << "," << I.b << "[";

    return s;
  }
};


/// Models ordered sets of half-open intervals
/// Operator + is union and operator * is intersection.
template <typename Q> class Intervals : public std::list< Interval<Q> > {
  typedef enum {A_STRICTLY_LESS, B_STRICTLY_LESS, A_LESS_OVERLAPPING, B_LESS_OVERLAPPING} overlap_status_t;

  static const Q min(const Q& x, const Q& y) { return x<y? x : y; }
  static const Q max(const Q& x, const Q& y) { return x>y? x : y; }

 public:
  typedef Interval<Q>  interval_t;
  typedef std::list<interval_t> intervallist_t;
  typedef typename intervallist_t::const_iterator const_iterator;
  typedef typename intervallist_t::iterator       iterator;

  Intervals(){}
  Intervals(const intervallist_t& intervals) : intervallist_t(intervals){}

  const_iterator begin() const;
  const_iterator end()   const;
  iterator begin();
  iterator end();

  overlap_status_t overlapping(const Interval<Q>& A, const Interval<Q>& B) const;
  
  
  Intervals operator *(const Intervals& B) const;
  Intervals operator +(const Intervals& B) const;
  Intervals operator -(const Intervals& B) const;
  Intervals operator +=(const Intervals& B) { *this = (*this)+B; return *this; }
  Intervals operator -=(const Intervals& B) { *this = (*this)-B; return *this; }
  Intervals& operator *=(const Intervals& B){ *this = (*this)*B; return *this; }

  friend std::ostream& operator<<(std::ostream& s, const Intervals<Q>& Is){ 
    bool first = true;
    for(const_iterator i=Is.begin();i!=Is.end();i++){
      if(first) first = false;
      else s << " + ";
      s << *i;

    }

    return s;
  }
};

#endif
