#ifndef SMALLVECTOR_H
# define SMALLVECTOR_H

/** Light weight class for small vectors:
    - There's no reason to use heavy Teuchos vectors for points -- you'll never want to distribute them anyway.
    - TVectors lack all the useful spatial operators and are awkward to work with.
*/
#include <string.h>
#include <iostream>

template <int dim, typename X = double, typename S = double> class SmallVector {
 public:
  X x[dim];

  SmallVector() { memset(x,0,sizeof(x)); }	
  SmallVector(const SmallVector& y) { for(size_t i=0;i<dim;i++) x[i] = y.x[i]; }
  SmallVector(const X *x_) { for(size_t i=0;i<dim;i++) x[i] = x_[i]; }
  SmallVector& operator += (const SmallVector& y) {     
    for(size_t i=0;i<dim;i++) x[i] += y.x[i];
    return *this;
  }
  SmallVector& operator -= (const SmallVector& y) {     
    for(size_t i=0;i<dim;i++) x[i] -= y.x[i];
    return *this;
  }
  SmallVector& operator *=(const SmallVector& y) {
    for(size_t i=0;i<dim;i++) x[i] *= y.x[i];
    return *this;
  }

  SmallVector& operator /=(const SmallVector& y) {
    for(size_t i=0;i<dim;i++) x[i] /= y.x[i];
    return *this;
  }
  SmallVector& operator *=(const S& lambda) {
    for(size_t i=0;i<dim;i++) x[i] *= lambda;
    return *this;
  }
  SmallVector& operator +=(const S& lambda) {
    for(size_t i=0;i<dim;i++) x[i] += lambda;
    return *this;
  }
  SmallVector& operator -=(const S& lambda) {
    for(size_t i=0;i<dim;i++) x[i] -= lambda;
    return *this;
  }
  SmallVector& operator /=(const S& lambda) {
    for(size_t i=0;i<dim;i++) x[i] /= lambda;
    return *this;
  }

 bool operator<(const SmallVector& Y) const {
   for(size_t i=0;i<dim;i++) {
     if     (x[i] < Y.x[i]) return true;
     else if(x[i] > Y.x[i]) return false;
   }
   return false;
 }

 bool operator==(const SmallVector& Y) const {
   for(size_t i=0;i<dim;i++) if(x[i] != Y.x[i]) return false;
   return true;
 }

  SmallVector operator - (const SmallVector& y) const {
    SmallVector r(x);
    r -= y;
    return r;
  }
  SmallVector operator + (const SmallVector& y) const {
    SmallVector r(x);
    r += y;
    return r;
  }

  SmallVector operator * (const X& lambda) const {
    SmallVector r(x);
    r *= lambda;
    return r;
  }
  SmallVector operator + (const X& lambda) const {
    SmallVector r(x);
    r += lambda;
    return r;
  }
  SmallVector operator - (const X& lambda) const {
    SmallVector r(x);
    r -= lambda;
    return r;
  }
  SmallVector operator / (const X& lambda) const {
    SmallVector r(x);
    r /= lambda;
    return r;
  }
  SmallVector operator * (const SmallVector& y) const {
    SmallVector r(x);
    r *= y;
    return r;
  }
  SmallVector operator / (const SmallVector& y) const {
    SmallVector r(x);
    r /= y;
    return r;
  }

  X dot(const SmallVector& y) const {
    double sum = 0;
    for(off_t i=0;i<dim;i++) sum += x[i]*y.x[i];
    return sum;
  }

  X norm() const { return sqrt(dot(*this)); }

 X operator[](const size_t i) const { return x[i]; }
 X& operator[](const size_t i) { return x[i]; }

 friend std::ostream& operator<<(std::ostream& F, const SmallVector& v){
    F << "[ ";
    for(size_t i=0;i<dim;i++) F << v.x[i] << " ";
    F << "]";

    return F;
 } 
};

#endif
