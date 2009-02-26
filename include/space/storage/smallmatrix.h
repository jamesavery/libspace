#ifndef SMALLMATRIX_H
# define SMALLMATRIX_H

/** Light weight class for small matrices:
    - There's no reason to use heavy Teuchos vectors for points -- you'll never want to distribute them anyway.
    - TVectors lack all the useful spatial operators and are awkward to work with.
*/
#include <string.h>
#include <assert.h>
#include <iostream>

template <typename X, typename S> class SmallMatrix {
 public:
  X *x;
  const size_t m,n, dim;

   SmallMatrix(const size_t m=1, const size_t n=1, const X *y=0) : m(m), n(n), dim(m*n) { 
    x = new X[dim]; 
    if(y != 0) memcpy(x,y,dim*sizeof(X));
    else       memset(x,0,dim*sizeof(X));
  }
  ~SmallMatrix(){ delete[] x; }

  SmallMatrix& operator += (const SmallMatrix& y) {     
    assert(m == y.m && n == y.n);
    for(size_t i=0;i<dim;i++) x[i] += y.x[i];
    return *this;
  }
  SmallMatrix& operator -= (const SmallMatrix& y) {     
    assert(m == y.m && n == y.n);
    for(size_t i=0;i<dim;i++) x[i] -= y.x[i];
    return *this;
  }
  SmallMatrix& operator /=(const SmallMatrix& y) {
    assert(m == y.m && n == y.n);
    for(size_t i=0;i<dim;i++) x[i] /= y.x[i];
    return *this;
  }
  SmallMatrix& operator *=(const S& lambda) {
    for(size_t i=0;i<dim;i++) x[i] *= lambda;
    return *this;
  }
  SmallMatrix& operator *=(const SmallMatrix& y) {
    assert(n == y.m);
    X *r = new X[m*y.n];
    for(size_t i=0;i<m;i++)
      for(size_t j=0;j<y.n;j++){
	X rij(0);
	for(size_t k=0;k<n;k++) rij += x[i*m+k]*y.x[k*n+j];
	r[i*m+j] = rij;
      }

    delete[] x;
    x = r;
    return *this;
  }

  SmallMatrix operator - (const SmallMatrix& y) const {
    SmallMatrix r(m,n,x);
    r -= y;
    return r;
  }
  SmallMatrix operator + (const SmallMatrix& y) const {
    SmallMatrix r(m,n,x);
    r += y;
    return r;
  }
  SmallMatrix operator / (const SmallMatrix& y) const {
    SmallMatrix r(x);
    r /= y;
    return r;
  }

  SmallMatrix operator * (const SmallMatrix& y) const {
    assert(n == y.m);
    SmallMatrix r(m,y.n);
    for(size_t i=0;i<m;i++)
      for(size_t j=0;j<y.n;j++){
	X rij(0);
	for(size_t k=0;k<n;k++) rij += x[i*m+k]*y.x[k*n+j];
	r.x[i*m+j] = rij;
      }
    return r;
  }


  SmallMatrix operator * (const S& lambda) const {
    SmallMatrix r(x);
    r *= lambda;
    return r;
  }

  friend std::ostream& operator<<(std::ostream& F, const SmallMatrix& v){
    for(size_t i=0;i<v.m;i++){
      F << "[ ";
      for(size_t j=0;j<v.m;j++) F << v.x[i*v.m+j] << " ";
      F << "]\n";
    }
    return F;
  }
};

#endif
