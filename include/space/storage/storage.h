#ifndef LIBDISC_STORAGE_H
# define LIBDISC_STORAGE_H

#include <sys/types.h>

#include <complex>
#include <vector>
#include <space/storage/smallvector.h>
#include <space/storage/smallmatrix.h>

/* f : X -> Y, where X,Y are metric spaces with scalars in S */
template <class ContainerClass, class DF, typename Y = double, typename S = double> 
class DiscreteFunction {
 public: 
 typedef ContainerClass  container_type;
 container_type coefficients;

  /* Iterator interface to coefficients */
 typedef typename container_type::iterator       iterator;
 typedef typename container_type::const_iterator const_iterator;

 /* Constructors */
 DiscreteFunction() {}
 DiscreteFunction(const container_type &cs) : coefficients(cs) {}

 /* Abstract interface */
 virtual DF  operator + (const DF& ) const = 0;
 virtual DF  operator - (const DF& ) const = 0;

 virtual DF& operator +=(const DF& ) = 0;
 virtual DF& operator -=(const DF& ) = 0;

 virtual DF  operator * (const S& ) const = 0;
 virtual DF& operator *=(const S& ) = 0;
 
 virtual DF& operator =(const S& s) = 0;

 /* XXX: No obvious way of doing multiplication without integration for dofs */
 // virtual DF& operator * (const DF& ) const = 0;
 // virtual DF& operator *=(const DF& ) = 0;

 /* XXX: No obvious way of doing pointwise addition without integration for dofs */
 /*  virtual DF& operator +=(const Y& ) = 0;  */
 /*  virtual DF& operator -=(const Y& ) = 0; */
 /*  virtual DF& operator + (const Y& ) const = 0; */
 /*  virtual DF& operator - (const Y& ) const = 0; */

 const_iterator begin() const { return coefficients.begin(); }
 const_iterator end()   const { return coefficients.end();   }
 iterator begin()             { return coefficients.begin(); }
 iterator end()               { return coefficients.end();   }
 virtual void resize(const size_t n) { coefficients.resize(n); } 
 size_t size() const           { return coefficients.size();  }

 /* Random access interface to coefficients. Perhaps this isn't necessary. */
 Y& operator[](const size_t i)       { return coefficients[i]; }
 Y  operator[](const size_t i) const { return coefficients[i]; }

 virtual void update(const double epsilon=0) { }

 virtual size_t number_nonzero()   const { 
   size_t nonzero = 0;
   for(size_t i=0;i<coefficients.size();i++) if(fabs(coefficients[i])>0) ++nonzero;
   return nonzero;
 }
 virtual size_t number_intervals() const { return 1; }
};


#include <space/storage/simple.h>
#include <space/storage/pfinterval.h>

 

#endif
