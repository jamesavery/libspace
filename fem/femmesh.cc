#include <math.h>
#include <libspace/mesh/femmesh.h>

/* Integration and inner products */
                                                                 
template <int dim>  double FEMMesh<dim>::Integrate(const RealValues& f) const 
{  /**< \int_\Omega f(x) dx */
  RealValues::const_iterator    c = f.coefficients.begin(), 
    end = f.coefficients.end(), w = weights.coefficients.begin();

  for(sum = 0; c != end; c++,w++) sum += (*c) * (*w);

  return sum;  
}


/** \f$\int_\Omega f(x)~dx\f$ */
template <int dim>  double FEMMesh<dim>::Integrate(const RealValues& f, const RealValues &g) const 
{  
  RealValues::const_iterator    fc = f.coefficients.begin(), 
    end = f.coefficients.end(), gc = g.coefficients.begin(),
    w = weights.coefficients.begin();

  for(sum = 0; c != end; fc++,gc++,w++) sum += (*fc) * (*gc) * (*w);

  return sum;  
}


