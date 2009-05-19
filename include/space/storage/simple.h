#ifndef LIBDISC_SIMPLE_STORAGE_H
# define LIBDISC_SIMPLE_STORAGE_H

/* TODO: Split in .h and .cc */
/* TODO: Maybe add complex DiscreteFunction. */
/* TODO: Maybe add vector valued DiscreteFunction. */

#define pointwise_operator_vector(op,code)					        \
  SelfType& operator op ## = (const SelfType& g) {                                 \
    typename ContainerClass::const_iterator cg = g.coefficients.begin();		\
    for(typename ContainerClass::iterator c = coefficients.begin();c != coefficients.end(); c++, cg++)\
      { const Y &rhs = *cg; code; }					\
    return *this; \
  }               \
                  \
  SelfType& operator op (const SelfType& g) const { \
    SelfType &f = *new SelfType(coefficients);      \
    typename ContainerClass::const_iterator cg = g.coefficients.begin();\
    for(typename ContainerClass::iterator c = f.coefficients.begin();c != f.coefficients.end(); c++,cg++) \
      { const Y &rhs = *cg; code; }					\
    return f; \
  }

#define pointwise_operator_scalar(op,code,S)\
  SelfType& operator op ## =(const S& rhs){ \
    for(typename ContainerClass::iterator c = coefficients.begin();c != coefficients.end(); c++)\
      code; \
    return *this; \
  } \
  SelfType& operator op (const S& rhs) const { \
    SelfType& f = *new SelfType(coefficients);\
    for(typename ContainerClass::iterator c = f.coefficients.begin();c != f.coefficients.end(); c++)\
      code;   \
    return f; \
  }
#define pointwise_operator(op,code,S) pointwise_operator_vector(op,code) pointwise_operator_scalar(op,code,S)


#include <iostream>
#include <vector>
#include <string>
#include <netcdfcpp.h>

using namespace std;

template  <typename ContainerClass=std::vector<double>, typename Y = double, typename S = double> 
  class PointFunction_Simple : public DiscreteFunction<ContainerClass,PointFunction_Simple<ContainerClass,Y,S>,
  Y,S> 
{
 private:
 typedef DiscreteFunction<ContainerClass,PointFunction_Simple<ContainerClass,Y,S>,Y,S> BaseType;
 typedef PointFunction_Simple<ContainerClass,Y,S> SelfType;
 public:
  using BaseType::coefficients;
  /* Construction */
  PointFunction_Simple() {}
  PointFunction_Simple(const ContainerClass &cs) : BaseType(cs) {} 

  /* Point wise operations */
 pointwise_operator(+, (*c) += rhs,Y);
 pointwise_operator(-, (*c) -= rhs,Y);
 pointwise_operator(*, (*c) *= rhs,S);

};



template  <typename ContainerClass=std::vector<double>, typename Y = double, typename S = double> 
  class FEFunction_Simple : public DiscreteFunction<ContainerClass,FEFunction_Simple<ContainerClass,Y,S>,
  Y,S> 
{
 private:
 typedef DiscreteFunction<ContainerClass,FEFunction_Simple<ContainerClass,Y,S>,Y,S> BaseType;
 typedef FEFunction_Simple<ContainerClass,Y,S> SelfType;
 public:
  using BaseType::coefficients;
 /* Construction */
 FEFunction_Simple() {}
 FEFunction_Simple(const ContainerClass &cs) : BaseType(cs) {}
 
 /* FE wise operations */
 pointwise_operator_vector(+, (*c) += rhs);
 pointwise_operator_vector(-, (*c) -= rhs);
 pointwise_operator_scalar(*, (*c) *= rhs,S);
};



#endif
