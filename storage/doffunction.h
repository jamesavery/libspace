#ifndef LIBDISC_DOFFUNCTION_STORAGE_H
# define LIBDISC_DOFFUNCTION_STORAGE_H

/* TODO: Split in .h and .cc */
/* TODO: Maybe add complex DiscreteFunction. */
/* TODO: Maybe add vector valued DiscreteFunction. */


#include <iostream>
#include <string>
#include <netcdfcpp.h>

using namespace std;


#define pointwise_operator(op, code)\
  DiscreteFunction_DOFFunction& operator op ## = (const DiscreteFunction_DOFFunction& g) { \
    std::vector<double>::const_iterator cg = g.coefficients.begin();				  \
    for(std::vector<double>::iterator c = coefficients.begin();c != coefficients.end(); c++, cg++)\
      { const double &rhs = *cg; code; }					\
    return *this; \
  }               \
                  \
  DiscreteFunction_DOFFunction& operator op (const DiscreteFunction_DOFFunction& g) const { \
    DiscreteFunction_DOFFunction &f = *new DiscreteFunction_DOFFunction(coefficients);      \
    std::vector<double>::const_iterator cg = g.coefficients.begin();\
    for(std::vector<double>::iterator c = f.coefficients.begin();c != f.coefficients.end(); c++,cg++) \
      { const double &rhs = *cg; code; }					\
    return f; \
  }           \
              \
  DiscreteFunction_DOFFunction& operator op ## =(const double& rhs){ \
    for(std::vector<double>::iterator c = coefficients.begin();c != coefficients.end(); c++)\
      code; \
    return *this; \
  } \
  DiscreteFunction_DOFFunction& operator op (const double& rhs) const { \
    DiscreteFunction_DOFFunction &f = *new DiscreteFunction_DOFFunction(coefficients);\
    for(std::vector<double>::iterator c = f.coefficients.begin();c != f.coefficients.end(); c++)\
      code;   \
    return f; \
  }


class DOFFunction : public DiscreteFunction<std::vector<double>,DOFFunction> {
 private:
  typedef DiscreteFunction<std::vector<double>,DOFFunction> Base;
 public:
  /* Construction */
  DOFFunction() {}
  DOFFunction(const std::vector<double> &cs) : Base(cs) {}

  /* Point wise operations */
  pointwise_operator(+, (*c) += rhs);
  pointwise_operator(-, (*c) -= rhs);
  pointwise_operator(*, (*c) *= rhs);

  /* Random access interface to coefficients */
  double& operator[](const size_t i)       { return coefficients[i]; }


  DOFFunction(NcFile &f, const string varname)        { load(f,varname);    }
  DOFFunction(const string path, const string varname){ load(path,varname); }

  bool write(NcFile &f, const string varname){
    /* TODO: Write and read multidimensional NetCDF files where dimensions are taken from discretization. */
    NcDim *index_dim =  f.add_dim("index",coefficients.size());
    NcVar *data = f.add_var(varname.c_str(), ncDouble,index_dim);
    data->put(&coefficients[0], coefficients.size());
    return true;
  }

  bool write(const string path, const string varname){
    NcFile file(path.c_str(),NcFile::Replace);
    return write(file,varname);
  }

};


#endif
