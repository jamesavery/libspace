#ifndef LIBDISC_SIMPLE_STORAGE_H
# define LIBDISC_SIMPLE_STORAGE_H

/* TODO: Split in .h and .cc */
/* TODO: Maybe add complex DiscreteFunction. */
/* TODO: Maybe add vector valued DiscreteFunction. */


#include <iostream>
#include <string>
#include <netcdfcpp.h>

using namespace std;

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

  /* Storage management */
  void resize(const size_t n){
    coefficients.resize(n);
  }

  /* NetCDF interface. Perhaps encapsulate to accomodate NetCDF, HDF5, VTK, etc. */
  bool load(NcFile &f, const string varname){
    if(!f.is_valid()){
      cerr << "Cannot read "<<varname<<": NetCDF-file is invalid." << endl;
      return false;
    } else {
      NcVar *data = f.get_var(varname.c_str());
      if(!data){ 
	cerr << "Cannot read "<<varname<<": NetCDF-file is valid, but variable was not found." << endl;
	return false;
      } else {
	NcDim *dim  = data->get_dim(0);
	coefficients.resize(dim->size());
	data->get(&coefficients[0],dim->size());
      }
    }
    return true;
  }
  bool load(const string path, const string varname){
    NcFile file(path.c_str(),NcFile::ReadOnly);
    return load(file,varname);
  }
  DiscreteFunction_Simple(NcFile &f, const string varname)        { load(f,varname);    }
  DiscreteFunction_Simple(const string path, const string varname){ load(path,varname); }

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
