#ifndef LIBDISC_COMPACT_STORAGE_H
# define LIBDISC_COMPACT_STORAGE_H


template  <typename Space, typename ContainerClass=std::vector<double>, typename Y = double, typename S = double> 
  class PointFunction_Compact : public DiscreteFunction<ContainerClass,PointFunction_Simple<ContainerClass,Y,S>,
  Y,S> 
{
 private:
  typedef DiscreteFunction<ContainerClass,PointFunction_Simple<ContainerClass,Y,S>,Y,S> BaseType;
  typedef PointFunction_Simple<Space,ContainerClass,Y,S> SelfType;

 public:
  using BaseType::coefficients;
  Space& space;

 PointFunction_Compact(Space& space, const ContainerClass& cs) : BaseType(cs) 
 {
   /* 1. Determine bounding box */
 }
};

#endif
