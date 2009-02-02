
template <class ContainerClass, class subDF, typename Q = double> class DiscreteFunction { 
 public:
  typedef ContainerClass container_type;
  container_type coefficients;

  /* Iterator interface to coefficients */
  typedef typename container_type::iterator       iterator;
  typedef typename container_type::const_iterator const_iterator;

  /* Constructors */
  DiscreteFunction() {}
  DiscreteFunction(const typename DiscreteFunction::container_type &cs) : coefficients(cs) {}

  /* Abstract interface */
  virtual DiscreteFunction& operator + (const subDF& ) const = 0;
  virtual DiscreteFunction& operator - (const subDF& ) const = 0;
  virtual DiscreteFunction& operator * (const subDF& ) const = 0;
  virtual DiscreteFunction& operator + (const Q& ) const = 0;
  virtual DiscreteFunction& operator - (const Q& ) const = 0;
  virtual DiscreteFunction& operator * (const Q& ) const = 0;

  virtual DiscreteFunction& operator +=(const subDF& ) = 0;
  virtual DiscreteFunction& operator -=(const subDF& ) = 0;
  virtual DiscreteFunction& operator *=(const subDF& ) = 0;
  virtual DiscreteFunction& operator +=(const Q& ) = 0;
  virtual DiscreteFunction& operator -=(const Q& ) = 0;
  virtual DiscreteFunction& operator *=(const Q& ) = 0;

  virtual Q& operator[](const size_t) = 0;
};

