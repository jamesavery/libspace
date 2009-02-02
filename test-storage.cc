#include <stdio.h>
#include <sys/types.h>

#include <libspace/storage/smallvector.h>
#include <libspace/storage/smallmatrix.h>
#include <libspace/storage/storage.h>
#include <math.h>
#include <complex>
#include <iostream>


using namespace std;

int main()
{
  double a_[3]  = {1,2,3};
  double b_[3]  = {3,2,1};
  complex<double> c_[3] = {complex<double>(3,1),2,1};
  double A_[9] = {1,2,3,
		  4,5,6,
		  7,8,9};

  typedef SmallVector<3>  R3;
  typedef SmallVector<3,std::complex<double> > C3;
  typedef SmallMatrix<double,double> RMxN;
  R3 a(a_), b(b_);
  C3 c(c_);
  RMxN A(3,3,A_);
  

  cout << (b*2.0) << endl;
  cout << (b/2.0) << endl;
  cout << (c/=c) << endl;
  cout << A+A << endl;


  return 0;
}
