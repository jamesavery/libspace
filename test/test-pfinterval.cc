#include <space/storage/interval.h>
#include <space/storage/pfinterval.h>

int main()
{
  typedef Interval<unsigned int>  interval;
  typedef Intervals<unsigned int> intervals;

  interval A(1,10), B(2,15);
  interval As_[3] = {interval(3,7), interval(9,12), interval(13,17)};
  interval Bs_[5] = {interval(1,2), interval(3,5), interval(6,9), interval(14,15), interval(16,17)};

  intervals As(list<interval>(&As_[0],&As_[3]));
  intervals Bs(list<interval>(&Bs_[0],&Bs_[5]));
  
  cout << A << "*"<<B<<" = " <<(B*A) << endl;
  cout << "As = " << As << endl;
  cout << "Bs = " << Bs << endl;
  cout << "As*Bs = " << (As*Bs) << endl;
  cout << "As+Bs = " << (As+Bs) << endl;

  double X_[] = {0,0,0,3,4,5,0,0,0,0,10,0,12,13,0,0,0};
  double Y_[] = {0,1,2,3,4,5,6,0,8,9,0,0,12,0,14,15,16};
  PFInterval X(vector<double>(&X_[0],&X_[17]));
  PFInterval Y(vector<double>(&Y_[0],&Y_[17]));

  cout << "X = "<<X<<endl;
  cout << "Y = "<<Y<<endl;
  cout << "X+Y = "<<X+Y<<endl;
  cout << "Y+X = "<<Y+X<<endl;
  cout << "X*Y = "<<X*Y<<endl;
  cout << "Y*X = "<<Y*X<<endl;

  return 0;
}
