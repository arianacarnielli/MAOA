#include<string>
#include<vector>

#ifndef PRP_
#define PRP_

using namespace std;

class PRP{
 public:

  int n;     // number of clients
  int l;     // number of time periods

  double u;  // unit production cost
  double f;   // fixed production setup cost
  double C;   // production capacity
  double Q;   // vehicle capacity
  int m;     // number of vehicles

  vector<pair<double,double> > xy; // clients coordinates
  vector<double> h;  // h_i unit inventory cost at node i
  vector<vector<double> > d; // d_it demand at customer i in period t
  vector<double> L; // maximum inventory level at node i
  vector<double> L0; //initial inventory at node i

  int dist;  // 1 if the instance is from Archetti et al. 
             // 2 if the instance is from Boudia et al.,
  // if dist==1
      // The transportation cost c_ij is 
      // INT((SQRT((X(I)-X(J))**2+(Y(I)-Y(J))**2))+.5)
  // if dist ==2:
       double mc;   // kilometric cost
       // The transportation cost c_ij is mc times the euclidian cost
  

  PRP(int n, int l);
  PRP(istream&);
  void write_screen_txt();
  void write_instance_file(string);
  double cost(int i, int j);

};

#endif