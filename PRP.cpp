#include<iostream>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<cmath>
#include "PRP.h"

using namespace std;


PRP::PRP(int nn, int ll){
  int i;

  n=nn;
  l=ll;

  xy.resize(n+1);
  h.resize(n+1);
  d.resize(n+1);
  for (i=1;i<=n;i++) d[i].resize(l);
  L.resize(n+1);
  L0.resize(n+1);

}

void PRP::write_screen_txt(){
  int i,t;
  if (dist==1)
    cout<<"Instance A"<<endl;
  else
    cout<<"Instance B"<<endl;

  cout<<"n="<<n<<" l="<<l<<endl;

  cout<<"u="<<u<<endl;
  cout<<"f="<<f<<endl;
  cout<<"C="<<C<<endl;
  cout<<"Q="<<Q<<endl;
  cout<<"m="<<m<<endl;

  if (dist==2) cout<<"mc="<<mc<<endl;

  cout<<"Depot 0 ("<<xy[0].first<<","<<xy[0].second<<") : ";
  cout<<"h[0]="<<h[0];
  cout<<"  L[0]="<<L[0];
  cout<<"  L0[0]="<<L0[0]<<endl;

  for (i=1;i<=n;i++){
    cout<<"Client "<<i<<" ("<<xy[i].first<<","<<xy[i].second<<") : ";
    cout<<" h["<<i<<"]="<<h[i];
    cout<<"  L["<<i<<"]="<<L[i];
    cout<<"  L0["<<i<<"]="<<L0[i];
    cout<<endl<<"d"<<i<<"= [";
    for (t=0;t<l-1;t++)
      cout<<d[i][t]<<", ";
    cout<<d[i][l-1]<<" ]"<<endl;
  }

}


void PRP::write_instance_file(string inst){
  ostringstream name;
  int i,t;

  name.str("");
  if (dist==1)   name<<"A_";
  else    name<<"B_";

  name<<setfill('0')<<setw(3)<<n<<"_"<<inst<<".prp";

  ofstream fic(name.str().c_str());

  fic<<"Type "<<dist<<endl;
  fic<<"n "<<n<<endl;
  fic<<"l "<<l<<endl;
  fic<<"u "<<u<<endl;
  fic<<"f "<<f<<endl;
  fic<<"C "<<C<<endl;
  fic<<"Q "<<Q<<endl;
  fic<<"m "<<m<<endl;
  if (dist==2) fic<<"mc "<<mc<<endl;
  for (i=0;i<=n;i++){
    fic<<i<<" ";
    fic<<xy[i].first<<" "<<xy[i].second<<" : ";
    fic<<"h "<< h[i];
    fic<<" L "<<L[i];
    fic<<" L0 "<<L0[i];
    fic<<endl;
  }
  fic<<"d"<<endl;
  for (i=1;i<=n;i++){
    fic<<i<<" ";
    for(t=0;t<l;t++)
      fic<<d[i][t]<<" ";
    fic<<endl;
  }


  fic.close();

}

PRP::PRP(istream &fic){
  string r;
  int i,t;

  fic>>r;  fic>>dist;
  fic>>r;  fic>>n;
  fic>>r; fic>>l;;

  xy.resize(n+1);
  h.resize(n+1);
  d.resize(n+1);
  for (i=1;i<=n;i++) d[i].resize(l);
  L.resize(n+1);
  L0.resize(n+1);

  fic>>r;  fic>>u;
  fic>>r;  fic>>f;
  fic>>r;  fic>>C;
  fic>>r;  fic>>Q;
  fic>>r;  fic>>m;
  if (dist==2) {fic>>r;fic>>mc;}
  for (i=0;i<=n;i++){
    fic>>r;
    fic>>xy[i].first; fic>>xy[i].second;
    fic>>r; fic>>r;
    fic>>h[i];
    fic>>r;fic>>L[i];
    fic>>r;fic>>L0[i];
  }
  fic>>r;
  for (i=1;i<=n;i++){
    fic>>r;
    for(t=0;t<l;t++)
      fic>>d[i][t];
  }

}

double PRP::cost(int i, int j) {
// if dist==1 the transportation cost c_ij is 
// INT((SQRT((X(I)-X(J))**2+(Y(I)-Y(J))**2))+.5)
// if dist ==2 the transportation cost c_ij is mc times the euclidian cost
    double euclid = sqrt(pow(xy[i].first - xy[j].first, 2) + pow(xy[i].second - xy[j].second, 2));
    if (dist == 2) {
        return mc * euclid;
    }
    else {
        return (double)((int)(euclid + 0.5));
    }
}