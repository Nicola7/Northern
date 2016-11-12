#include <iostream>
#include <cmath>
#include<iomanip>
double SDR (double R,double z);
double SDZ (double R,double z);
double stormerR(double R,double z,double v,double h);
double stormerv(double R,double z,double v,double y,double h);
double stormerZ(double R,double z,double y,double h);
double stormery(double R,double z,double v,double y,double h);
int main (void){
  const double h=0.1;
  double R=stormerR(1,1,1,h);
  double Z=stormerZ(1,1,1,h);
  double v=stormerv(1,1,1,1,h);
  double y=stormery(1,1,1,1,h);
  for (int ii=0;ii<3;++ii){
    double r,z,u,i;
    r=R;    z=Z;    u=v;    i=y;
    R=stormerR(r,z,u,h);
    Z=stormerZ(r,z,u,h);
    v=stormerv(r,z,u,v,h);
    y=stormery(r,z,u,v,h);
    std::cout<<R<<'\t'<<Z<<'\n';
  }
  return 0;
}
double SDR(double R,double z){
  double y=(-0.5);
  double r=std::hypot(R,z);
  r*=r;
  double t=(((2*y)/R)+(R/(std::pow(r,3))));
  t*=(((2*y)/(R*R))+((3*(R*R))/(std::pow(r,5)))-(1.0/(std::pow(r,3))));
  return t;
}
double SDZ (double R,double z){
  double y=(-0.5);
  double r=std::hypot(R,z);
  r*=r;
  double t=((3*R*z)/(std::pow(r,5)));
  t*=(((2*y)/R)+((R)/(std::pow(r,3))));
  return t;
}
double stormerR(double R,double z,double v,double h){
  double n=((((SDR(R,z))*(h*0.5))+v)*h);
  n+=R;
  return n;
}
double stormerZ(double R,double z,double y,double h){
  double n=((((SDZ(R,z))*(h*0.5))+y)*h);
  n+=R;
  return n;

}
double stormerv(double R,double z,double v,double y,double h){
  double n=stormerR(R,z,v,h);
  double t=stormerZ(R,z,y,h);
  n=((SDR(n,z))*(h*0.5));
  n+=((h*0.5)*(SDR(R,z)));
  n+=v;
  return n;
}
double stormery(double R,double z,double v,double y,double h){
  double n=stormerR(R,z,v,h);
  double t=stormerZ(R,z,y,h);
  n=((SDZ(n,z))*(h*0.5));
  n+=((h*0.5)*(SDZ(R,z)));
  n+=y;
  return n;
}
