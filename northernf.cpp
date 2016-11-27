#include <iostream>
#include <cmath>
#include<iomanip>
double  rinicial (double Ri, double Zi);
double  Qinicial (double Ri, double ri);
double  Rpinicial (double Qi);
double  Zpinicial (double Qi);
double SDR (double R,double z, const double &l);
double SDZ (double R,double z,const double &l);
double stormerR(double R,double z,const double &l,double v,const double &h);
double stormerZ(double R,double z,const double &l,double y,const double &h);
double stormerv(double R,double z,const double &l,double v,double y,const double &h);
double stormery(double R,double z,const double &l,double v,double y,const double &h);
double phiprima (double R, double z);
double rungekutta(double Rb,double zb,double Rm,double zm,double Ra,double za,const double &h, const double &l); 
double cartex(double R,double phi);
double cartey(double R,double phi);
int main (void){
  std::cout<<std::fixed;std::setprecision(7);
  /*std::cout<<"s"<<'\t'<<'\t'<<"X"<<'\t'<<'\t'<<"Y"<<'\t'<<'\t'<<"Z"<<'\n';*/
  const double l=-0.5; double h=0.00001;  int n=(0.3)/h; double w=(2*h);
  double Ri=0.257453;  double Zi=0.314687;
  std::cout<<"0"<<'\t'<<'\t'<<cartex(Ri,0)<<'\t'<<cartey(Ri,0)<<'\t'<<Zi<<'\n';
  double ri=rinicial(Ri,Zi);
  double Qi=Qinicial(Ri,ri);
  double Rpi=Rpinicial(Qi);
  double Zpi=Zpinicial(Qi);
  double R=stormerR(Ri,Zi,l,Rpi,h);
  double Z=stormerZ(Ri,Zi,l,Zpi,h);
  double v=stormerv(Ri,Zi,l,Rpi,Zpi,h);
  double y=stormery(Ri,Zi,l,Rpi,Zpi,h);
  for (int ii=0;ii<=(n-2);++ii){
    double r,z,u,i,uu,kk;
    uu=r; kk=z;
    r=R;    z=Z;    u=v;    i=y;
    R=stormerR(r,z,l,u,h);
    Z=stormerZ(r,z,l,u,h);
    v=stormerv(r,z,l,u,i,h);
    y=stormery(r,z,l,u,i,h);
    if (ii==0){
      double phi=rungekutta(Ri,Zi,r,z,R,Z,w,l);
      std::cout<<h*(ii+2)<<'\t'<<cartex(R,phi)<<'\t'<<cartey(R,phi)<<'\t'<<Z<<'\n';
    }
    if (ii%2==0 && ii!=0){
      double phi=rungekutta(uu,kk,r,z,R,Z,w,l);
      std::cout<<h*(ii+2)<<'\t'<<cartex(R,phi)<<'\t'<<cartey(R,phi)<<'\t'<<Z<<'\n';
    }
  }
  return 0;
}
double  rinicial (double Ri, double Zi){
  return std::hypot(Ri,Zi);
}
double  Qinicial (double Ri, double ri){
  double Qi=(((-1)/Ri)+(Ri/(std::pow(ri,3))));
  Qi*=Qi;Qi=1-Qi;
  return Qi;
}
double  Rpinicial (double Qi){
  return ((std::sqrt(Qi))*(std::cos((5*(M_PI))/4)));
}
double Zpinicial (double Qi){
  return ((std::sqrt(Qi))*(std::sin((5*(M_PI))/4)));
}
double SDR(double R,double z,const double &l){
  double r=std::hypot(R,z);
  double t=(((2*l)/R)+(R/(std::pow(r,3))));
  t*=(((2*l)/(R*R))+((3*(R*R))/(std::pow(r,5)))-(1.0/(std::pow(r,3))));
  return t;
}
double SDZ (double R,double z,const double &l){
  double r=std::hypot(R,z);
  double t=((3*R*z)/(std::pow(r,5)));
  t*=(((2*l)/R)+((R)/(std::pow(r,3))));
  return t;
}
double stormerR(double R,double z,const double &l,double v,const double &h){
  double n=((((SDR(R,z,l))*(h*0.5))+v)*h);
  n+=R;
  return n;
}
double stormerZ(double R,double z,const double &l,double y,const double &h){
  double n=((((SDZ(R,z,l))*(h*0.5))+y)*h);
  n+=z;
  return n;

}
double stormerv(double R,double z,const double &l,double v,double y,const double &h){
  double n=stormerR(R,z,l,v,h);
  double t=stormerZ(R,z,l,y,h);
  n=((SDR(n,t,l))*(h*0.5));
  n+=((h*0.5)*(SDR(R,z,l)));
  n+=v;
  return n;
}
double stormery(double R,double z,const double &l,double v,double y,const double &h){
  double n=stormerR(R,z,l,v,h);
  double t=stormerZ(R,z,l,y,h);
  n=((SDZ(n,t,l))*(h*0.5));
  n+=((h*0.5)*(SDZ(R,z,l)));
  n+=y;
  return n;
}
double phiprima (double R, double z,const double &l){
  double t=std::hypot(R,z);
  t=std::pow(t,3);
  t=(1/R)*(((2*l)/R)+(R/t));
  return t;
}
double rungekutta(double Rb,double zb,double Rm,double zm,double Ra,double za,const double &h,const double &l){
  double k1=phiprima(Rb,zb,l);
  double k2=phiprima(Rm,zm,l);
  double k3=phiprima(Rm,zm,l);
  double k4=phiprima(Ra,za,l);
  double y=k1+((h/6)*(k1+(2*k2)+(2*k3)+k4));
  return y;
}
double cartex(double R,double phi){
  return (R*(std::cos(phi)));
}
double cartey(double R,double phi){
  return (R*(std::sin(phi)));
}
