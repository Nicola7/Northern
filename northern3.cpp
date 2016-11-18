#include <iostream>
#include <cmath>
#include<iomanip>
double  rinicial (double Ri, double Zi);
double  Qinicial (double Ri, double ri);
double  Rpinicial (double Qi);
double  Zpinicial (double Qi);
double SDR (double R,double z);
double SDZ (double R,double z);
double stormerR(double R,double z,double v,double h);
double stormerZ(double R,double z,double y,double h);
double stormerv(double R,double z,double v,double y,double h);
double stormery(double R,double z,double v,double y,double h);
double phiprima (double R, double z);
double rungekutta(double Rb,double zb,double Rm,double zm,double Ra,double za,double h); 
double cartex(double R,double phi);
double cartey(double R,double phi);
int main (void){
  std::cout<<std::fixed;std::setprecision(7);
  /*std::cout<<"s"<<'\t'<<'\t'<<"X"<<'\t'<<'\t'<<"Y"<<'\t'<<'\t'<<"Z"<<'\n';*/
  double h=0.0000001;  int n=(0.3)/h; double w=(2*h);
  double Ri=0.257453;  double Zi=0.314687;
  std::cout<<"0"<<'\t'<<'\t'<<cartex(Ri,0)<<'\t'<<cartey(Ri,0)<<'\t'<<Zi<<'\n';
  double ri=rinicial(Ri,Zi);
  double Qi=Qinicial(Ri,ri);
  double Rpi=Rpinicial(Qi);
  double Zpi=Zpinicial(Qi);
  double R=stormerR(Ri,Zi,Rpi,h);
  double Z=stormerZ(Ri,Zi,Zpi,h);
  double v=stormerv(Ri,Zi,Rpi,Zpi,h);
  double y=stormery(Ri,Zi,Rpi,Zpi,h);
  /*todo lo que este arriba de este comentario es para calcular el valor de las condiciones iniciales y el "for" es para calcular los valores de la serie*/
  for (int ii=0;ii<=(n-2);++ii){
    double r,z,u,i,uu,kk;
    uu=r; kk=z;
    r=R;    z=Z;    u=v;    i=y;
    R=stormerR(r,z,u,h);
    Z=stormerZ(r,z,u,h);
    v=stormerv(r,z,u,i,h);
    y=stormery(r,z,u,i,h);
    if (ii==0){
      double phi=rungekutta(Ri,Zi,r,z,R,Z,w);
      std::cout<<h*(ii+2)<<'\t'<<cartex(R,phi)<<'\t'<<cartey(R,phi)<<'\t'<<Z<<'\n';
    }
    if (ii%2==0 && ii!=0){
      double phi=rungekutta(uu,kk,r,z,R,Z,w);
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
/*la funcion SDR calcula el valor de la segunda derivada de R dados uno valores R y z*/
double SDR(double R,double z){
  double y=(-0.5);
  double r=std::hypot(R,z);
  double t=(((2*y)/R)+(R/(std::pow(r,3))));
  t*=(((2*y)/(R*R))+((3*(R*R))/(std::pow(r,5)))-(1.0/(std::pow(r,3))));
  return t;
}
/*la funcion SDZ calcula el valor de la segunda derivada de z dados ciertos valores de R y z*/
double SDZ (double R,double z){
  double y=(-0.5);
  double r=std::hypot(R,z);
  double t=((3*R*z)/(std::pow(r,5)));
  t*=(((2*y)/R)+((R)/(std::pow(r,3))));
  return t;
}
/* La funcion stormerR calcula el valor de R para s+h donde s es el valor en el que se evaluo la R y z que se tienen coo argumentos de la funcion*/
double stormerR(double R,double z,double v,double h){
  double n=((((SDR(R,z))*(h*0.5))+v)*h);
  n+=R;
  return n;
}
/*identico que la funcion anterior pero para z*/
double stormerZ(double R,double z,double y,double h){
  double n=((((SDZ(R,z))*(h*0.5))+y)*h);
  n+=z;
  return n;

}
/*esta funcion stormerv calcula el valor de la primera derivada de R para s+h donde s es el valor que R y Z que se tienen como argumentos de la funcion*/
double stormerv(double R,double z,double v,double y,double h){
  double n=stormerR(R,z,v,h);
  double t=stormerZ(R,z,y,h);
  n=((SDR(n,t))*(h*0.5));
  n+=((h*0.5)*(SDR(R,z)));
  n+=v;
  return n;
}
/*igual que la anterior funcion pero todo para z*/
double stormery(double R,double z,double v,double y,double h){
  double n=stormerR(R,z,v,h);
  double t=stormerZ(R,z,y,h);
  n=((SDZ(n,t))*(h*0.5));
  n+=((h*0.5)*(SDZ(R,z)));
  n+=y;
  return n;
}
double phiprima (double R, double z){
  double t=std::hypot(R,z);
  t=std::pow(t,3);
  t=(1/R)*(((-1)/R)+(R/t));
  return t;
}
double rungekutta(double Rb,double zb,double Rm,double zm,double Ra,double za,double h){
  double k1=phiprima(Rb,zb);
  double k2=phiprima(Rm,zm);
  double k3=phiprima(Rm,zm);
  double k4=phiprima(Ra,za);
  double y=k1+((h/6)*(k1+(2*k2)+(2*k3)+k4));
  return y;
}
double cartex(double R,double phi){
  return (R*(std::cos(phi)));
}
double cartey(double R,double phi){
  return (R*(std::sin(phi)));
}
