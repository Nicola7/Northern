#include <iostream>
#include <cmath>
#include<iomanip>
double SDR (double R,double z);
double SDZ (double R,double z);
double stormerR(double R,double z,double v,double h);
double stormerZ(double R,double z,double y,double h);
double stormerv(double R,double z,double v,double y,double h);
double stormery(double R,double z,double v,double y,double h);
int main (void){
  std::cout<<std::fixed;std::setprecision(7);
  std::cout<<"s"<<'\t'<<"R"<<'\t'<<"z"<<'\n';
  double h=0.001;  int n=(0.3)/h;
  double Ri=0.257453;  double Zi=0.314687;
  std::cout<<"0"<<'\t'<<Ri<<'\t'<<Zi<<'\n';
  double ri=std::hypot(Ri,Zi);
  double Qi=(((-1)/Ri)+(Ri/(std::pow(ri,3))));
  Qi*=Qi;Qi=1-Qi;
  double Rpi=((std::sqrt(Qi))*(std::cos((5*(M_PI))/4)));
  double Zpi=((std::sqrt(Qi))*(std::sin((5*(M_PI))/4)));
  double R=stormerR(Ri,Zi,Rpi,h);
  double Z=stormerZ(Ri,Zi,Zpi,h);
  std::cout<<h<<'\t'<<R<<'\t'<<Z<<'\n';
  double v=stormerv(Ri,Zi,Rpi,Zpi,h);
  double y=stormery(Ri,Zi,Rpi,Zpi,h);
  /*todo lo que este arriba de este comentario es para calcular el valor de las condiciones iniciales y el "for" es para calcular los valores de la serie*/
  for (int ii=0;ii<=(n-2);++ii){
    double r,z,u,i;
    r=R;    z=Z;    u=v;    i=y;
    R=stormerR(r,z,u,h);
    Z=stormerZ(r,z,u,h);
    v=stormerv(r,z,u,i,h);
    y=stormery(r,z,u,i,h);
    std::cout<<h*(ii+2)<<'\t'<<R<<'\t'<<Z<<'\n';
  }
  return 0;
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
