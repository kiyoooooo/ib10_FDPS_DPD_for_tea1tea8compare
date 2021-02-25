#include<iostream>
#include<fstream>
#include<unistd.h>
#include<sys/stat.h>
#include <stdio.h>
#include<vector>
#include<string.h>//文字列の代入に使う 



struct vec3{
  double x,y,z;

  vec3() : x(0), y(0), z(0) {}

  vec3& operator +=(const vec3 &bb){
    x += bb.x;
    y += bb.y;
    z += bb.z;
    return *this;
}

  vec3 operator /(const double &bb){
    vec3 v = *this;
    v.x /= bb;
    v.y /= bb;
    v.z /= bb;
    return v;
  }
};


struct data{
  int step;
  vec3 pressure;
};

int main(){

  
  std::vector<data> aa;
  //  data *aa = new data[line];
  
  vec3 temppres = vec3();

  //open reading file
  FILE *fpi0;
  if((fpi0=fopen("pressure.dat","r")) == NULL){
    printf("fpi0ERROR\n");
    return 0;
  }

  //open writing file
  /*  FILE *fpo0;
  if((fpo0=fopen("avepressure.dat","w"))=NULL){
    printf("fpo0ERROR\n");
    return 0;
    }*/
  while(fscanf(fpi0,"%lf %lf %lf\n",&temppres.x,&temppres.y,&temppres.z)!=EOF){
    aa.push_back(data());
    aa.back().pressure.x = temppres.x;
    aa.back().pressure.y = temppres.y;
    aa.back().pressure.z = temppres.z;
  }

  vec3 sumpressure;
  for(int i=30000; i<50000; i++){
    sumpressure += aa[i].pressure;
  }
  vec3 avepressure = sumpressure / 20000;

  std::cout<< avepressure.x <<" "<< avepressure.y <<" "<< avepressure.z <<std::endl;
  return 0;
}
