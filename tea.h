#include <stdio.h>
#include <stdlib.h>
#include <math.h>


unsigned long long teacode(double v){
union Uni
{
  double aa;
  //long long y;
  unsigned int bb[2];
  unsigned long long cc;
}uni;
//unsigned int sum=0,delta=0x9e3779b9,n=1;

  uni.aa=v;
  uni.cc = uni.cc << 32;
  //  unsigned int y=uni.bb[0];
  //uni.bb[0]=uni.bb[1];
  //uni.bb[1]=y;
  //y=uni.bb[0];
  //z=uni.bb[1];

  //printf("%.20lf %.20lf %16.llu %16.llu\n",uni.aa,v[1],y,z);

  /*while(n-->0){
    sum+=delta;
    //uni.bb[0]+=((uni.bb[1]<<4)+k[0])^(uni.bb[1]+sum)^((uni.bb[1]>>5)+k[1]);
    //uni.bb[1]+=((uni.bb[0]<<4)+k[2])^(uni.bb[0]+sum)^((uni.bb[0]>>5)+k[3]);
    uni.bb[0]+=((uni.bb[1]<<4)+3)^(uni.bb[1]+sum)^((uni.bb[1]>>5)+4);
    uni.bb[1]+=((uni.bb[0]<<4)+5)^(uni.bb[0]+sum)^((uni.bb[0]>>5)+6);
    }*/
  
  // g[0]=(double)v[0];
  //g[0]=g[0]/pow10;

  //uni.bb=y;
  //  v=(double)uni.cc;
  //uni.bb=z;
  //v[1]=(double)z;
  //printf("%lf\n",v);
  //return v; 
  return uni.cc;
}
