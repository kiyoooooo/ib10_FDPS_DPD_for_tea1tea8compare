#include<iostream>
#include<fstream>
#include<unistd.h>
#include<sys/stat.h>
#include <stdio.h>
#include<vector>
#include<string.h>//文字列の代入に使う
struct data{ 
  char name[10];
  char type[10];
  int id;
  //  char type[5];
  int idsub;
  double x;
  double y;
  double z;
};

int main(){

  char title[40];
  int num,step;
  double boxdhx,boxdhy,boxdhz;

  //open reading file                                                                           
  FILE *fpi0;
  if((fpi0=fopen("nextplace.gro","r"))==NULL){
    printf("fpi0ERROR\n");
    return 0;
  }
  //open writing file
  FILE*fpo0;
  if((fpo0=fopen("sortnextplace.gro","w"))==NULL){
    printf("fpo0ERROR\n");
    return 0;
  }

  while(fscanf(fpi0,"%d\n%d",&step,&num) !=EOF){
    data *aa;
    int initidsub,initid;
    char initname[5],inittype[5];
    double initx,inity,initz;
    aa = new data[num];
    fprintf(fpo0,"%d\n%d\n",step,num);

    //std::cout<<"test0"<<std::endl;

    for(int i = 0; i < num; i++){
      //      fscanf(fpi0,"%5d-%5s%5s%5d%8.3f%8.3f%8.3f",&idsub,&name[i],&type[i],&id[i],&x[i],&y[i],&z[i]);
      //std::cout<<i<<std::endl;
      fscanf(fpi0,"%d %s %s %d %lf %lf %lf",	\
             &initidsub,			\
             &initname,				\
             &inittype,				\
             &initid,				\
             &initx,				\
             &inity,				\
             &initz);

      // aa[17839].idsub=5;
      //      std::cout<<"id = "<<initidsub<<std::endl;
      //aa[17840].idsub=initid;
      //      std::cout<<"ok1"<<std::endl;
      aa[initidsub].idsub=initidsub;
      //std::cout<<"ok2"<<std::endl;
      strcpy(aa[initidsub].name,initname);
      strcpy(aa[initidsub].type,inittype);
      aa[initidsub].id=initidsub;
      aa[initidsub].x=initx;
      aa[initidsub].y=inity;
      aa[initidsub].z=initz;
    }
    fscanf(fpi0,"%lf %lf %lf",&boxdhx,&boxdhy,&boxdhz);
    //        std::cout<<boxdhz<<std::endl;
    for(int i = 0; i < num; i++){
      //std::cout<<"i="<<i<<std::endl;
      //      std::cout<<"           aa[i].id = "<<aa[i].type<<std::endl;
      //            std::cout<<aa[i].name<<std::endl;
      fprintf(fpo0,"%9d%-5s%5s%9d%8.3f%8.3f%8.3f\n",\
	      aa[i].idsub,\
	      aa[i].name,\
	      aa[i].type,\
	      aa[i].id,\
	      aa[i].x,\
	      aa[i].y,\
	      aa[i].z);
      /*      fprintf(fpo0,"%5d",aa[i].idsub);
      fprintf(fpo0,"%-5s",aa[i].name);
      fprintf(fpo0,"%5s",aa[i].type);
      fprintf(fpo0,"%5d",aa[i].id);
      fprintf(fpo0,"%8.3f",aa[i].x);
      fprintf(fpo0,"%8.3f",aa[i].y);
      fprintf(fpo0,"%8.3f\n",aa[i].z);*/
    }
    fprintf(fpo0,"%f %f %f\n",boxdhx,boxdhy,boxdhz);
  }



  return 0;
}
