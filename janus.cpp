//
// vdwtest.cpp
//
// simple short-range MD test program
//
// Jun Makino March 9 2015
//
// Known problem as of Mar 13 2015
//   -- acc and phi are consistent only for the case of m=1
// This has been fixed as of Mar 15 2015

#include<iostream>
#include<fstream>
#include<unistd.h>
#include<sys/stat.h>
#include<particle_simulator.hpp>
#include <stdio.h>//add by kiyoshiro for FILE *fpo2
#ifdef DISSIPATIVE_RANDOM
#include "random_number.h"
#include "tea.h"//for32TEA(1)
#include<limits.h>//for32TEA(1)
#endif

#include <cmath>

#define LOWEST_LIMIT 1e-4
#define CALC_PRESSURE
#define UNI_SET
const PS::F64 dj = 0.3;

// patch informations
#if 1
const int Npatch = 3;
const PS::F64vec3 patch[Npatch] = {PS::F64vec3( 0.0,           1.0, 0.0),
				   PS::F64vec3( 0.86602540378,-0.5, 0.0),
				   PS::F64vec3(-0.86602540378,-0.5, 0.0)};
//const PS::F64 coef_r = 396.; // repulsive coefficient
const PS::F64 coef_r[7][7] = {{ 25.00, 24.41, 24.22, 29.55, 40.53, 40.57, 64.53},
			      { 24.41, 25.00, 22.97, 27.32, 42.85, 44.00, 66.40},
			      { 24.22, 22.97, 25.00, 28.89, 40.73, 43.84, 59.56},
                              { 29.55, 27.32, 28.89, 25.00, 31.51, 11.92, 56.82},
			      { 40.53, 42.85, 40.73, 31.51, 25.00, 33.34, 3.290},
			      { 40.57, 44.00, 43.84, 11.92, 33.34, 25.00, 44.62},
			      { 64.53, 66.40, 59.56, 56.82, 3.290, 44.62, 25.00}}; // repulsive coefficient
//const PS::F64 coef_a[Npatch][Npatch] ={{ 220., 220., 220.},
//				       { 220., 220., 220.},
//				       { 220., 220., 220.}}; // attractive coefficient of patch
//const PS::F64 coef_v = 0.5;  // exponent of f(n,n,r), nu

//const PS::F64 tm[Npatch] = {45.0 / 180.0 * M_PI,
//			    45.0 / 180.0 * M_PI,
//			    45.0 / 180.0 * M_PI}; // theta_m
#endif



#ifdef DISSIPATIVE_RANDOM
const PS::F64 gamma_dpd = 4.5;
const PS::F64 rn_max = sqrt(3.0);
#endif
//------------------


template <class T>
T P0(const T a){
  return a;
}
template <class T>
T P1(const T a){
  return T(-a.y, a.x, a.w,-a.z);
}
template <class T>
T P2(const T a){
  return T(-a.z,-a.w, a.x, a.y);
}
template <class T>
T P3(const T a){
  return T(-a.w, a.z,-a.y, a.x);
}

/*
template <class T>
PS::F64mat4 S(const T a){
  return PS::F64mat4(a.x,-a.y,-a.z,-a.w,
		     a.y, a.x, a.w,-a.z,
		     a.z,-a.w, a.x, a.y,
		     a.w, a.z,-a.y, a.x);
		     }*/

class BinaryHeader{
public:
  PS::S32 n_ptcl;
  BinaryHeader(const PS::S32 n_ptcl_) : n_ptcl(n_ptcl_){}

  void writeBinary(FILE* fp) const {
    fwrite(&n_ptcl,sizeof(PS::S32),1,fp);
  }
  PS::S32 readBinary(FILE* fp){
    fread(&n_ptcl,sizeof(PS::S32),1,fp);
    return n_ptcl;
  }
};

class CDVHeader{
public:
  PS::F64vec s,e;
  CDVHeader(){}
  CDVHeader(PS::F64vec _s,PS::F64vec _e) : s(_s),e(_e) {}
  PS::S32 readAscii(FILE * fp){
    return 0;
  }
  void writeAscii(FILE* fp) const{

    fprintf(fp, "'box_sx=%lf,box_sy=%lf,box_sz=%lf,box_ex=%lf,box_ey=%lf,box_ez=%lf\n",
	    s.x,s.y,s.z, e.x,e.y,e.z);

    fprintf(fp, "'r0=0.35\n");
    fprintf(fp, "'r1=0.20\n");
    fprintf(fp, "'r2=0.1 c2=(256,256,0)\n");
  }
};

class Force{
public:
  PS::F64vec3 force;
  PS::F64vec3 midforce;//修正速度ベルレ法で利用
  PS::F64vec3 torque;
  PS::F64 pot;
#ifdef CALC_PRESSURE
  PS::F64 pressure_kinetic[9] = {};//for calculation of pressure
  PS::F64 pressure_rf[9] = {};//for calculation of pressure
#endif
  void clear(){
    //    force    = 0.0;
    midforce = 0.0;
    torque   = 0.0;
    pot      = 0.0;
#ifdef CALC_PRESSURE
    //for calclation of pressure
    for(PS::S32 i = 0; i < 9; i++){
      pressure_kinetic[i] = 0.0;
      pressure_rf[i] = 0.0;
    }
#endif
  }
};

class FP{
public:
  PS::S32 id;
  PS::S32 type;
  PS::S32 cid[7];//added by kiyoshiro
  PS::F64 mass;
  PS::F64vec3 pos;
  PS::F64vec3 vel;
  PS::F64vec3 midvel;//修正速度ベルレ法で利用
  PS::F64vec3 force;
  PS::F64vec3 midforce;//修正速度ベルレ法で利用
  PS::S32 bond_num;
  PS::S32 bond_pair_id[3];//最終的に使いたい結合数が焼く3だから静的に3個確保しておく．サンプルではせいぜい2結合ほど．

  PS::F64vec3 torque;
#ifdef CALC_PRESSURE
  PS::F64 pressure_kinetic[9];//for calculation of pressure
  PS::F64 pressure_rf[9];//for calculation of pressure
#endif
#ifdef DISSIPATIVE_RANDOM
  PS::S32 seed;
#endif
  PS::F64 pot;
  PS::F64 search_radius;
  PS::F64 getRsearch() const {
    return this->search_radius;
  }
  PS::F64vec getPos() const { return pos; }
  void copyFromForce(const Force & f){
    //    force  = f.force;
    midforce = f.midforce;
    torque   = f.torque;
    pot      = f.pot;
#ifdef CALC_PRESSURE
    //for calclation of pressure
    for(PS::S32 i = 0; i < 9; i++){
      pressure_kinetic[i] = f.pressure_kinetic[i];
      pressure_rf[i] = f.pressure_rf[i];
    }
#endif
  }

  // writeXXX must be a const member function
  void writeAscii(FILE* fp) const {
    fprintf(fp, "%d %d %lf %lf %lf\n",
	    id, type, pos.x, pos.y, pos.z/*, angle.x, angle.y, angle.z, angle.w*/);
    
  }
  void readAscii(FILE* fp){
    fscanf(fp, "%d %d %lf %lf %lf\n",
	   &id, &type, &pos.x, &pos.y, &pos.z);
  }

  void writeBinary(FILE* fp) const {
    fwrite(this,sizeof(FP),1,fp);
  }
  void readBinary(FILE* fp){
    fread(this,sizeof(FP),1,fp);
  }


  void IntegrateBeforeForceCalc(const PS::F64 dt,const PS::F64 lh){
    //    const PS::F64 dth = dt * 0.5;
    const PS::F64 lambda = 0.65;//0.50;
    //    vel += force * dth;
    //    pos += vel * dt;
    //修正速度ベルレ法
    pos += dt * vel + 0.50 * dt *dt * force;
    midvel = vel + lambda * dt * force;

    // range must be [-lh,lh)

    //{
    while(pos.x >= lh) pos.x -= 2.0*lh;
    while(pos.y >= lh) pos.y -= 2.0*lh;
    while(pos.x < -lh) pos.x += 2.0*lh;
    while(pos.y < -lh) pos.y += 2.0*lh;
    //}

    //{
    while(pos.z >= lh) pos.z -= 2.0*lh;
    while(pos.z < -lh) pos.z += 2.0*lh;
    //}


    //    if(type != 1) angle = RichardsonMethod(angle,angvel,dth);
  }

  void IntegrateAfterForceCalc(const PS::F64 dt){
    //    const PS::F64 dth = dt * 0.5
    //vel += force * dth;
    vel += 0.50 * dt * (force + midforce);
    force = midforce;
#ifdef CALC_PRESSURE
    pressure_kinetic[0] = vel.x * vel.x;
    pressure_kinetic[1] = vel.x * vel.y;
    pressure_kinetic[2] = vel.x * vel.z;
    pressure_kinetic[3] = vel.y * vel.x;
    pressure_kinetic[4] = vel.y * vel.y;
    pressure_kinetic[5] = vel.y * vel.z;
    pressure_kinetic[6] = vel.z * vel.x;
    pressure_kinetic[7] = vel.z * vel.y;
    pressure_kinetic[8] = vel.z * vel.z;
#endif
    //    if(type != 1) angvel += torque * dth;
  }

  void CalcWallForce(){

  }
};

class EPI{
public:
  PS::S32 id;
  PS::S32 type;
  //  PS::S32 cid;//added by kiyoshiro
  PS::F64vec3 pos;
  //  Quaternion angle;
#ifdef DISSIPATIVE_RANDOM
  PS::F64vec3 vel;
  PS::F64vec3 midvel;//修正速度ベルレ法で利用
  PS::S32 seed;
#endif

  PS::F64vec3 getPos() const { return pos;}
  void copyFromFP(const FP & fp){
    pos = fp.pos;
#ifdef DISSIPATIVE_RANDOM
    vel = fp.vel;
    midvel = fp.midvel;
    seed = fp.seed;
#endif
    //    angle = fp.angle;
    type  = fp.type;
    id = fp.id;
  }
};

class EPJ{
public:
  PS::S32 id;
  PS::S32 type;
  PS::S32 cid;
  PS::F64vec3 pos;
  //  Quaternion angle;
#ifdef DISSIPATIVE_RANDOM
  PS::F64vec3 vel;
  PS::F64vec3 midvel;//修正速度ベルレ法で利用
  PS::S32 seed;
#endif
  PS::F64 search_radius;
  //  PS::S32 bond_num;//EPJは不必要？だと思われる．
  //  PS::S32 bond_pair_id[3];
  void copyFromFP(const FP & fp){ 
    pos = fp.pos;
#ifdef DISSIPATIVE_RANDOM
    vel = fp.vel;
    midvel = fp.midvel;
    seed = fp.seed;
#endif
    //    angle = fp.angle;
    type = fp.type;
    //cid = fp.cid;//added by kiyoshiro
    id = fp.id;
    search_radius = fp.search_radius;
  }
  PS::F64 getRSearch() const{
    return this->search_radius;
  }
  PS::F64vec getPos() const { return pos; }
  void setPos(const PS::F64vec & pos_new){
    pos = pos_new;
  }
  PS::S64 getId() const {return id;}
  //PS::F64 getCharge() const { return charge; }
};


struct CalcForceEpEp{
#ifdef DISSIPATIVE_RANDOM
  const PS::F64 dt;
  const PS::F64 temperature;
#endif

  CalcForceEpEp(
#ifdef DISSIPATIVE_RANDOM
		const PS::F64 _dt,
		const PS::F64 _temperature
#endif
		)
#ifdef DISSIPATIVE_RANDOM
    : dt(_dt), temperature(_temperature)
#endif
  {}
  void operator () (const EPI * ep_i,
		    const PS::S32 n_ip,
		    const EPJ * ep_j,
		    const PS::S32 n_jp,
		    Force * force){
    const PS::F64 ph = M_PI*0.5;



#ifdef DISSIPATIVE_RANDOM
    const PS::F64 sqrtdti = 1.0 / sqrt(dt);
    const PS::F64 sigma_dpd = sqrt(2.0 * temperature * gamma_dpd);
    //static XORShift rn;
    static TEA rn;
#endif
    const PS::F64 dij[2][2] = {{1.0,(1.0+dj)*0.5},{(1.0+dj)*0.5,dj}};
    for(int i=0;i<n_ip;i++){
      const PS::F64vec3 ri = ep_i[i].pos;
      //      const Quaternion  ai = ep_i[i].angle;
      const PS::S64 type_i = ep_i[i].type;
      PS::F64vec3 force_i  = 0.0;
      PS::F64vec3 force_temp  = 0.0;
      PS::F64vec3 torque_i = 0.0;
      PS::F64 pot_i = 0.0;
      PS::F64 sqrt3 = sqrt(3);
#ifdef CALC_PRESSURE
      PS::F64 pressure_rf_temp[9] = {};
#endif
      for(int j=0;j<n_jp;j++){
	const PS::F64vec3 rj = ep_j[j].pos;
	const PS::F64vec3 dr = ri - rj;
	const PS::F64 r2 = dr*dr;
	const PS::F64 d = 1.0;//dij[type_i][ep_j[j].type];
#ifdef DISSIPATIVE_RANDOM
	const PS::S32 seed_i = ep_i[i].seed;
	if(r2 > d*d || 1e-20f > r2) continue;
#else
	if(r2 > d*d || r2 == 0.0) continue;
#endif
	const PS::F64 rinv = 1.0 / sqrt(r2);
	const PS::F64 r2i = rinv*rinv;
	const PS::F64 r   = r2*rinv;
	// repulsive force
	force_i += dr * (coef_r[type_i][ep_j[j].type] * (d - r) * rinv);
	pot_i   += 0.5 * coef_r[type_i][ep_j[j].type] * (d - r) * (d - r);
#ifdef CALC_PRESSURE
	force_temp = dr * (coef_r[type_i][ep_j[j].type] * (d - r) * rinv);
#endif	
#if 1//#ifdef DISSIPATIVE_RANDOM
	// dissipative force
	const PS::F64vec dv = ep_i[i].midvel - ep_j[j].midvel;//midに変更
	const PS::F64 wij = 1.0 - r;
	const PS::F64 fd = gamma_dpd * wij*wij * (dv*dr) * rinv * rinv;
	force_i -= fd * dr;
	//#ifdef CALC_PRESSURE
        //force_temp -= fd * dr;
	//#endif
	// random force
	//const PS::F64 tmp = rn.drand((float)ep_i[i].vel.x,(float)ep_j[j].vel.x);//forNomalTEA
	//---------NEWRANDwithoutTEA
	/*	unsigned long long temp0 = teacode(ep_i[i].vel.x);
        unsigned long long temp1 = teacode(ep_j[j].vel.x);
        unsigned long long temp2 = temp0 ^ temp1;
	//std::cout<<temp2<<" "<<temp1<<" "<<temp0<<std::endl;
	PS::F64 newrand = (double)temp2;

	newrand = newrand / ULLONG_MAX;
	//       	std::cout<<newrand<<std::endl;
	const PS::F64 tmp = newrand * 2 * sqrt(3) - sqrt(3);*/
	//std::cout<<tmp<<std::endl;
	//---------NEWRANDwithoutTEA

	//---------NEWRANDwithoutTEA2                                                               
        unsigned long long temp0 = teacode(ep_i[i].vel.x+ep_j[j].vel.x);
        //unsigned long long temp1 = teacode(ep_j[j].vel.x);
        //unsigned long long temp2 = temp0 ^ temp1;
        //std::cout<<temp2<<" "<<temp1<<" "<<temp0<<std::endl;                                      
	PS::F64 newrand = (double)temp0;
        newrand = newrand / ULLONG_MAX;
        //              std::cout<<newrand<<std::endl;                                              
        const PS::F64 tmp = newrand * 2 * sqrt(3) - sqrt(3);
        //std::cout<<tmp<<std::endl;                                                                
        //---------NEWRANDwithoutTEA2

	//const PS::F64 tmp = (fabs(teacode(fabs(ep_i[i].vel.x - ep_j[j].vel.x)))/ULLONG_MAX) * 2 * sqrt3 - sqrt3;//for32TEA(1)
	//tmp = tmp * 2 * sqrt(3) - sqrt(3);
	const PS::F64 fr = sigma_dpd * wij * tmp * rinv * sqrtdti;
	force_i += fr * dr;
	//#ifdef CALC_PRESSURE
        //force_temp += fr * dr;
	//#endif
#endif
	if(type_i == 1 || ep_j[j].type == 1) continue;
	// attractive force and torque
	//	const Quaternion  aj = ep_j[j].angle;

#ifdef CALC_PRESSURE
	//for calculation of pressure
	pressure_rf_temp[0] += dr.x * force_temp.x;
	pressure_rf_temp[1] += dr.x * force_temp.y;
	pressure_rf_temp[2] += dr.x * force_temp.z;
	pressure_rf_temp[3] += dr.y * force_temp.x;
	pressure_rf_temp[4] += dr.y * force_temp.y;
	pressure_rf_temp[5] += dr.y * force_temp.z;
	pressure_rf_temp[6] += dr.z * force_temp.x;
	pressure_rf_temp[7] += dr.z * force_temp.y;
	pressure_rf_temp[8] += dr.z * force_temp.z;
#endif

      }
      force[i].pot    = pot_i;
      force[i].midforce  = force_i;//midに変更
      force[i].torque = torque_i;

#ifdef CALC_PRESSURE
      //for calculation of pressure                                                                
      force[i].pressure_rf[0] = pressure_rf_temp[0];
      force[i].pressure_rf[1] = pressure_rf_temp[1];
      force[i].pressure_rf[2] = pressure_rf_temp[2];
      force[i].pressure_rf[3] = pressure_rf_temp[3];
      force[i].pressure_rf[4] = pressure_rf_temp[4];
      force[i].pressure_rf[5] = pressure_rf_temp[5];
      force[i].pressure_rf[6] = pressure_rf_temp[6];
      force[i].pressure_rf[7] = pressure_rf_temp[7];
      force[i].pressure_rf[8] = pressure_rf_temp[8];
#endif

    }
  }
};

//by kiyoshiro
template<class Ttree>
void CalcIntraForce(FP &fp,Ttree &tree, const PS::F64vec box_size){
  EPJ ep;
  PS::F64 BOND_LENGTH = 0.86;
  PS::F64 BOND_COEF = 4.0;
  PS::F64vec box_size2i;

  if(fp.bond_num > 0){
    for(PS::S32 i=0; i < fp.bond_num; i++){
      EPJ* tmp = tree.getEpjFromId(fp.bond_pair_id[i]);
      if(tmp == NULL){
	//	std::cout<<"fp.bond_pair_id["<<i<<"] = "<<fp.bond_pair_id[i]<<std::endl;
	//	std::cout<<"fp.bond_num = "<<fp.bond_num<<std::endl;
	//std::cout<<"fp.id = "<<fp.id<<std::endl;
      }
      assert(tmp != NULL);
      ep = *tmp;
      while(fp.pos.x - ep.pos.x >= box_size.x) ep.pos.x += 2.0 * box_size.x;
      while(fp.pos.x - ep.pos.x < -box_size.x) ep.pos.x -= 2.0 * box_size.x;
      while(fp.pos.y - ep.pos.y >= box_size.y) ep.pos.y += 2.0 * box_size.y;
      while(fp.pos.y - ep.pos.y < -box_size.y) ep.pos.y -= 2.0 * box_size.y;
      while(fp.pos.z - ep.pos.z >= box_size.z) ep.pos.z += 2.0 * box_size.z;
      while(fp.pos.z - ep.pos.z < -box_size.z) ep.pos.z -= 2.0 * box_size.z;
      PS::F64vec drij = {fp.pos-ep.pos};
      PS::F64     rij = {sqrt(drij*drij)};
      const PS::F64 dij = rij - BOND_LENGTH;
      PS::F64vec midforce_i = (BOND_COEF * dij / rij) * drij;//midに変更
      fp.pot += 0.5*BOND_COEF*(dij*dij);
      //      fp.midforce -= (BOND_COEF * dij / rij) * drij;//midに変更
      fp.midforce -= midforce_i;//midに変更

#ifdef CALC_PRESSURE
      //for calculation of pressure
      fp.pressure_rf[0] -= drij.x * midforce_i.x;
      fp.pressure_rf[1] -= drij.x * midforce_i.y;
      fp.pressure_rf[2] -= drij.x * midforce_i.z;
      fp.pressure_rf[3] -= drij.y * midforce_i.x;
      fp.pressure_rf[4] -= drij.y * midforce_i.y;
      fp.pressure_rf[5] -= drij.y * midforce_i.z;
      fp.pressure_rf[6] -= drij.z * midforce_i.x;
      fp.pressure_rf[7] -= drij.z * midforce_i.y;
      fp.pressure_rf[8] -= drij.z * midforce_i.z;
#endif
    }
  }
#if 0
    EPJ* tmp = tree.getEpjFromId(1);
    assert(tmp != NULL);
    ep = *tmp;
    //
    while(fp.pos.x - ep.pos.x >= box_size.x) ep.pos.x += 2.0 * box_size.x;
    while(fp.pos.x - ep.pos.x < -box_size.x) ep.pos.x -= 2.0 * box_size.x;
    while(fp.pos.y - ep.pos.y >= box_size.y) ep.pos.y += 2.0 * box_size.y;
    while(fp.pos.y - ep.pos.y < -box_size.y) ep.pos.y -= 2.0 * box_size.y;
    while(fp.pos.z - ep.pos.z >= box_size.z) ep.pos.z += 2.0 * box_size.z;
    while(fp.pos.z - ep.pos.z < -box_size.z) ep.pos.z -= 2.0 * box_size.z;
    
    PS::F64vec drij =  {fp.pos-ep.pos};
    PS::F64     rij = {sqrt(drij*drij)};
    const PS::F64 dij = rij - BOND_LENGTH;
    fp.pot += 0.5*BOND_COEF*(dij*dij);
    fp.force -= (BOND_COEF * dij / rij) * drij;
  }
  else {
    EPJ* tmp = tree.getEpjFromId(0);
    assert(tmp != NULL);
    ep = *tmp;
    while(fp.pos.x - ep.pos.x >= box_size.x) ep.pos.x += 2.0 * box_size.x;
    while(fp.pos.x - ep.pos.x < -box_size.x) ep.pos.x -= 2.0 * box_size.x;
    while(fp.pos.y - ep.pos.y >= box_size.y) ep.pos.y += 2.0 * box_size.y;
    while(fp.pos.y - ep.pos.y < -box_size.y) ep.pos.y -= 2.0 * box_size.y;
    while(fp.pos.z - ep.pos.z >= box_size.z) ep.pos.z += 2.0 * box_size.z;
    while(fp.pos.z - ep.pos.z < -box_size.z) ep.pos.z -= 2.0 * box_size.z;

    PS::F64vec drij = {fp.pos-ep.pos};

    PS::F64     rij = {sqrt(drij*drij)};
    const PS::F64 dij = rij - BOND_LENGTH;
    fp.pot += 0.5*BOND_COEF*(dij*dij);
    fp.force -= (BOND_COEF * dij / rij) * drij;
  }
#elif 0
  EPJ ep[3];
  PS::S32 myid=0;
  for(PS::S32 i=0;i<3;i++){
    if(fp.cid[i] == fp.id){
      ep[i].copyFromFP(fp);
      myid = i;
      //std::cout<<"!!!"<<fp.id<<std::endl;
    }else{
      assert(fp.cid[i] >= 0);
      //std::cout<<fp.id<<std::endl;
      EPJ* tmp = tree.getEpjFromId(fp.cid[i]);
      assert(tmp != NULL);
      ep[i] = *tmp;
      // shift ep if ep is image particle in image cell
      while(fp.pos.x - ep[i].pos.x >= 0.5*box_size.x) ep[i].pos.x += box_size.x;
      while(fp.pos.x - ep[i].pos.x < -0.5*box_size.x) ep[i].pos.x -= box_size.x;
      while(fp.pos.y - ep[i].pos.y >= 0.5*box_size.y) ep[i].pos.y += box_size.y;
      while(fp.pos.y - ep[i].pos.y < -0.5*box_size.y) ep[i].pos.y -= box_size.y;
      while(fp.pos.z - ep[i].pos.z >= 0.5*box_size.z) ep[i].pos.z += box_size.z;
      while(fp.pos.z - ep[i].pos.z < -0.5*box_size.z) ep[i].pos.z -= box_size.z;
    }
  }
  PS::F64vec drij[2] = {ep[1].pos-ep[0].pos, ep[2].pos-ep[0].pos};
  PS::F64     rij[2] = {sqrt(drij[0]*drij[0]),sqrt(drij[1]*drij[1])};

  // bond
  PS::F64 BOND_LENGTH = 0.5;
  PS::F64 BOND_COEF = 4.0;

  if(myid == 0){ // case that fp is Oxygen
    const PS::F64 dij = rij[0] - BOND_LENGTH;
    const PS::F64 dkj = rij[1] - BOND_LENGTH;
    fp.pot += 0.5*BOND_COEF*(dij*dij + dkj*dkj);
    fp.force += (2.0*BOND_COEF * dij / rij[0]) * drij[0];
    fp.force += (2.0*BOND_COEF * dkj / rij[1]) * drij[1];
  }else{
    int bond_id = myid - 1;
    const PS::F64 dij = rij[bond_id] - BOND_LENGTH;
    fp.pot += 0.5*BOND_COEF*dij*dij;
    fp.force -= (2.0*BOND_COEF*dij/rij[bond_id]) * drij[bond_id];
  }
#endif
}

//template<class Ttree>
//void CalcIntraForce(FP &fp,Ttree &tree, const PS::F64vec box_size){
/*
template<class Ttree>
void CalcPressure(FP &fp, Ttree &tree, const PS::F64 volume){
  PS::F64 kinetic_loc[9];
  PS::F64 rf_loc[9];

  //for



  }*/
//end by kiyoshiro

void MakeFaceCubicCenter(const PS::S32 n_tot,
			 PS::F64 *&mass,
			 PS::F64vec *&pos,
			 PS::F64vec *&vel,
			 //			 Quaternion *&angle,
			 // PS::F64vec *&angvel,
			 const double density,
			 const int seed = 0){
  //static const double PI = atan(1.0) * 4.0;
  PS::MTTS mt;
  double cell_size = pow((double)n_tot/density,1./3.);
  int nunit = 1;
#ifdef UNI_SET
  while(4*nunit*nunit*nunit < n_tot) nunit++;
#endif
  //if (n_tot != 4*nunit*nunit*nunit){
  //std::cerr << "MakeFaceCubicCenter: n_tot and 4*nunit^3 must be the same. "
  //    << n_tot << "!= " << 4*nunit*nunit*nunit <<std::endl;
    // PS::Abort();
  //}
  mt.init_genrand(PS::Comm::getRank()*PS::Comm::getNumberOfThread()+PS::Comm::getThreadNum());


#ifdef UNI_SET
  double unit_size = cell_size/(double)nunit;
  double ush = unit_size * 0.5;
  PS::F64vec unit[4];
  unit[0].x = 0.0; unit[1].x = ush; unit[2].x = 0.0; unit[3].x = ush;
  unit[0].y = 0.0; unit[1].y = ush; unit[2].y = ush; unit[3].y = 0.0;
  unit[0].z = 0.0; unit[1].z = 0.0; unit[2].z = ush; unit[3].z = ush;

  int ip=0;
  for(int i=0; i<nunit; i++){
    for(int j=0; j<nunit; j++){
      for(int k=0; k<nunit; k++){
        for(int l=0; l<4; l++){
	  //      for(int l=0; l<1; l++){                                                              
          /*      pos[ip].x = i*unit_size + unit[l].x + 0.1*ush;                                     
          pos[ip].y = j*unit_size + unit[l].y + 0.1*ush;                                             
          pos[ip].z = k*unit_size + unit[l].z + 0.1*ush;*/
          pos[ip].x = i*unit_size + 0.1*unit[l].x + 0.1*ush;
          pos[ip].y = j*unit_size + 0.1*unit[l].y + 0.1*ush;
          pos[ip].z = k*unit_size + 0.1*unit[l].z + 0.1*ush;

          ip++;
        }
      }
    }
  }
#endif
  std::cout<<ip<<std::endl;










#ifdef RANDOM_SET  
  //  PS::F64 random = 0;
  PS::F64 random_x = 0;
  PS::F64 random_y = 0;
  PS::F64 random_z = 0;


  int ip=0;
  for(int i=0; i*12<n_tot;i++){
    random_x = (mt.genrand_res53() - 0.50) * cell_size;
    random_y = (mt.genrand_res53() - 0.50) * cell_size;
    random_z = (mt.genrand_res53() - 0.50) * cell_size;
    for(int l=0; l<12; l++){
      pos[ip].x = random_x+l*0.80;
      pos[ip].y = random_y;
      pos[ip].z = random_z;
      ip++;
    }
  }
#endif



  assert(ip == n_tot);

  for(int i=0; i<n_tot; i++){
    mass[i] = 1.0;
    const double v_max = 0.1;
    do {
      vel[i][0] = (2. * mt.genrand_res53() - 1.) * v_max;
      vel[i][1] = (2. * mt.genrand_res53() - 1.) * v_max;
      vel[i][2] = (2. * mt.genrand_res53() - 1.) * v_max;
    }while(vel[i] * vel[i] >= v_max * v_max);
  }

  for(int i=0; i<n_tot; i++){
  }


  for(int i=0; i<n_tot; i++){
    pos[i].x -= 0.5*cell_size;
    pos[i].y -= 0.5*cell_size;
    pos[i].z -= 0.5*cell_size;

    //{
    while(pos[i].x >= 0.5*cell_size) pos[i].x -= cell_size;
    while(pos[i].y >= 0.5*cell_size) pos[i].y -= cell_size;
    while(pos[i].z >= 0.5*cell_size) pos[i].z -= cell_size;
    while(pos[i].x < -0.5*cell_size) pos[i].x += cell_size;
    while(pos[i].y < -0.5*cell_size) pos[i].y += cell_size;
    while(pos[i].z < -0.5*cell_size) pos[i].z += cell_size;
    //}


  }

  PS::F64vec cm_vel = 0.0;
  double  cm_mass = 0.0;
  for(int i=0; i<n_tot; i++){
    cm_vel += mass[i] * vel[i];
    cm_mass += mass[i];
  }
  cm_vel /= cm_mass;
  for(int i=0; i<n_tot; i++){
    vel[i] -= cm_vel;
  }
}



template<class Tpsys>
void SetParticles(Tpsys & psys,
		  const PS::S32 n_tot,
		  const double density,
		  const double temperature,
		  const double solvent_ratio){
  /*#if 0 // particles are generated on each rank (need bug fix)
  PS::F64    *mass   = new PS::F64[n_tot];
  PS::F64vec *pos    = new PS::F64vec[n_tot];
  PS::F64vec *vel    = new PS::F64vec[n_tot];


  MakeFaceCubicCenter(n_tot, mass, pos, vel, density);

  PS::S32 n_proc = PS::Comm::getNumberOfProc();
  PS::S32 rank = PS::Comm::getRank();

  PS::S32 n_loc = n_tot / n_proc;
  PS::S32 i_h = n_loc*rank;
  if(n_loc%n_proc > rank) i_h += rank;
  if(n_loc%n_proc > rank) n_loc++;

  psys.setNumberOfParticleLocal(n_loc);
  for(int i=0; i<n_loc; i++){
    const int id = i + i_h;
    //assert(id < n_tot);
    psys[i].mass   = mass[id];
    psys[i].pos    = pos[id];
    psys[i].vel    = vel[id];

    psys[i].id     = id;
#if 1 //by kiyoshiro 
    if(id % 7 == 0){
      psys[i].type = 0;
    }
    else{
      psys[i].type = 1;
    }
#else //end by kiyoshiro
    psys[i].type   = 0;
#endif
    psys[i].search_radius = 3.0;
  }
  const PS::S32 n_solvent = (PS::S32)(n_loc * solvent_ratio);
  int n = 0;
  while(n<n_solvent){
    const PS::S32 pivot = rand()%n_loc;
    if(psys[pivot].type == 0){
      psys[pivot].type = 1;

      n++;

      //std::cout << n << std::endl;
    }
  }



  if(mass   != nullptr) delete [] mass;
  if(pos    != nullptr) delete [] pos;
  if(vel    != nullptr) delete [] vel;

  #else*/ // all particles are generated on rank 0
  if(PS::Comm::getRank() == 0){
    PS::F64    *mass   = new PS::F64[n_tot];
    PS::F64vec *pos    = new PS::F64vec[n_tot];
    PS::F64vec *vel    = new PS::F64vec[n_tot];


    MakeFaceCubicCenter(n_tot, mass, pos, vel,/* angle, angvel,*/ density);


    psys.setNumberOfParticleLocal(n_tot);
    for(int i=0; i<n_tot; i++){
      psys[i].mass   = mass[i];
      psys[i].pos    = pos[i];
      psys[i].vel    = vel[i];

      psys[i].id     = i;
      if(i==165811)std::cout<<"psys["<<i<<"].id = "<<i<<std::endl;
      //by kiyoshiro
      for(int c=0;c<3;c++){
	const int offset = (i / 3)*3;
	psys[i].cid[c] = offset + c;
      }
      //end by kiyoshiro
      //by kiyoshiro "for bond_num and bond_pair_id"
      //1 bond
      if(psys[i].id % 12 == 0){
	psys[i].bond_num = 1;
	psys[i].bond_pair_id[0] = psys[i].id + 1; 
      }
      
      else if(psys[i].id % 12 == 6 ||
	      psys[i].id % 12 == 11){
	psys[i].bond_num = 1;
	psys[i].bond_pair_id[0] = psys[i].id - 1; 
      }
      //2 bond
      else if(psys[i].id % 12 == 1 ||
	      psys[i].id % 12 == 3 ||
	      psys[i].id % 12 == 4 ||
	      psys[i].id % 12 == 5 ||
	      psys[i].id % 12 == 8 ||
	      psys[i].id % 12 == 9 ||
	      psys[i].id % 12 == 10){
	psys[i].bond_num = 2;
	psys[i].bond_pair_id[0] = psys[i].id + 1; 
	psys[i].bond_pair_id[1] = psys[i].id - 1; 
      }
      else if(psys[i].id % 12 == 7){
	psys[i].bond_num = 2;
        psys[i].bond_pair_id[0] = psys[i].id + 1;
        psys[i].bond_pair_id[1] = psys[i].id - 5;
      }
      //3 bond
      else if(psys[i].id % 12 == 2){
	psys[i].bond_num = 3;
	psys[i].bond_pair_id[0] = psys[i].id + 1; 
	psys[i].bond_pair_id[1] = psys[i].id - 1; 
	psys[i].bond_pair_id[2] = psys[i].id + 5; 
      }
      else{
	std::cout<<"error"<<std::endl;
      }

      //end by kiyoshiro
      if(psys[i].id % 12 == 0)psys[i].type = 5;
      else if(psys[i].id % 12 == 1)psys[i].type = 4;
      else if(psys[i].id % 12 == 2)psys[i].type = 3;
      else if(psys[i].id % 12 == 3)psys[i].type = 1;
      else if(psys[i].id % 12 == 4)psys[i].type = 2;
      else if(psys[i].id % 12 == 5)psys[i].type = 0;
      else if(psys[i].id % 12 == 6)psys[i].type = 1;
      else if(psys[i].id % 12 == 7)psys[i].type = 3;
      else if(psys[i].id % 12 == 8)psys[i].type = 0;
      else if(psys[i].id % 12 == 9)psys[i].type = 0;
      else if(psys[i].id % 12 == 10)psys[i].type = 0;
      else if(psys[i].id % 12 == 11)psys[i].type = 1;

      //psys[i].id >= 1200 is water17040,19392
#ifdef RANDOM_SET
      if(psys[i].id >= 19392)psys[i].bond_num = 0;
      if(psys[i].id >= 19392)psys[i].type = 6;
#endif
#ifdef UNI_SET
      if(psys[i].id >= 18360)psys[i].bond_num = 0;//18360POPC17% 18900POPC17.5% 19116POPC17.7% 19440POPC18% 21600POPC20% 32400POPC30% 43200POPC40%
      if(psys[i].id >= 18360)psys[i].type = 6;
#endif
      //psys[i].type   = 0;
      psys[i].search_radius = 6.0;
    }
    
    const PS::S32 n_solvent = (PS::S32)(n_tot * solvent_ratio);
    int n = 0;
    while(n<n_solvent){
      const PS::S32 pivot = rand()%n_tot;
      if(psys[pivot].type == 0){
	psys[pivot].type = 1;
	n++;

	//std::cout << n << std::endl;
      }
    }

    if(mass   != nullptr) delete [] mass;
    if(pos    != nullptr) delete [] pos;
    if(vel    != nullptr) delete [] vel;
    //psys.setNumberOfParticleLocal(n_tot);//add temporary kiyoshiro
    /*  }else{
	psys.setNumberOfParticleLocal(0); is it necesary?*/
  }

  //#endif
  ScaleVelocity(psys,temperature);
}

template<class Tpsys>
void RemoveTotalMomentum(Tpsys &system){
  const PS::S32 n_loc = system.getNumberOfParticleLocal();
  PS::F64vec cm_vel_loc = 0.0;
  PS::F64  cm_mass_loc = 0.0;
  for(int i=0; i<n_loc; i++){
    cm_vel_loc += system[i].mass * system[i].vel;
    cm_mass_loc += system[i].mass;
  }
  PS::F64vec cm_vel=0.0;
  PS::F64 cm_mass=0.0;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
  MPI::COMM_WORLD.Allreduce(&cm_vel_loc.x, &cm_vel.x, 1, PS::GetDataType<PS::F64>(), MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&cm_vel_loc.y, &cm_vel.y, 1, PS::GetDataType<PS::F64>(), MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&cm_vel_loc.z, &cm_vel.z, 1, PS::GetDataType<PS::F64>(), MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&cm_mass_loc, &cm_mass, 1, PS::GetDataType<PS::F64>(), MPI::SUM);
#else
  cm_vel = cm_vel_loc;
  cm_mass = cm_mass_loc;
#endif

  cm_vel /= cm_mass;
  for(int i=0; i<n_loc; i++){
    system[i].vel -= cm_vel;
  }
}

template<class Tpsys>
void ScaleVelocity(Tpsys & system,const PS::F64 T){
  const PS::S32 natom_local = system.getNumberOfParticleLocal();
  PS::F64 etra_loc = 0.0;
  PS::F64vec erot_loc = 0.0;
  for(PS::S32 i=0; i<natom_local; i++){
    etra_loc += system[i].mass * system[i].vel * system[i].vel;
  }
  etra_loc *= 0.5;
  erot_loc *= 0.5;
  PS::S32 natom;
  PS::F64 etra;
  PS::F64vec erot;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
  MPI::COMM_WORLD.Allreduce(&etra_loc, &etra, 1, PS::GetDataType<PS::F64>(), MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&erot_loc.x, &erot.x, 1, PS::GetDataType<PS::F64>(), MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&erot_loc.y, &erot.y, 1, PS::GetDataType<PS::F64>(), MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&erot_loc.z, &erot.z, 1, PS::GetDataType<PS::F64>(), MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&natom_local, &natom, 1, PS::GetDataType<PS::S32>(), MPI::SUM);
#else
  etra = etra_loc;
  erot = erot_loc;
  natom = natom_local;
#endif
    const PS::F64 s_tra = sqrt(1.5*natom*T / etra);
  //  const PS::F64 s_tra = sqrt(1.5*natom*T / etra);
      for(PS::S32 i=0;i<natom_local;i++) system[i].vel *= s_tra;
#if 1
  const PS::F64 s_rotx = sqrt(0.5*natom*T / erot.x);
  const PS::F64 s_roty = sqrt(0.5*natom*T / erot.y);
  const PS::F64 s_rotz = sqrt(0.5*natom*T / erot.z);
  for(PS::S32 i=0;i<natom_local;i++){
  }
#endif
  RemoveTotalMomentum(system);
}


//for calc momentum=0 
template<class Tpsys>
void checkmomentum0(const Tpsys & system){
  //for calc momentum=0                                                                              
  PS::F64 vel_sum_x_loc = 0;
  PS::F64 vel_sum_y_loc = 0;
  PS::F64 vel_sum_z_loc = 0;
  PS::F64 vel_sum_x = 0;
  PS::F64 vel_sum_y = 0;
  PS::F64 vel_sum_z = 0;

  const PS::S32 nbody = system.getNumberOfParticleLocal();
  for(PS::S32 i=0; i<nbody; i++){
    //for calc momentum=0     

    vel_sum_x_loc += system[i].vel.x;
    vel_sum_y_loc += system[i].vel.y;
    vel_sum_z_loc += system[i].vel.z;
  }
  MPI::COMM_WORLD.Allreduce(&vel_sum_x_loc, &vel_sum_x, 1, PS::GetDataType<PS::S32>(), MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&vel_sum_y_loc, &vel_sum_y, 1, PS::GetDataType<PS::S32>(), MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&vel_sum_z_loc, &vel_sum_z, 1, PS::GetDataType<PS::S32>(), MPI::SUM);
  std::cout << vel_sum_x << " " << vel_sum_y << " " << vel_sum_z << std::endl;
  
  
}






template<class Tpsys>
void CalcKineticEnergy(const Tpsys & system,
		       PS::F64 & etra,
		       PS::F64 & erot,
		       PS::S32 & deg_free){
  PS::F64 etra_loc = 0.0;
  PS::F64 erot_loc = 0.0;
  PS::S32 deg_free_loc = 0;
  //for calc momentum=0
  /*PS::F64 vel_sum_x_loc = 0;
  PS::F64 vel_sum_y_loc = 0;
  PS::F64 vel_sum_z_loc = 0;
  PS::F64 vel_sum_x = 0;
  PS::F64 vel_sum_y = 0;
  PS::F64 vel_sum_z = 0;*/


  const PS::S32 nbody = system.getNumberOfParticleLocal();
  for(PS::S32 i=0; i<nbody; i++){
    etra_loc += 0.5 * system[i].mass * system[i].vel * system[i].vel;

    //for calc momentum=0
    /*    vel_sum_x_loc += system[i].vel.x;
    vel_sum_y_loc += system[i].vel.y;
    vel_sum_z_loc += system[i].vel.z;*/
    // intertia is unit {1.0,0.0,0.0; 0.0,1.0,0.0; 0.0,0.0,1.0}
    // A(q) I^-1_b A^T(q) = I (unit)
    /*
    const PS::F64mat3asym A = Rotate(system[i].angle) * Trans(Rotate(system[i].angle));
    printf(" %lf %lf %lf\n %lf %lf %lf\n %lf %lf %lf\n\n",
	   A.xx,A.xy,A.xz,A.yx,A.yy,A.yz,A.zx,A.zy,A.zz);
    */
    if(system[i].type == 0){

      deg_free_loc += 6;
    }else{
      deg_free_loc += 3;
    }
  }
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
  MPI::COMM_WORLD.Allreduce(&etra_loc, &etra, 1, PS::GetDataType<PS::F64>(), MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&erot_loc, &erot, 1, PS::GetDataType<PS::F64>(), MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&deg_free_loc, &deg_free, 1, PS::GetDataType<PS::S32>(), MPI::SUM);
  //for calc momentum=0
  /*MPI::COMM_WORLD.Allreduce(&vel_sum_x_loc, &vel_sum_x, 1, PS::GetDataType<PS::S32>(), MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&vel_sum_y_loc, &vel_sum_y, 1, PS::GetDataType<PS::S32>(), MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&vel_sum_z_loc, &vel_sum_z, 1, PS::GetDataType<PS::S32>(), MPI::SUM);
  std::cout << vel_sum_x << " " << vel_sum_y << " " << vel_sum_z << std::endl;*/
#else
  etra = etra_loc;
  erot = erot_loc;
  deg_free = deg_free_loc;
#endif
}

template<class Tpsys>
void CalcMomentum(const Tpsys & system,
		  PS::F64vec & mom){
  PS::F64vec mom_loc = 0.0;
  const PS::S32 nbody = system.getNumberOfParticleLocal();
  for(PS::S32 i=0; i<nbody; i++){
    mom_loc += 0.5 * system[i].mass * system[i].vel;
  }
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
  MPI::COMM_WORLD.Allreduce(&mom_loc.x, &mom.x, 1, PS::GetDataType<PS::F64>(), MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&mom_loc.y, &mom.y, 1, PS::GetDataType<PS::F64>(), MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&mom_loc.z, &mom.z, 1, PS::GetDataType<PS::F64>(), MPI::SUM);
#else
  mom = mom_loc;
#endif
}

template<class Tpsys>
void CalcEnergy(const Tpsys & system,
                PS::F64 & etot,
                PS::F64 & etra,
		PS::F64 & erot,//use as temperature
                PS::F64 & epot,
                const bool clear=true){
  if(clear){
    etot = etra = erot = epot = 0.0;
  }
  PS::F64 etot_loc = 0.0;
  PS::F64 etra_loc = 0.0;
  PS::F64 erot_loc = 0.0;
  PS::F64 epot_loc = 0.0;
  PS::S32 numprocs;
  const PS::S32 nbody = system.getNumberOfParticleLocal();
  for(PS::S32 i=0; i<nbody; i++){
    etra_loc += 0.5 * system[i].mass * system[i].vel * system[i].vel;
    // intertia is unit {1.0,0.0,0.0; 0.0,1.0,0.0; 0.0,0.0,1.0}
    // A(q) I^-1_b A^T(q) = I (unit)
    /*
    const PS::F64mat3asym A = Rotate(system[i].angle) * Trans(Rotate(system[i].angle));
    printf(" %lf %lf %lf\n %lf %lf %lf\n %lf %lf %lf\n\n",
	   A.xx,A.xy,A.xz,A.yx,A.yy,A.yz,A.zx,A.zy,A.zz);
    */

    epot_loc += system[i].pot;
  }
  epot_loc *= 0.5;
  etot_loc = etra_loc /*+ erot_loc*/ + epot_loc;
  erot_loc = etra_loc / nbody / 3 * 2;//temperature
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
  MPI::COMM_WORLD.Allreduce(&etot_loc, &etot, 1, PS::GetDataType<PS::F64>(), MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&epot_loc, &epot, 1, PS::GetDataType<PS::F64>(), MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&etra_loc, &etra, 1, PS::GetDataType<PS::F64>(), MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&erot_loc, &erot, 1, PS::GetDataType<PS::F64>(), MPI::SUM);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  erot/=(PS::F64)numprocs;
#else
  etot = etot_loc;
  epot = epot_loc;
  etra = etra_loc;
  erot = erot_loc;
#endif
}

int main(int argc, char *argv[]){
  PS::Initialize(argc, argv);


  std::cout<<std::setprecision(15);
  std::cerr<<std::setprecision(15);
  PS::F64 theta = 0.5;
  const PS::S32 n_leaf_limit = 8;
#ifdef RANDOM_SET
  long long int n_tot = 107808;//6912;//1492992;//107808;//97556;//2916;//186624;//108000;//23328;//108000;//62500;//2916;//1372;//256*3;864
#endif
#ifdef UNI_SET
  long long int n_tot = 108000;
#endif
  PS::F64 density     = 3.0;
  PS::F64 temperature = 1.0;//0.8;
  PS::F64 solvent_ratio = 0.0;

  PS::F64 dt = 0.050;//0.0001;

  PS::S32 n_group_limit = 64;
  PS::S32 nstep      = 50000;//300000;
  PS::S32 nstep_eq   = 1000;//10000;
  PS::S32 nstep_snp  = 100;
  PS::S32 nstep_diag = 100;
  //  std::cout<<PS::Comm::getRank()<<std::endl;
  //add by kiyoshiro{
  //  PS::S32 fprank = PS::Comm::getRank();
  //if(fprank == 0){
  std::ofstream fpo1;//("placs.xyz");
  FILE *fpo2;
  FILE *fpo3;//(pressure.dat);
  //}
  std::string input_file = "";

  //char sinput[1024];
  char dir_name[1024];
  int c;
  sprintf(dir_name,"./result");
  while((c=getopt(argc,argv,"o:N:d:T:s:e:S:D:t:c:r:n:i:h")) != -1){
    switch(c){
    case 'o':
      sprintf(dir_name,"%s",optarg);
      break;
    case 'N':
      n_tot = atoi(optarg);
      std::cerr<<"n_tot="<<n_tot<<std::endl;
      break;
    case 'd':
      density = atof(optarg);
      std::cerr<<"density="<<density<<std::endl;
      break;
    case 'T':
      temperature = atof(optarg);
      std::cerr<<"temperature="<<temperature<<std::endl;
      break;
    case 's':
      nstep = atoi(optarg);
      std::cerr<<"nstep="<<nstep<<std::endl;
      break;
    case 'e':
      nstep_eq = atoi(optarg);
      std::cerr<<"nstep_eq="<<nstep_eq<<std::endl;
      break;
    case 'S':
      nstep_snp = atoi(optarg);
      std::cerr<<"nstep_snp="<<nstep_snp<<std::endl;
      break;
    case 'D':
      nstep_diag = atoi(optarg);
      std::cerr<<"nstep_diag="<<nstep_diag<<std::endl;
      break;
    case 't':
      dt = atof(optarg);
      std::cerr<<"dt="<<dt<<std::endl;
      break;
    case 'r':
      solvent_ratio = atof(optarg);
      std::cerr<<"solvent_ratio="<<solvent_ratio<<std::endl;
      break;
    case 'n':
      n_group_limit = atoi(optarg);
      std::cerr<<"n_group_limit="<<n_group_limit<<std::endl;
      break;
    case 'i':
      input_file = optarg;
      std::cerr<<"input_file="<<input_file<<std::endl;
      break;
    case 'h':
      std::cerr<<"N: n_tot (default: 1000)"<<std::endl;
      std::cerr<<"d: number density (default: 1.05)"<<std::endl;
      std::cerr<<"T: temperature (default: 0.8)"<<std::endl;
      std::cerr<<"s: number of steps (default: 1000)"<<std::endl;
      std::cerr<<"S: time step for snapshot(default: 100)"<<std::endl;
      std::cerr<<"D: time step for diag(default: 100)"<<std::endl;
      std::cerr<<"e: number of steps for equilibration(default: 1000)"<<std::endl;
      std::cerr<<"r: ratio of sovent (default: 0.0)"<<std::endl;
      std::cerr<<"o: dir name of output (default: ./result)"<<std::endl;
      std::cerr<<"i: checkpoint file name"<<std::endl;
      std::cerr<<"t: time step (default: 0.0001)"<<std::endl;
      std::cerr<<"n: n_group_limit (default: 64.0)"<<std::endl;
      return 0;
    }
  }

  PS::F64 boxdh = 0.5*powf((double)n_tot/density,1./3.);
  PS::F64 volume = 8 * boxdh * boxdh *boxdh;
  if(PS::Comm::getRank()==0){
    fprintf(stderr, "boxdh = %lf\n",boxdh);
    //add by kiyoshiro
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL) && defined(PARTICLE_SIMULATOR_MPI_PARALLEL)   
    PS::S32 numprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    std::cout<<"process_num="<<numprocs<<std::endl;
    std::cout << "OpenMP : On, threads =" << omp_get_max_threads() << std::endl;
#endif
    //end by kiyoshiro  
  }
  struct stat st;
  if(stat(dir_name, &st) != 0) {
    PS::S32 rank = PS::Comm::getRank();
    PS::S32 ret_loc, ret=0;
    if(rank == 0)
      ret_loc = mkdir(dir_name, 0777);
    PS::Comm::broadcast(&ret_loc, ret);
    if(ret == 0) {
      if(rank == 0)
	fprintf(stderr, "Directory \"%s\" is successfully made.\n", dir_name);
    } else {
      fprintf(stderr, "Directory %s fails to be made.\n", dir_name);
      PS::Abort();
      exit(0);
    }
  }

  std::ofstream fout_eng;
  std::ofstream fout_tcal;
  if(PS::Comm::getRank() == 0){
    char sout_de[1024];
    char sout_tcal[1024];
    sprintf(sout_de, "%s/t-de.dat", dir_name);
    sprintf(sout_tcal, "%s/t-tcal.dat", dir_name);
    std::cerr<<sout_de<<std::endl;
    std::cerr<<sout_tcal<<std::endl;
    fout_eng.open(sout_de);
    fout_tcal.open(sout_tcal);
  }

  PS::ParticleSystem<FP> system_janus;
  system_janus.initialize();

  PS::S32 n_grav_glb = n_tot;
  if(input_file == ""){
    SetParticles(system_janus, n_tot, density, temperature,solvent_ratio);
    //if(PS::Comm::getRank()==0) fprintf(stderr,"Particles are generated!\n");
  }else{
    BinaryHeader header(n_tot);
    system_janus.readParticleBinary(input_file.c_str(),header);
    //    if(PS::Comm::getRank()==0) fprintf(stderr,"Particles are read from %s!\n",input_file.c_str());
  }

  const PS::F64 coef_ema = 0.3;
  PS::DomainInfo dinfo;
  dinfo.initialize(coef_ema);


  dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
  dinfo.setPosRootDomain(PS::F64vec(-boxdh,-boxdh,-boxdh),
			 PS::F64vec( boxdh, boxdh, boxdh));

  dinfo.collectSampleParticle(system_janus);
  dinfo.decomposeDomain();
  system_janus.exchangeParticle(dinfo);
  PS::TreeForForceShort<Force, EPI, EPJ>::Scatter tree_janus;
  tree_janus.initialize(n_grav_glb, theta, n_leaf_limit, n_group_limit);
  /*#ifdef DISSIPATIVE_RANDOM
  for(int i=0;i<system_janus.getNumberOfParticleLocal();i++) system_janus[i].seed = rand();
  #endif*/
  tree_janus.calcForceAllAndWriteBack(CalcForceEpEp(
#ifdef DISSIPATIVE_RANDOM
						    dt,temperature
#endif
						    ),  system_janus, dinfo);

#if 1
  for(int i=0;i<system_janus.getNumberOfParticleLocal();i++){
    CalcIntraForce(system_janus[i],tree_janus,boxdh);
  }
#endif

  PS::F64 Epot0, Etra0, Erot0, Etot0, Epot1, Etra1, Erot1, Etot1;
  //ScaleVelocity(system_janus,temperature);
  CalcEnergy(system_janus, Etot0, Etra0, Erot0, Epot0);
  //if(PS::Comm::getRank() == 0) printf("Etot = %lf, Epot = %lf, Etra = %lf, Erot = %lf\n",Etot0,Epot0,Etra0,Erot0);

  PS::S32 snp_id = 0;
  PS::S32 time_snp = 0;
  PS::S32 time_diag = 0;
  bool isInitialized = false;
  PS::F64 time_sys = -dt * nstep_eq;

  PS::F64 Epot_ave = 0.0, Ekin_ave = 0.0;
  int n_loc = system_janus.getNumberOfParticleLocal();
  int n_glo = system_janus.getNumberOfParticleGlobal();
  PS::F64 time_integ = 0.0;
  PS::F64 time_intra = 0.0;//by kiyoshiro for intra time

  PS::S32 numprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  
  PS::S32 rank = PS::Comm::getRank();
  //#ifdef CALC_PRESSURE
  //PS::F64 pressure_kinetic[9];
  //for(PS::S32 k = 0; k < 9; k++){
  //  pressure_kinetic[k] = 0.0;
  //}
  //#endif

  //////////////////////////main loop///////////////////////////////
  for(int s=-nstep_eq;s<nstep;s++){
    //add by kiyoshiro{
    if(s==-9900){//search_radiusを変更する
      for(PS::S32 kk=0;kk<n_loc;kk++){
	system_janus[kk].search_radius=4.5;
      }
    }
    //output file start
    if(s % 1000 == 0){
      //        if(s == 1){

      if(rank == 0){ 
	fpo2 = fopen( "nextplace.gro" , "a" );
	if( fpo2 == NULL ) {
	  printf( "ファイルオープンエラー\n" );
	  return -1;
	}
		
	fprintf(fpo2,"%d\n%d\n",s,n_glo);
	fclose(fpo2);
      }
    
      for(PS::S32 r = 0; r < numprocs; r++){
	PS::Comm::barrier();
	usleep(3 * 10000);
	if(rank == r){
      
	  fpo2 = fopen( "nextplace.gro" , "a" );
	  if( fpo2 == NULL ) {
	    printf( "ファpイルオープンエラー\n" );
	    return -1;
	  }
	  for(PS::S32 i=0;i<n_loc;i++){
	    if(system_janus[i].type == 6){
	      fprintf(fpo2,"%9d%-5s%5s%9d%8.3f%8.3f%8.3f\n",	\
		      system_janus[i].id,			\
		      "H","H",				\
		      system_janus[i].id,			\
		      system_janus[i].pos.x+boxdh,			\
		      system_janus[i].pos.y+boxdh,			\
		      system_janus[i].pos.z+boxdh);
	    }else if(system_janus[i].type == 0 || system_janus[i].type == 1 || system_janus[i].type == 2){
	      fprintf(fpo2,"%9d%-5s%5s%9d%8.3f%8.3f%8.3f\n",	\
		      system_janus[i].id,			\
		      "C","C",				\
		      system_janus[i].id,			\
		      system_janus[i].pos.x+boxdh,			\
		      system_janus[i].pos.y+boxdh,			\
		      system_janus[i].pos.z+boxdh);
	    }else{
              fprintf(fpo2,"%9d%-5s%5s%9d%8.3f%8.3f%8.3f\n",    \
                      system_janus[i].id,                       \
                      "N","N",                              \
                      system_janus[i].id,                       \
                      system_janus[i].pos.x+boxdh,                    \
                      system_janus[i].pos.y+boxdh,                    \
                      system_janus[i].pos.z+boxdh);
            }
	  }
	  fclose(fpo2);
	}
      }
      
      PS::Comm::barrier();
      if(PS::Comm::getRank() == 0){
	fpo2 = fopen( "nextplace.gro" , "a" );
	if( fpo2 == NULL ) {
	  printf( "ファイルオープンエラー\n" );
	  return -1;
	}
	fprintf(fpo2,"%8.3f%8.3f%8.3f\n",boxdh*2,boxdh*2,boxdh*2);
	fclose(fpo2);
      }
    }
    
    //output file end 
    MPI_Barrier( MPI_COMM_WORLD ); 
#if 1
    if(s < 0){
      ScaleVelocity(system_janus,temperature);
    }
#endif
    //    if(s%1000==0) RemoveTotalMomentum(system_janus);
    PS::Timer timer;
    timer.reset();
    timer.start();
    if(s == time_snp){
      CDVHeader cdvh(PS::F64vec(-boxdh,-boxdh,-boxdh),PS::F64vec(boxdh,boxdh,boxdh));
      char filename[256];
      sprintf(filename, "%s/%05d.cdv", dir_name, snp_id);
      system_janus.writeParticleAscii(filename, cdvh);

      BinaryHeader bh(n_tot);
      sprintf(filename,"%s/checkpoint", dir_name);

      system_janus.writeParticleBinary(filename,bh);

      snp_id++;
      time_snp += nstep_snp;
    }



    //if(!isInitialized && s == 0){
      CalcEnergy(system_janus, Etot0, Etra0, Erot0, Epot0);
      if(PS::Comm::getRank() == 0) fprintf(stderr,"step = %d, Etot0 = %lf, Epot0 = %lf, Etra0 = %lf, temperature = %lf\n",s,Etot0,Epot0,Etra0,Erot0);
      isInitialized = true;
      //    }


      PS::F64 time_offset = PS::GetWtime();
      for(int i=0;i<n_loc;i++){
	system_janus[i].IntegrateBeforeForceCalc(dt,boxdh);
      }
      time_integ += PS::GetWtime() - time_offset;


      if(s%10 == 0){
	dinfo.collectSampleParticle(system_janus);
	dinfo.decomposeDomain();
    }
    //    std::cout<<"test2"<<std::endl;
    system_janus.exchangeParticle(dinfo);
    //    std::cout<<"test3"<<std::endl;
    n_loc = system_janus.getNumberOfParticleLocal();
    //    std::cout<<"test4"<<std::endl;
    //#ifdef DISSIPATIVE_RANDOM
    //    for(int i=0;i<system_janus.getNumberOfParticleLocal();i++) system_janus[i].seed = rand();
    //#endif

    tree_janus.calcForceAllAndWriteBack
      (CalcForceEpEp(
#ifdef DISSIPATIVE_RANDOM
		     dt,temperature
#endif
		     ),  system_janus, dinfo);
    //    std::cout<<"test5"<<std::endl;
    PS::Comm::barrier();
#if 1
    time_offset = PS::GetWtime();
    for(int i=0;i<system_janus.getNumberOfParticleLocal();i++){
      CalcIntraForce(system_janus[i],tree_janus,boxdh);
    }
    time_intra += PS::GetWtime() - time_offset;
#endif
    PS::Comm::barrier();

    time_offset = PS::GetWtime();
    for(int i=0;i<n_loc;i++){
      system_janus[i].IntegrateAfterForceCalc(dt);
    }
    time_integ += PS::GetWtime() - time_offset;

    //for calc momentum=0
    //checkmomentum0(system_janus);
#ifdef CALC_PRESSURE
  const PS::S32 n_loc = system_janus.getNumberOfParticleLocal();
  PS::F64vec pressure_loc = 0.0;
  PS::F64vec pressure_gro = 0.0;
  for(PS::S32 k = 0; k < n_loc; k++){
    pressure_loc.x += system_janus[k].pressure_kinetic[0] + 0.5 * system_janus[k].pressure_rf[0];
    pressure_loc.y += system_janus[k].pressure_kinetic[4] + 0.5 * system_janus[k].pressure_rf[4];
    pressure_loc.z += system_janus[k].pressure_kinetic[8] + 0.5 * system_janus[k].pressure_rf[8];
  }

  pressure_loc /= volume;
  MPI::COMM_WORLD.Allreduce(&pressure_loc.x, &pressure_gro.x, 1, PS::GetDataType<PS::F64>(), MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&pressure_loc.y, &pressure_gro.y, 1, PS::GetDataType<PS::F64>(), MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&pressure_loc.z, &pressure_gro.z, 1, PS::GetDataType<PS::F64>(), MPI::SUM);

  if(rank == 0){
    fpo3 = fopen( "pressure.dat" , "a" );
    if( fpo3 == NULL ) {
      printf( "ファイルオープンエラー\n" );
      return -1;
    }

    fprintf(fpo3,"%lf %lf %lf\n",pressure_gro.x,pressure_gro.y,pressure_gro.z);
    fclose(fpo3);
  }
#endif




    CalcEnergy(system_janus, Etot1, Etra1, Erot1, Epot1);
    PS::F64vec mom;
    CalcMomentum(system_janus, mom);
    if(s>=0){
      Epot_ave += Epot1;
      Ekin_ave += Etra1 + Erot1;
    }

    if(s == time_diag) {
      if(PS::Comm::getRank() == 0){
	fout_eng<< " " << time_sys
		<< " " << Epot1
		<< " " << Etra1
		<< " " << Erot1;
	fout_eng << " " << (Etot1-Etot0)/Etot0
		 << " " << mom.x
		 << " " << mom.y
		 << " " << mom.z
		 << std::endl;
	/*
	fprintf(stderr, "%10.7f %lf %lf %lf %+e\n",
		time_sys, Epot1, Etra1, Erot1, (Etot1 - Etot0) / Etot0);
	*/

	time_diag += nstep_diag;
      }

    }
    time_sys += dt;
    /*
    //add by kiyoshiro{
    fpo1<<system_janus.getNumberOfParticleGlobal()<<std::endl;
    fpo1<<s<<std::endl;
    for(PS::S32 i=0;i<system_janus.getNumberOfParticleGlobal();i++){
      fpo1<<"H"<<" "<<system_janus[i].pos.x<<" "<<system_janus[i].pos.y<<" "<<system_janus[i].pos.z<<std::endl;
    }
    //}
    */
  }

  const PS::TimeProfile tp = system_janus.getTimeProfile() + dinfo.getTimeProfile() + tree_janus.getTimeProfile();      
  
  printf("TotalTime= %lf\n\n",tp.getTotalTime() + time_integ + time_intra);
  
  printf("kick_and_drift= %lf\n",time_integ);
  printf("calc_intra= %lf\n", time_intra);
  tp.dump();
  PS::Finalize();
  return 0;
}
