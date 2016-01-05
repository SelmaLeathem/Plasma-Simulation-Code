/*
 * File: particle.cpp
 * --------------------
 * This file implements the particle.cpp interface.
 *  
 * 4/15
 *
 */

#include "particle.h"
#include <algorithm>
#include <math.h>
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <fstream> 
#include <string.h>


particle::particle(double dt_p, double vt_p, double q_over_m_p, double qc_p,double mass_p,double dxp_p):dt(dt_p), vt(vt_p), q_over_m(q_over_m_p), qc(qc_p),mass(mass_p),dxp(dxp_p)
{
   count_reduce = 0;
   int i;
  
   q = new double[NO_OF_CELLS +1]();
   if(!q){
     printf("Error allocating memory for q\n");
     exit(1);
   }

   for(i=0;i<(NO_OF_CELLS+1);i++)
      *(q+i) = 0.0;

   v_limit = new double[NO_GROUPS]();
   if(!v_limit){
     printf("Error allocating memory for v_limit\n");
     exit(1);
   }

   group_limit = new int[NO_GROUPS]();
   if(!group_limit){
     printf("Error allocating memory for group_limit\n");
     exit(1);
   }

   new_group_size = new int[NO_GROUPS]();
   if(!new_group_size){
     printf("Error allocating memory for new_group_size\n");
     exit(1);
   }

   group_weight = new double[NO_GROUPS]();
   if(!group_weight){
     printf("Error allocating memory for group_weight\n");
     exit(1);
   }

   j_group = new int[NO_GROUPS]();
   if(!j_group){
     printf("Error allocating memory for j_group\n");
     exit(1);
   }

   Ngroup = new int[NO_GROUPS]();
   if(!Ngroup){
     printf("Error allocating memory for Ngroup\n");
     exit(1);
   }

   Nremove = new int[NO_GROUPS]();
   if(!Nremove){
     printf("Error allocating memory for Nremove\n");
     exit(1);
   }

}

particle::~particle()
{
   delete [] q;
   delete [] v_limit;
   delete [] group_weight;
   delete [] group_limit;
   delete [] new_group_size;
   delete [] j_group;
   delete [] Ngroup;
   delete [] Nremove;

}

void particle::set_vt(double vt){
     this->vt=vt;
}

void particle::set_q_over_m(double q_over_m){
     this->q_over_m=q_over_m;
}

void particle::set_mass(double mass){
     this->mass=mass;
}

void particle::set_dxp(double dxp){
     this->dxp=dxp;
}

void particle::set_init_no(int init_particle_no){
     this->init_particle_no=init_particle_no;
}

void particle::set_no_inject(int no_to_inject){
     this->no_inject = no_to_inject;
}

void particle::set_dt(double dt){
     this->dt=dt;
}

void particle::set_qc(double qc){
     this->qc=qc;
}
 
void particle::set_no_groups(int val){
     this->no_groups = val;
}

void particle::set_v_limit(double val,int element){
     v_limit[element]=val;
}

void particle::set_group_limit(int val,int element){
     group_limit[element]=val;
}

void particle::set_new_group_size(int val,int element){
     new_group_size[element]=val;
}

void particle::set_size_to_check(int val){
     this->size_to_check=val;
}

int particle::get_init_no(void){
     return this->init_particle_no;
}

void particle::randomize_vector(void){
     random_shuffle(particleObj.begin(), particleObj.end());
}

void particle::create_particles(void){
   int i;
   for(i=0;i<init_particle_no;i++){
       particleObj.push_back(particleVec());
   }
}

void particle::is_equall_to(vector<particleVec> &obj1, vector<particleVec> &obj2,int i, int j){

    obj1[i].x = obj2[j].x;
    obj1[i].v = obj2[j].v;
    obj1[i].vy = obj2[j].vy;
    obj1[i].vz = obj2[j].vz;
    obj1[i].weight = obj2[j].weight;
}


void particle::load_particles(void){
     int i=0;
     double x,vx,vy,vz,weight;
     ifstream in;

     if(!in){
       cout<<"Unable to open file"<<endl;
       exit(1);
     }

     char filename[80]="dataDump_";
     strcat(filename,type);
     strcat(filename,".txt"); 
  
     in.open(filename);

     while(!in.eof()){
       in>>x>>vx>>vy>>vz>>weight;
       particleObj.push_back(particleVec(x,vx,vy,vz,weight));
    }

   in.close(); 
}

/*
 * Implementation Notes: initiate_positions
 * -------------------------------------------  
 * Particles are placed evenly along the x-axis at a distance dxp from each 
 * other.
 */
void particle::initiate_positions(void)
{
  int i;
  for (i=0; i< init_particle_no; i++)
     particleObj[i].x = ((double)i)*dxp;
}

/*
 * Implementation Notes: get_particle_velocity
 * ---------------------------------------------
 * This method is a temporary place holder for a more accurate way to 
 * calculate the particle velocity taken randomly off a particle velocity
 * distribution. In this case Maxwellian.
 */
double particle::get_particle_velocity(void)
{
    double R1, R2;
    double v = 0.0;

    do{
       R1 =(double)rand()/((double)RAND_MAX);
       R2 =(double)rand()/((double)RAND_MAX);
         v = vt*cos(2.0*PI*R2)*sqrt(fabs(log(R1)));
    }while( fabs(v) > vt*10000);

    return v;
}

/*place holder for more accurate 3d method*/
void particle::initiate_velocities(void)
{
   int i;

   for(i=0; i< init_particle_no; i++){
      particleObj[i].v = get_particle_velocity();
      particleObj[i].vy = get_particle_velocity();
      particleObj[i].vz = get_particle_velocity();
   } 
}

void particle::inject_more_particles(void){
   int i;
   double random_number; 

   for(i=0;i<no_inject;i++){
       particleObj.push_back(particleVec());
       random_number =(double)rand()/((double)RAND_MAX);
       particleObj[(particleObj.size()-1)].x = PUT_PARTICLE_AT;
       particleObj[(particleObj.size()-1)].v = get_particle_velocity();
       particleObj[(particleObj.size()-1)].vy = get_particle_velocity();
       particleObj[(particleObj.size()-1)].vz = get_particle_velocity();
   }
}

void particle::initialize_map(void){
   int i;

   for(i=0;i<no_groups;i++){
      j_group[i] = 0;
      Ngroup[i] = 0;
      Nremove[i] = 0;
      group_weight[i]=1.0;
   }
}

/* Implementation Notes: reduce_population
 * -----------------------------------------
 * See the implementation notes for the contain_particle_growth method for
 * more information.
 */
void particle::reduce_population(void){
   int i,g, n = particleObj.size();
   int particles_to_remove=0;
   int R;
   int paws=0;

  /* Calculate the number of particles in each velocity group */
   for(i=0; i< n; i++){

   /* Ngroup holds the number of particles in each velocity group */
       g=0; 
       while( fabs(particleObj[i].v) > v_limit[g] ) {
           g++;
       }
 
       Ngroup[g]+=1; 
   }

   /* Calculate Nremove, the number of particles to delete from each velocity
    * group.
    */
   for(g=0;g<no_groups;g++){
      if(Ngroup[g] > group_limit[g])
            Nremove[g]= Ngroup[g] - new_group_size[g];
   }

   for(g=0;g<no_groups;g++){
   /* The total number of particles to remove from the particle vector */
      particles_to_remove+=Nremove[g];
   }
    
    i=0;
   
    while( i < particles_to_remove){

     /* Randomly select a particle from the particle vector */
      R =rand()%(particleObj.size()-1);
                 
      /* Determine which velocity group the particle vector element at R 
         belongs to. */ 
      g=0; 
      while( fabs(particleObj[R].v) > v_limit[g] ){
            g++;
      }

      /* j_group is a counter that keeps track of the number of particles
       * deleted from that group 
       */
      /* Remove vector particle element R if j_group < Nremove */
      if( j_group[g] < Nremove[g] ){
            j_group[g]++;
            swap(particleObj[R],particleObj[(particleObj.size()-1)]);
            particleObj.pop_back();
            i++;
      }
   }//while i<particles_to_remove

   /* Calculate the particle weight factor of each velocity group */

   if(particles_to_remove >0){
      for(g=0;g<no_groups;g++){
         if(j_group[g]>0){
          group_weight[g]= (double)Ngroup[g]/( (double)(Ngroup[g] - j_group[g]));
          }
      }
 
      n=particleObj.size();

   /* Increase particle wieght of all remaining particles */

      for(i=0;i<n;i++){
          g=0;
          while(fabs(particleObj[i].v) > v_limit[g]){
            g++;
          }
        particleObj[i].weight = (particleObj[i].weight)*group_weight[g];
      }

  } //if particles_to_remove >0
}

/*
 * Implementation Notes: contain_particle_growth
 * -----------------------------------------------
 * Limit the size of the vector particle arrays by reducing the particle number
 * according to the scheme set out by E E Kunhardt et al:
 * -Journal of Computational Physics,vol 67,(1986) 279
 * -----------------------------------------------------
 * This scheme involves dividing the particles up into velocity groups,
 * and decreasing the size of each group once the upper group size is reached.
 * Each remaining particle is then rewieghted according to the amount that the 
 * group is reduced by. The number of particle velocity groups and the upper 
 * size limit of each group is determined by the user. 
 */
void particle::contain_particle_growth(void){
   if( particleObj.size() > size_to_check ){
        initialize_map();
        reduce_population();
        count_reduce++;
   }
}

/*
 * Implementation Notes: collect_charge
 * --------------------------------------
 * Calculates charge on the grid using the first order scheme presented in
 * C K Birdsall and A B Langdon, Plasma Physics via Computer Simulation, 
 * Institute of Physics Publishing, chapter 2.
 */
void particle::collect_charge(void)
{
   int i,k,index=0,jp,paws=0;
   double w=0.0;
   double particle_weight;
   double qc_dx = qc/(2.0*DX);
   int no_particles = particleObj.size();

   for(k=0; k< NO_OF_CELLS+1; k++){
       *(q+k) = 0.0;
   }

   for(i=0; i < no_particles; i++){

      jp = (int)(fabs(particleObj[i].x/DX));

      if(jp < NO_OF_CELLS && jp >= 0){
         w = particleObj[i].x/DX - (double)jp;
         particle_weight = particleObj[i].weight; 

         q[jp]+= particle_weight*(1.0 - w);
         q[jp+1]+= particle_weight*w;
      } //if jp
   }

   *(q+NO_OF_CELLS)=0.0;
   for(k=0; k< NO_OF_CELLS; k++){
      *(q+k) = qc_dx*(*(q+k));
   }

   *q = 2.0*(*q);
   *(q+NO_OF_CELLS - 1) = 2.0*(*(q+NO_OF_CELLS-1));
   *(q+NO_OF_CELLS) = 0.0;
}

/* 
 * Implementation Notes: collect_extraCharge
 * -------------------------------------------
 * Implements the same algorithm used in collect_charge but only iterates
 * over the newer ions, rather than the entire vector.
 */
void particle::collect_extraCharge(int old_size)
{
   int i,k,index=0,jp,paws=0;
   double w=0.0;
   double particle_weight;
   double qc_dx = qc/(2.0*DX);
   int current_size = particleObj.size();

   *(q+NO_OF_CELLS)=0.0;

   for(k=0; k< NO_OF_CELLS; k++){
       *(q+k) = (*(q+k))/qc_dx;
   }

   for(i=old_size; i < current_size; i++){

      jp = (int)(fabs(particleObj[i].x/DX));

      if(jp < NO_OF_CELLS && jp >= 0){
          w = particleObj[i].x/DX - (double)jp;
          particle_weight = particleObj[i].weight; 

          q[jp]+= particle_weight*(1.0 - w);
          q[jp+1]+= particle_weight*w;
      }
   }

   for(k=0; k< NO_OF_CELLS; k++){
       *(q+k) = qc_dx*(*(q+k));
   }
   *q = 2.0*(*q);
   *(q+NO_OF_CELLS - 1) = 2.0*(*(q+NO_OF_CELLS - 1));
   *(q+NO_OF_CELLS)=0.0;
}

/*
 * Implementation Notes: move_particles
 * ----------------------------------------
 * Calculates the new particle velocities and positions due to the electric
 * and optional magnetic fields using the leap-frog and  
 * accelerate-rotate-accelerate methods presented in:
 * C K Birdsall and A B Langdon, Plasma Physics via Computer Simulation, 
 * Institute of Physics Publishing, chapter 4.
 *    
 * The electric field on the grid is interpolated to the electric-field
 * at the particle according to Birdsall and Langdon, chapter 2.
 */                                                                                                                                               
void particle::move_particles(double *Egrid){
   int paws=0;
   int i=0,jp,w;
   int no_particles = particleObj.size();
   double coef = 0.5*dt*q_over_m;
   double omega_dt = dt*q_over_m*MAGNETIC_FIELD;
   double cos_angle = cos(ANGLE);
   double sin_angle = sin(ANGLE);
   double cos_omega = cos(omega_dt);
   double sin_omega = sin(omega_dt);
   double vperp,vpar,vperp_after,ux; //v perpendicular and parallel to B


    for(i=0;i<no_particles;i++){
      jp = (int)(fabs(particleObj[i].x/DX));

      if(jp < NO_OF_CELLS && jp >= 0){
         w = particleObj[i].x/DX - (double)jp;
  
      /* ux is used to monitor large velocity increases */
      ux=particleObj[i].v;


        /***  - 1/2 accelerate - ***/

        particleObj[i].v += coef*( (1.0 - w)*Egrid[jp] + w*Egrid[jp+1] );

        /*** -rotate around magnetic field- ***/

        if(fabs(MAGNETIC_FIELD) > 0.0 ){

          /*perpendicular and parallel velocities*/

          vpar = cos_angle*particleObj[i].v + sin_angle*particleObj[i].vy;
          vperp = -sin_angle*particleObj[i].v + cos_angle*particleObj[i].vy;

         /* rotate perpendicular components */
      
          vperp_after=vperp*cos_omega +sin_omega*particleObj[i].vz;
          particleObj[i].vz= -vperp*sin_omega +cos_omega*particleObj[i].vz;

         /* get vx and vy from vperp_after and vpar */
   
          particleObj[i].v = cos_angle*vpar - sin_angle*vperp_after;
          particleObj[i].vy = sin_angle*vpar + cos_angle*vperp_after;
        }


        /***  - 1/2 accelerate - ***/

         particleObj[i].v += coef*( (1.0 - w)*Egrid[jp] + w*Egrid[jp+1] );

     /* Check for over-heating */    
      if(fabs(particleObj[i].v) > 7.9e7 && particleObj[i].x >=DX && 
                 particleObj[i].x <= (LENGTH-DX)){

           cout<<"in move, v>7.9e7,v before move= "<<ux<<" v after= "<<particleObj[i].v<<" x= "<<particleObj[i].x<<" jp= "<<jp<<" weight= "<<particleObj[i].weight<<endl;
           cout<<"i = "<<i<<" E= "<<Egrid[jp]<<" Ej+1 "<<Egrid[jp+1]<<endl; 
           cout<<"This is typically due to overheating. Need to ensure the"
                 " electron time-step is less than the 1/plasma_frequency. The" 
                 " electron density also needs to be less than the density" 
                 " used in the plasma frequency to help avoid over-heating."  
                " Things to try include reducing the macro-particle size and/or"
                 " the number of allowed ionization events."<<endl;
           cout<<"Hit enter to move on. "<<endl;
           cin>>paws; 
      }

        /* Move particle along the x-axis */
        particleObj[i].x += dt*particleObj[i].v;

    } //if jp
  } //for i
}

void particle::set_type(char *type_in){
     type = type_in;
}

void particle::dump_data(void){
     int i;
     ofstream out;

     if(!out){
       cout<<"Unable to open file"<<endl;
       exit(1);
     }

     char filename[80]="dataDump_";
     strcat(filename,type);
     strcat(filename,".txt"); 
  
     out.open(filename);

     for(i=0;i<particleObj.size();i++){
       out<<particleObj[i].x<<" "<<particleObj[i].v<<" "<<particleObj[i].vy<<" "<<particleObj[i].vz<<" "<<particleObj[i].weight<<endl;
    }

    out.close();
}

