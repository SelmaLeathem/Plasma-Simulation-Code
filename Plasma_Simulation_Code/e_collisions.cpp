/*
 * File: e_collisions.cpp
 * -----------------------
 * This file implements the e_collisions.h interface. User defined collision
 * cross-section methods are declared in "e_coll_input_declare" and user
 * defined cross-section methods are implemented in "e_coll_input_methods".
 *
 * 6/15
 *
 */

#include "e_collisions.h"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>


e_collisions::e_collisions(int if_electron_in,double dt_in):collisions(if_electron_in,dt_in)
{
}

e_collisions::~e_collisions()
{
}

void e_collisions::set_neutral_density(double density){
      neutral_density=density;
}

void e_collisions::set_vtn(double neutral_velocity){
      vtn=neutral_velocity;
}

void e_collisions::set_ion_mass(double mass){
      ion_mass= mass;
}

void e_collisions::set_constant_elastic(double val){
      constant_elastic = val;
}

void e_collisions::set_constant_ionization(double val){
      constant_ionization = val;
}

void e_collisions::set_crossSec_elastic(double (e_collisions::*func1)(void)){
      sigma_elastic = func1;
}

void e_collisions::set_crossSec_ionization(double (e_collisions::*func1)(void)){
      sigma_ionization= func1;
}

void e_collisions::set_ionization_threshold(double val){
      threshold = val;
}

void e_collisions::calc_energy(double vx,double vy,double vz){
     double u2;
     u2 = vx*vx + vy*vy +vz*vz;
     energy = (0.5/CHARGE)*ELECTRON_MASS*u2;
}

/*  
 * Implementation Notes: sigma_const_elastic
 * -------------------------------------------- 
 * For a constant value or no collision need to set constant_elastic to a value
 * or zero before passing this method to the set_crossSec_elastic method.
 *
 * e.g. 
 * object.set_constant_elastic(0.0);
 * object.set_crossSec_elastic(&e_collisions::sigma_const_elastic);
 */
double e_collisions::sigma_const_elastic(void){
        return constant_elastic;
}

/*  
 * Implementation Notes: sigma_const_ionization
 * --------------------------------------------- 
 * For a constant value or no collision need to set constant_ionization to a 
 * value or zero before passing this method to the set_crossSec_ionization
 * method.
 *
 * e.g. 
 * object.set_constant_ionization(1.0e-20);
 * object.set_crossSec_ionization(&e_collisions::sigma_const_ionization)
 */
double e_collisions::sigma_const_ionization(void){
        return constant_ionization;
}

/*
 * Implementation Notes: sigma_Ar_elastic
 * ---------------------------------------
 * Calculates the electron-Argon energy dependent cross-section for
 * energies in the range: threshold to 6 keV. The formulae are a collection
 * excel data fits to data presented in:
 * -Ashok Jain et al,Physica Scripta,vol 41,(1990)321
 */
double e_collisions::sigma_Ar_elastic(void){
   double sigma;

   if ( energy <= 0.35 )
  
      sigma = 28019.0*pow(energy,6.0)- 32798.0*pow(energy,5.0)
              + 15048.0*pow(energy,4.0) - 34312.0*pow(energy,3.0)
              + 4084.0*pow(energy,2.0) - 249.0*energy + 7.576;

   else if( energy <= 15.0 && energy > 0.35)

      sigma = - 0.013*pow(energy,3.0)
              + 0.232*pow(energy,2.0) + 1.066*energy + 0.085;

   else if( energy  > 15.0 ) //stictly only valid to 6k

      sigma = 126.2*pow(energy,-0.69);

   return (fabs(sigma*1.0e-20));
}

/*
 * Implementation Notes: sigma_Ar_ionize
 * ---------------------------------------
 * Calculates the electron-Argon energy dependent cross-section for
 * energies in the range: threshold to 6 keV. The formulae are a collection
 * excel data fits to data presented in:
 * -Ashok Jain et al,Physica Scripta,vol 41,(1990)321
 *
 * Decrease sigma by factor of 0.1 since currently having issues
 * with excessive ionization in the sheaths causing large localized
 * electric fields.
 */
double e_collisions::sigma_Ar_ionize(void){
   double sigma;

   if(energy < 15.755)
      sigma = 0.0;
   
   else if (energy <=100.0 && energy >= 15.755)
      
       sigma = (-3.0e-7)*pow(energy,4.0)+ (8.0e-5)*pow(energy,3.0)
             - 0.008*pow(energy,2.0) + 0.416*energy - 4.758;

   else if (energy > 100.0 )  //only valid to 6k

      sigma = 130.0*pow(energy,-0.74);

   return fabs(sigma*1.0e-20);
}

  /* User defined collision cross-section methods: */
   #include "input_files/e_coll_input_methods"

/*
 * Implementation Notes: ionization 
 * ----------------------------------
 * Implements ionization collisions using the method outlined in :
 *  -Kenichi Nanbu et al,IEEE Trans. Plasma Sci.,vol 28,(2000) 971
 *   
 * In addition to calculating the post-collision velocities, this method
 * adds the new particles to the particle vectors.
 */   
void e_collisions::ionization(int i, vector<particleVec> &particleObj_e, vector<particleVec> &particleObj_i,double energy, double coll_factor)
{
   int index,paws=0;
   int n;
   double weight = particleObj_e[i].weight; 
   const double a = 10.3;
   double eta0,eta1,ejected_energy,R,Wx,Wy,Wz;
   double v_before, vy_before, vz_before;
   double v_temp_ion, vy_temp_ion, vz_temp_ion;
   double nvx = 0.0, nvy = 0.0, nvz = 0.0;
   double ve=0.0,vye=0.0,vze=0.0;
   double xe=particleObj_e[i].x;
   double energyFacOrig,Efactor;
   double energyFacNew,EfactorNew;
   double energyOrig,energyNew,u;

   eta0 = 2.0 - 100.0/(energy+10.0);
   eta1 = (energy - threshold)/2.0 - eta0;
   R = (double)rand()/((double)RAND_MAX);

   ejected_energy = eta0+a*tan(R*(atan(eta1/a)+atan(eta0/a))-atan(eta0/a)); 

  /*The neutral velocities are chosen randomly from an array that is 
   *populated before the simulation starts.
   */
   index= rand()%(INIT_NO_NEUTRALS - 1);

   nvx = neutral_vx[index];
   nvy = neutral_vy[index];
   nvz = neutral_vz[index];
  
   ve = particleObj_e[i].particleVec::v;
   vye = particleObj_e[i].particleVec::vy;
   vze = particleObj_e[i].particleVec::vz;
  
  /* Center of mass velocity */
   Wx = (ion_mass*nvx+ELECTRON_MASS*ve)/(ion_mass+ELECTRON_MASS);
   Wy = (ion_mass*nvy+ELECTRON_MASS*vye)/(ion_mass+ELECTRON_MASS);
   Wz = (ion_mass*nvz+ELECTRON_MASS*vze)/(ion_mass+ELECTRON_MASS);

   Efactor=((1.0-(threshold+ejected_energy)/energy));
   if(Efactor < 0.0) cout<<" sqrt of neg number in ionize "<<endl;
   energyFacOrig = sqrt(Efactor);

   v_before = ve*energyFacOrig;
   vy_before = vye*energyFacOrig;
   vz_before = vze*energyFacOrig;


   particleObj_e[i].particleVec::v = v_before;
   particleObj_e[i].particleVec::vy = vy_before;
   particleObj_e[i].particleVec::vz = vz_before;
   
   energyOrig = energy - (threshold + ejected_energy);

  /*Apply equations 56 to calculate the impacting electron's post-collision
   velocity*/ 

   //scattered_velocity(i,particleObj_e, ejected_energy ,coll_factor);
   scattered_velocity(i,particleObj_e, ejected_energy ,coll_factor,YES);
   
   EfactorNew = ejected_energy/energy;
   if(EfactorNew < 0.0) cout<<" sqrt of neg number New in ionize "<<endl;
   energyFacNew = sqrt(EfactorNew);

   v_before=ve*energyFacNew;
   vy_before=vye*energyFacNew;
   vz_before=vze*energyFacNew;

   particleObj_e.push_back(particleVec(xe,v_before,vy_before,vz_before,weight));

   n=particleObj_e.size();

  /*Apply equations 56 to calculate the new electron's post-collision velocity*/

   //scattered_velocity((n-1),particleObj_e,energyOrig,coll_factor);
   scattered_velocity((n-1),particleObj_e, energyOrig,coll_factor,YES);

  /*Calculate the new ion's postcollision velocity*/ 

  v_temp_ion= Wx -(ELECTRON_MASS/ion_mass)*((particleObj_e[i].v-Wx) + (particleObj_e[n-1].v - Wx) );

  vy_temp_ion= Wy -(ELECTRON_MASS/ion_mass)*((particleObj_e[i].vy-Wy) + (particleObj_e[n-1].vy - Wy) );

  vz_temp_ion= Wz-(ELECTRON_MASS/ion_mass)*((particleObj_e[i].vz-Wz) + (particleObj_e[n-1].vz - Wz) );
  
particleObj_i.push_back(particleVec(xe,v_temp_ion,vy_temp_ion,vz_temp_ion,weight) );
   
}

/*
 * Implementation Notes: collision
 * ---------------------------------
 * Determines if a collision occurs, and calls the collision methods if it 
 * does. This is done in a few steps.
 *
 * Originally followed the following method:
 *1)Determines if a collision occurs by comparing a random number to the 
 *  total collision probability as is described in the following articles:
 *  -V. Vahedi et al,Computer Physics Communications,vol 87,(1995) 179
 *  -C K Birdsall, IEEE Trans. Plasma Sci.,vol 19,(1991) 65
 *  -Kenichi Nanbu et al,IEEE Trans. Plasma Sci.,vol 28,(2000) 971
 *
 *2)If a collision occurs the type of collision is determined by comparing 
 *  each individual cross-section to random numbers as described in the
 *  articles listed in step one.
 *
 *3)If a collision occurs the relevant collision method to implement the 
 *  the collision is called with the coll_factor set to the negative
 *  of the reduced mass.
 *  See Kenichi Nanbu and Takizuka for more information about the coll_factor. 
 *  -Tomonori Takizuka et al,Journal of Computational Physics,vol 25,(1977)205
 *
 * 
  *The code which decides whether or not an electron collision occurs
  *is still a work in progress due to the need to limit the number of
  *ionization events. Currently an elastic collision occurs if a random
  *number between zero and one is less than the collision probability and
  *likewise for ionization. In the case of ionization a further condition
  *that the number of ionization events is less than a user defined limit
  *is also enforced.
 * 
 */
void e_collisions::collision(vector<particleVec> &particleObj_e, vector<particleVec> &particleObj_i)
{
   int count_escat=0,count_ionize=0;
   int paws=0;
   int i=0, particle_no;
   double probability;
   double sigma_eEl, sigma_eIoniz; /*elastic and ionization cross-sections*/
   double sigma_total; /*sum of elastic and ionization cross-sections*/
   double R;
   double reduced_mass,u,coll_factor,ux;
   double offset = 1.0e-6;
   double Pmax=1.0;/*Maximum probability set to one since problem with too
                     many ionization events*/

   particle_no = particleObj_e.size();

   reduced_mass = ion_mass/(ELECTRON_MASS + ion_mass); 

   while(i< particle_no){

    if(particleObj_e[i].x >=0.0 && particleObj_e[i].x < LENGTH){

      calc_energy(particleObj_e[i].v,particleObj_e[i].vy,particleObj_e[i].vz);
 
      u= sqrt( (2.0*CHARGE*this->energy)/ELECTRON_MASS);
      
      ux = particleObj_e[i].v;

   /* Calculate collision cross-sections */
      sigma_eEl = (this->*sigma_elastic)();
      sigma_eIoniz = (this->*sigma_ionization)();
      sigma_total = sigma_eIoniz + sigma_eEl;
      sigma_total+=offset; 

   /* probability = (1.0 -exp(-(u*dt*neutral_density*sigma_total)))/Pmax;*/

      R = ((double)(rand()%10000000))/1.0e7;
      probability = (1.0 -exp(-(u*dt*neutral_density*sigma_eEl)))/Pmax;

    /*   
      if(R < probability){

        if( (sigma_eEl < sigma_eIoniz && R < sigma_eEl/sigma_total)
            || (sigma_eIoniz < sigma_eEl && R > sigma_eIoniz/sigma_total))
     */
         if( R < probability )
         {
               /*elastic collision*/
               coll_factor = -reduced_mass;
               scattered_velocity(i,particleObj_e,this->energy,coll_factor,YES);
               count_escat += 1;   //used to monitor the number of collisions
               
         }

      /*  else */

       R = ((double)(rand()%10000000))/1.0e7;
       probability = (1.0 -exp(-(u*dt*neutral_density*sigma_eIoniz)))/Pmax;
         
      /* In light of the fact that an upper limit on the number of ionization
       * events is used the electron vector is randomly shuffled once each
       * time-step
       */ 
       if( R < probability && count_ionize< IONIZE_LIMIT )
       {
                 /*ionizing collision*/
                 coll_factor = -reduced_mass;
                 ionization(i,particleObj_e,particleObj_i,this->energy,coll_factor);
                 count_ionize += 1;  //used to monitor the number of collisions
        } // if R < ionize probability

   /* } */ /*goes with if R > */

   /* Check for over-heating */
      if(fabs(particleObj_e[i].v) > 7.9e7 && particleObj_e[i].x >=DX &&
           particleObj_e[i].x <= (LENGTH-DX)){
               cout<<"in collisions, v>7.9e7,v before coll = "<<ux<<" v after= "<<particleObj_e[i].v<<" x= "<<particleObj_e[i].x<<" weight= "<<particleObj_e[i].weight<<endl;
           cout<<"energy = "<<energy<<endl; 
           cout<<"i = "<<i<<endl;
           cout<<"This is typically due to overheating. Need to ensure the"
                 " electron time-step is less than the 1/plasma_frequency. The"
                 " electron density also needs to be less than the density" 
                 " used in the plasma frequency to help avoid over-heating." 
                " Things to try include reducing the macro-particle size and/or"
                 " the number of allowed ionization events."<<endl;   
           cout<<"Hit enter to move on "<<endl;
           cin>>paws;
        }

     }/* if x within boundaries */
     i+=1;
   } /*goes with while loop*/

}
