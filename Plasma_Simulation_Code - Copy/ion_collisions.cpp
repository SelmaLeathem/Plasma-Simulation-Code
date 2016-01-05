/*
 * File: ion_collisions.cpp
 * -----------------------
 * This file implements the ion_collisions.h interface. User defined collision
 * cross-section methods are declared in "ion_coll_input_declare" and user
 * defined cross-section methods are implemented in "ion_coll_input_methods".
 *
 * 6/15
 *
 */

#include "ion_collisions.h"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>


ion_collisions::ion_collisions(int if_electron_in, double dt_in, bool CM_elastic, bool CM_CX):collisions(if_electron_in,dt_in),center_mass_El(CM_elastic),center_mass_CX(CM_CX)
{
}

ion_collisions::~ion_collisions()
{
}

void ion_collisions::set_neutral_density(double density){
      neutral_density=density;
}

void ion_collisions::set_vtn(double neutral_velocity){
      vtn=neutral_velocity;
}

void ion_collisions::set_ion_mass(double mass){
     ion_mass= mass;
}

void ion_collisions::set_CMtrue_elastic(bool val){
     center_mass_El=val;
}

void ion_collisions::set_CMtrue_CX(bool val){
     center_mass_CX=val;
}

void ion_collisions::set_crossSec_elastic(double (ion_collisions::*func1)(void)){
     sigma_elastic = func1;
}

void ion_collisions::set_crossSec_CX(double (ion_collisions::*func1)(void)){
     sigma_CX= func1;
}

void ion_collisions::set_constant_elastic(double val){
     constant_elastic =val;
}

void ion_collisions::set_constant_CX(double val){
     constant_cx =val;
}

void ion_collisions::calc_CM_energy(double vx,double vy,double vz){
     double ux,uy,uz;

     ux=(ion_mass*vx+ion_mass*vtn);
     uy=(ion_mass*vy+ion_mass*vtn);
     uz=(ion_mass*vz+ion_mass*vtn);

     CMenergy=(0.5/CHARGE)*(1.0/(2.0*ion_mass))*(ux*ux +uy*uy +uz*uz);
}

void ion_collisions::calc_energy(double vx,double vy,double vz){
     double u2;

     u2 = vx*vx + vy*vy +vz*vz;
     energy = (0.5/CHARGE)*ion_mass*u2;
}

/*
 * Implementation Notes: sigma_const_elastic
 * --------------------------------------------
 * For a constant value or no collision need to set constant_elastic to a value
 * or zero before passing this method to the set_crossSec_elastic method.
 *
 * e.g.
 * object.set_constant_elastic(0.0);
 * object.set_crossSec_elastic(&ion_collisions::sigma_const_elastic);
 */
double ion_collisions::sigma_const_elastic(void){
      return constant_elastic;
}

/*
 * Implementation Notes: sigma_const_cx
 * --------------------------------------
 * For a constant value or no collision need to set constant_cx to a value
 * or zero before passing this method to the set_crossSec_CX method.
 *
 * e.g.
 * object.set_constant_CX(value);
 * object.set_crossSec_CX(&ion_collisions::sigma_const_cx);
 */
double ion_collisions::sigma_const_cx(void){
      return constant_cx;
}

/*
 * Implementation Notes: sigma_Ar_elastic
 * ----------------------------------------
 * Calculates the ion-atom energy dependent cross-section with formulae
 * pinched from Berkeley's pdp1 code.
 */
double ion_collisions::sigma_Ar_elastic(void){
  
   if (energy > 4.0)
       return (1.8e-19 + 4.0e-19/sqrt(energy));

   return(-2.0e-19*sqrt(energy) + 7.8e-19);
}

/*
 * Implementation Notes: sigma_Ar_cx
 * ----------------------------------------
 * Calculates the ion-atom energy dependent cross-section with formulae
 * pinched from Berkeley's pdp1 code.
 */
double ion_collisions::sigma_Ar_cx(void){
   
   if(energy > 4.0)
      return(2.0e-19 + 5.5e-19/sqrt(energy));

   return(-2.95e-19*sqrt(energy) + 10.65e-19);
}

/* Include user defined cross-section method implementations */
  #include "input_files/ion_coll_input_methods"

/*
 * Implementation Notes: charge_exchange_velocity
 * ------------------------------------------------
 * Implement charge_exhange collision by giving ions a post-collision velocity
 * equall to the neutral incident speed as described in:
 * -V Vahedi et al,Computer Physics Communications,vol 87,(1995)179
 */
void ion_collisions::charge_exchange_velocity(int i, vector<particleVec> &particleObj)
{
   int index;
   double R;

   R = (double)rand()/((double)RAND_MAX);
   index= rand()%(INIT_NO_NEUTRALS - 1);

   particleObj[i].particleVec::v = neutral_vx[index];
   particleObj[i].particleVec::vy = neutral_vy[index];
   particleObj[i].particleVec::vz = neutral_vz[index];

}

/*
 * Implementation Notes: collision
 * ---------------------------------
 * Determines if a collision occurs, and calls the collision methods if it
 * does. This is done in a few steps.
 *
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
 *  the collision is called with the coll_factor set to -reduced mass
 *  
 *  See Kenichi Nanbu and Takizuka for more information about the coll_factor.
 *  -Tomonori Takizuka et al,Journal of Computational Physics,vol 25,(1977)205
 *
 * Need to verify that have implemented tests to find out if collision occurs
 * correctly.
 *
 * Currently implement an elastic collision with an atom by scattering the ion
 * using the scattering formulae in Nanbu and Takizuka which serves as a 
 * temporary place holder until this is properly looked into.
 *
 */
void ion_collisions::collision(vector<particleVec> &particleObj)
{
   int count_iscat=0,count_cx=0;
   int i=0, particle_no;
   double probability;
   double sigma_iEl, sigma_iCx, sigma_total;
   double R;
   double u;
   double reduced_mass,coll_factor;
   double Pmax=1.0; // 1.0 - exp(-1.0);

   particle_no = particleObj.size();

  /*treat ion mass as same wieght as atom */
   reduced_mass = 1.0/2.0;

   while(i < particle_no)
   {
     if(particleObj[i].x >=0.0 && particleObj[i].x < LENGTH){

      calc_energy(particleObj[i].v,particleObj[i].vy,particleObj[i].vz);

      u= sqrt( (2.0*CHARGE*this->energy)/ion_mass);

      if( center_mass_CX==true || center_mass_El==true)
         calc_CM_energy(particleObj[i].v,particleObj[i].vy,particleObj[i].vz);

      sigma_iEl = (this->*sigma_elastic)();
      sigma_iCx = (this->*sigma_CX)();

     /* Charge exchange collisions */

        R = ((double)(rand()%10000000))/1.0e7;
        probability = 1.0 -exp(-(u*dt*neutral_density*sigma_iCx))/Pmax;

      if(R < probability){
               charge_exchange_velocity(i,particleObj);
               count_cx += 1;  //used to monitor the number of collisions
      } /*goes with if R > */

/* 
 * Currently implement an elastic collision with an atom by scattering the ion
 * using the scattering formulae in Nanbu and Takizuka which serves as a 
 * temporary place holder until this is properly looked into.
 */
       R = ((double)(rand()%10000000))/1.0e7;
       probability = 1.0 -exp(-(u*dt*neutral_density*sigma_iEl))/Pmax;

      if(R < probability){
               coll_factor = -reduced_mass;
               scattered_velocity(i,particleObj,this->energy,coll_factor,YES);
               count_iscat += 1;  //used to monitor the number of collisions
      } /*goes with if R > */
      
      } /* if x within boundaries */
      i+=1;
   } /* goes with while loop */

}
