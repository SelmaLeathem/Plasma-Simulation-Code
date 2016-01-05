/*
 * File: collisions.cpp
 * ---------------------
 * This file implements the collisions.h interface.
 *
 * 6/15
 */

#include "collisions.h"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>

collisions::collisions(int if_electron_in,double dt_in):if_electron(if_electron_in),dt(dt_in)
{

   neutral_vx = new double[INIT_NO_NEUTRALS]();
   if (!neutral_vx){
      printf("Error allocating memory for neutral_vx\n");
      exit(1);
   }
   
   neutral_vy = new double[INIT_NO_NEUTRALS]();
   if (!neutral_vy){
      printf("Error allocating memory for neutral_vy\n");
      exit(1);
   }

   neutral_vz = new double[INIT_NO_NEUTRALS]();
   if (!neutral_vz){
      printf("Error allocating memory for neutral_vz\n");
      exit(1);
   }
}

collisions::~collisions()
{
    delete [] neutral_vx;
    delete [] neutral_vy;
    delete [] neutral_vz;
}

/* 
 * Implementation Notes: calc_neutral_velocity
 * --------------------------------------------
 * This is a place-holder for a more accurate way of taking
 * velocities off a Maxwellian. 
 *
 * Arrays of neutral velocities are calculated before the simulation starts.
 * Neutral velocities are randomly chosen from these arrays during the 
 * simulation in collision calculations.
 */
void collisions::calc_neutral_velocity(void)
{
   int i;
   double R1,R2;

   for (i=0; i< INIT_NO_NEUTRALS; i++)
   {
    do{
       R1 = (double)rand()/((double)RAND_MAX);
       R2 = (double)rand()/((double)RAND_MAX);

       neutral_vx[i] = vtn*cos(PI*R2)*sqrt(fabs(log(R1)));
    }while(fabs(neutral_vx[i]) > 10000.0*vtn );

   do{

      R1 = (double)rand()/((double)RAND_MAX);
      R2 = (double)rand()/((double)RAND_MAX);
      neutral_vy[i] = vtn*cos(PI*R2)*sqrt(fabs(log(R1)));
      neutral_vz[i] = vtn*sin(PI*R2)*sqrt(fabs(log(R1)));
    }while(fabs(neutral_vy[i]) > 10000.0*vtn || fabs(neutral_vz[i]) > 10000.0*vtn);
   }

}

/*
 * Implementation Notes: scattered_velocity
 * ------------------------------------------
 * Calculates the post-collision velocites after a scattering event, as 
 * presented in:
 *  -Kenichi Nanbu,IEEE Trans. Plasma Sci.,vol 28,(2000)971
 *  -Tomonori Takizuka et al,Journal of Computational Physics,vol 25,(1977)205
 * 
 * This method is called for elastic scattering postcollision velocities, and
 * post ionization scattering.
 *
 * Note: Energy is in eV.
 */
void collisions::scattered_velocity(int i,vector<particleVec> &particleObj, double energy, double coll_factor,int yes_or_no)
{
   int index;
   double R,cosChi,sinChi,cosPhi,sinPhi;
   double nvx=0.0,nvy=0.0,nvz=0.0;
   double ux,uy,uz,u,uperp;
   double two_pi = 2.0*PI;
   double min_vel = 1.0;
   double minEnergy=1.0e-2;

   R = (double)rand()/((double)RAND_MAX);

   if (energy < minEnergy)
      cosChi = (double)if_electron*(1.0/(4.0*PI)) + (double)(1-if_electron)*sqrt(1.0 - R);
   else 
      cosChi = (double)if_electron*(2.0 + energy - 2.0*pow((1.0 +energy),R))/energy + (double)(1-if_electron)*sqrt(1.0 - R);

   sinChi = sqrt(1.0 - cosChi*cosChi);

   cosPhi = cos(two_pi*R);
   sinPhi = sin(two_pi*R);

   index= rand()%(INIT_NO_NEUTRALS - 1);
  
   nvx = neutral_vx[index]*((double)yes_or_no);
   nvy = neutral_vy[index]*((double)yes_or_no);
   nvz = neutral_vz[index]*((double)yes_or_no);

   ux = particleObj[i].particleVec::v - nvx;
   uy = particleObj[i].particleVec::vy - nvy;
   uz = particleObj[i].particleVec::vz - nvz;

   uperp = sqrt(uy*uy + uz*uz);

  if(uperp > min_vel){

    u = sqrt( ux*ux + uperp*uperp);

    particleObj[i].particleVec::vz += coll_factor*(-ux*(uz/uperp)*sinChi*cosPhi - (uy/uperp)*u*sinChi*sinPhi  + uz*(1.0 - cosChi));

    particleObj[i].particleVec::vy += coll_factor*( -(uy/uperp)*ux*sinChi*cosPhi
      - (uz/uperp)*u*sinChi*sinPhi + uy*(1.0 - cosChi));

    particleObj[i].particleVec::v += coll_factor*( uperp*sinChi*cosPhi
       + ux*(1.0 - cosChi));
  }

  else{
    u = ux;

    particleObj[i].particleVec::vz += coll_factor*(-u*sinChi*cosPhi);
    particleObj[i].particleVec::vy += coll_factor*( -u*sinChi*sinPhi);
    particleObj[i].particleVec::v += coll_factor*( u*(1.0-cosChi));
  }

}
