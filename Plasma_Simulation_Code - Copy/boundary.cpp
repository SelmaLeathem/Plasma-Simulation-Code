/*
 * File: boundary.cpp
 * --------------------
 * This file implements the boundary.h interface.
 *
 * 11/15
 *
 */

#include "boundary.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>

boundary::boundary(double dt_in):dt(dt_in){

   wall_offset = 1.0e-6;

}

boundary::~boundary(){

}

void boundary::set_max_yield(double val){
       max_yield=val;
}


void boundary::set_E0(double val){
       E0=val;
}

void boundary::set_Emax0(double val){
       Emax0=val;
}

void boundary::set_incident_mass(double val){
       incident_mass=val;
}

void boundary::set_vt(double val){
       vt=val;
}

void boundary::set_threshold(double val){
       threshold=val;
}

void boundary::set_work_function(double val){
       work_function=val;
}

void boundary::set_ion_proportion(double val){
       ion_proportion=val/100.0;
}

void boundary::set_implement_e_boundary(void (boundary::*func1)(vector<particleVec> &electronObj, vector<particleVec> &ionObj)){
      implement_e_boundary = func1;
}

void boundary::set_implement_i_boundary(void (boundary::*func1)(vector<particleVec> &electronObj, vector<particleVec> &ionObj)){
      implement_i_boundary = func1;
}

void boundary::calc_energy(double vx,double vy,double vz){
     double u2;
     u2 = vx*vx + vy*vy +vz*vz;
     energy = (0.5/CHARGE)*incident_mass*u2;
}

void boundary::calc_angle(double vx){
    double offset=1.0e-6;
    double u2= 2.0*energy/incident_mass;
 
    angle= acos(vx/(sqrt(u2)+offset));
}

void boundary::set_e_impact_yield(void (boundary::*func1)(void)){
               e_impact_see_yield = func1;
}

void boundary::set_i_impact_yield(void (boundary::*func1)(void)){
                i_impact_see_yield = func1;
}

void boundary::set_reflection_yield(void (boundary::*func1)(void)){
                calc_reflection_yield=func1;
}

void boundary::set_backscatter_yield(void (boundary::*func1)(void)){
                calc_backscatter_yield=func1;
}

/* User defined boundary yield methods */
  #include "boundary_input_methods"

void boundary::reflect_electron(int current_particle,vector<particleVec> &particleObj){
   double R;
   R = (double)(rand())/((double)RAND_MAX);
   particleObj[current_particle].v = -particleObj[current_particle].v;
/*
  if (particleObj[current_particle].x > LENGTH_MINUS_DX){
      particleObj[current_particle].x = LENGTH_MINUS_DX -
      R*dt*fabs(particleObj[current_particle].v);
*/
  if (particleObj[current_particle].x > LENGTH){
      particleObj[current_particle].x = LENGTH -
      R*dt*fabs(particleObj[current_particle].v);
  }
  else if (particleObj[current_particle].x < 0.0)
      particleObj[current_particle].x = R*dt*fabs(particleObj[current_particle].v);
}

/*
 * Implementation Notes: secondary_vel
 * --------------------------------------
 * Calculates the position and velocity of backscattered and true secondary
 * electrons using the formulae presented in:
 *   -Ralf Krimke et al, J. Phys. D: Appl. Phys.,vol 29,(1996),378
 *
 */
void boundary::secondary_vel(int current_particle,double velocity,vector<particleVec> &particleObj){
   int sign=1,paws;
   double R,Rphi,Rchi,phi,chi;

//   if (particleObj[current_particle].x > LENGTH_MINUS_DX) sign = -1;
   if (particleObj[current_particle].x > LENGTH) sign = -1;

   R = (double)(rand())/((double)RAND_MAX);
   Rphi = (double)(rand())/((double)RAND_MAX);
   Rchi = (double)(rand())/((double)RAND_MAX);
   chi=asin(sqrt(Rchi));
   phi= 2.0*PI*Rphi;
   particleObj[current_particle].v = ((double)sign)*fabs(velocity*cos(chi));
   particleObj[current_particle].vy = velocity*sin(chi)*cos(phi);
   particleObj[current_particle].vz = velocity*sin(chi)*sin(phi);
   particleObj[current_particle].x = R*dt*fabs(particleObj[current_particle].v);

   if (sign < 0){
//  particleObj[current_particle].x = LENGTH_MINUS_DX-particleObj[current_particle].x;
         particleObj[current_particle].x = LENGTH-particleObj[current_particle].x;
   }

}

/*
 * Implementation Notes: secondary_electrons
 * ---------------------------------------------
 * Electron secondary velocity, 'velocity' corresponds to 'v' in:
 *   -Ralf Krimke et al, J. Phys. D: Appl. Phys.,vol 29,(1996),378
 */
void boundary::secondary_electrons(int current_particle,vector<particleVec> &particleObj){
   int last_particle;
   double R; 
   double velocity;

   particleObj.push_back(particleVec());
   last_particle  = particleObj.size()-1;

   particleObj[last_particle].weight = particleObj[current_particle].weight;
   particleObj[last_particle].x = particleObj[current_particle].x;

   R = (double)(rand())/((double)RAND_MAX);

   velocity= sqrt((2.0*R*CHARGE*(threshold-2.0*work_function)/ELECTRON_MASS));
   secondary_vel(last_particle,velocity,particleObj);
}

void boundary::absorb(vector<particleVec> &particleObj){
     int i=0;
     int n = particleObj.size();
     
     while( i < n){
         //if(particleObj[i].x < 0.0 || particleObj[i].x > LENGTH_MINUS_DX ){
         if(particleObj[i].x < 0.0 || particleObj[i].x >= LENGTH ){

              swap(particleObj[i],particleObj[(n-1)]);
              particleObj.pop_back();
              //i--;
         }
         else
           i++;
     n = particleObj.size();
     //i++;
     }
}

void boundary::e_boundary(vector<particleVec> &particleObj,vector<particleVec> &ionsObj){
   int sign=1,paws=1,no_inject=0;
   int i=0,j,count_ref=0,count_scat=0,count_sec=0,count_absorb=0;
   int n = particleObj.size();
   double R,velocity;
   double yield_max=0.0;
   double offset= 1.0e-6;

   while( i < n){

     //if(particleObj[i].x < 0.0 || particleObj[i].x > LENGTH_MINUS_DX ){
     if(particleObj[i].x < 0.0 || particleObj[i].x >= LENGTH ){

        calc_energy(particleObj[i].v,particleObj[i].vy,particleObj[i].vz);
        calc_angle(particleObj[i].v);

        (this->*e_impact_see_yield)();
        R = (double)(rand()%10000000)/1.0e7;

        if(this->yield > yield_max) yield_max=this->yield;

        if ( R < this->yield ){

              (this->*calc_reflection_yield)();
              (this->*calc_backscatter_yield)();
               no_inject = (int)floor(fabs((this->yield-this->reflection_yield-this->backscatter_yield)))+1;

        R = (double)(rand()%10000000)/1.0e7;

        if ( R < this->reflection_yield/(this->yield+offset)){
               reflect_electron(i,particleObj);
               count_ref++;
               i++;
         }

       else if (R > this->reflection_yield/this->yield && R < (this->reflection_yield + this->backscatter_yield)/(this->yield+offset) ){
               R = (double)(rand()%10000)/10000.0;
               velocity = R*sqrt((2.0*CHARGE*this->energy)/ELECTRON_MASS);
               secondary_vel(i,velocity,particleObj);
               count_scat++;
               i++;
        }

        else if (R >(this->reflection_yield+this->backscatter_yield)/(this->yield+offset)){
               for (j=0;j<no_inject; j++){
                    secondary_electrons(i,particleObj);
                    n = particleObj.size();
               }
               swap(particleObj[i],particleObj[(n-1)]);
               particleObj.pop_back();
               //i--;
               count_sec++;
         }
        } //if R < yield 
       
        else{
              /*absorb*/ 
              swap(particleObj[i],particleObj[(n-1)]);
              particleObj.pop_back();
              //i--;
              count_absorb++;
        }


     } // if particle outside boundary 
     else{
        i++;
     }

     n = particleObj.size();
     //i++;
 
     if(i>0 && (particleObj[i-1].x < 0.0 || particleObj[i-1].x > LENGTH )){
       cout<<"in while in electron boundary particle outside, particleObj.x= "<<particleObj[i-1].x<<" length= "<<LENGTH<<" i-1= "<<(i-1)<<" n= "<<particleObj.size()<<endl;
       cout<<"vel = "<<particleObj[i-1].v<<endl;
       cout<<"yield= "<<this->yield<<" R= "<<R<<endl;

    cout<<" yield max= "<<yield_max<<" no reflec= "<<count_ref<<" no scat= "<<count_scat<<" no sec= "<<count_sec<<" no absorb= "<<count_absorb<<endl;
      cin>>paws;
     }

   } //while loop

/*    cout<<" yield max= "<<yield_max<<" no reflec= "<<count_ref<<" no scat= "<<count_scat<<" no sec= "<<count_sec<<" no absorb= "<<count_absorb<<endl;
*/
}

/*
 * Implementation Notes: i_seconondary_vel
 * -----------------------------------------
 * Calculates the position and velocity of backscattered and true secondary
 * electrons using the formulae presented in:
 *   -Ralf Krimke et al, J. Phys. D: Appl. Phys.,vol 29,(1996),378
 *
 */
void boundary::i_secondary_vel(int last_particle,vector<particleVec> &electrons,double x, double weight){
   int no_inject=0,sign=1;
   double R;
   double velocity;
   double Rphi,Rchi,phi,chi;

   //if (x > LENGTH_MINUS_DX) sign=-1;
   if (x > LENGTH) sign=-1;

   R = (double)(rand())/((double)RAND_MAX);
   velocity= sqrt((2.0*R*CHARGE*(threshold-2.0*work_function)/ELECTRON_MASS));

   Rphi = (double)(rand())/((double)RAND_MAX);
   Rchi = (double)(rand())/((double)RAND_MAX);
   
   chi=asin(sqrt(Rchi));
   phi= 2.0*PI*Rphi;

   electrons[last_particle].v = ((double)sign)*fabs(velocity*cos(chi));
   electrons[last_particle].vy = velocity*sin(chi)*cos(phi);
   electrons[last_particle].vz = velocity*sin(chi)*sin(phi);
   electrons[last_particle].weight= weight;
   electrons[last_particle].x = R*dt*fabs(electrons[last_particle].v);

   if (sign < 0)
       electrons[last_particle].x = LENGTH -electrons[last_particle].x;
       //electrons[last_particle].x = LENGTH_MINUS_DX -electrons[last_particle].x;

}

void boundary::ion_boundary(vector<particleVec> &electrons,vector<particleVec> &ions){
   int no_inject=0,paws=1;
   int i=0,j,last_particle;
   int n = ions.size();
   double R;

   while( i < n){

     //if(ions[i].x < 0.0 || ions[i].x > LENGTH_MINUS_DX ){
     if(ions[i].x < 0.0 || ions[i].x >= LENGTH ){

         calc_energy(ions[i].v,ions[i].vy,ions[i].vz);

         (this->*i_impact_see_yield)();
        
         R = (double)(rand()%10000000)/1.0e7;
          
         if ( R < this->yield ){
             no_inject = (int)floor(this->yield)+1;
             for (j=0;j<no_inject; j++){
                electrons.push_back(particleVec());//verify adds to end
                last_particle  = electrons.size()-1;
                i_secondary_vel(last_particle,electrons,ions[i].x,ions[i].weight);
             }
          }  //if R<yield
                  
          swap(ions[i],ions[(n-1)]);
          ions.pop_back();
     //     i--;

      } //if x outside boundary 
      else{
         i++;
      }

      //i++;
      n=ions.size();
   
   } //while i<n

}

/*
 * Implementation Notes: get_particle_velocity
 * ---------------------------------------------
 * This method is a temporary place holder for a more accurate way to
 * calculate the particle velocity taken randomly off a particle velocity
 * distribution. In this case Maxwellian.
 */
double boundary::get_particle_velocity(double thermal_velocity)
{
    double R1, R2;
    double v = 0.0;

    do{
       R1 =(double)rand()/((double)RAND_MAX);
       R2 =(double)rand()/((double)RAND_MAX);
         v = thermal_velocity*cos(2.0*PI*R2)*sqrt(fabs(log(R1)));
    }while( fabs(v) > thermal_velocity*10000);

    return v;
}

void boundary::plop_particles(particleVec &elecElement, particleVec &ionElement,double vte, double vti){
     double random_number;

     random_number=(double)rand()/((double)RAND_MAX);

     elecElement.x = PLACE_AT; 
     ionElement.x = PLACE_AT; 

     elecElement.v = get_particle_velocity(vte);
     ionElement.v = get_particle_velocity(vti);

     elecElement.vy = get_particle_velocity(vte);
     ionElement.vy = get_particle_velocity(vti);

     elecElement.vz = get_particle_velocity(vte);
     ionElement.vz = get_particle_velocity(vti);
}

/*
 * Implementation Notes: steady_e_boundary
 * ---------------------------------------------
 * Implements a type of boundary that results in a rapid steady state, by
 * inserting an electron-ion pair in the simulation domain everytime an
 * electron arrives at the surface. Ionization should be turned off and
 * the ion time-step should equall the electron value. The particle number
 * is constant so the desired number of simulation particles should be given
 * the particle::init_particle_no parameter.
 */  
void boundary::steady_e_boundary(vector<particleVec> &electronObj,vector<particleVec> &ionObj){
     int i,j,k,n=electronObj.size();
     int n_ion=ionObj.size();
     int count_LHS=0,i_before;
     int count_RHS=0;
     int amount_to_plopLHS;
     int amount_to_plopRHS;
     int klast_lhs=0,klast_rhs=0;
     int no_plopped=0;


     for(i=0;i<n;i++){
        if(electronObj[i].x < 0.0 )
           count_LHS++;
        else if(electronObj[i].x >= LENGTH)
           count_RHS++;
     }

       amount_to_plopLHS = (int)((ion_proportion)*((double)count_LHS));
       amount_to_plopRHS = (int)((ion_proportion)*((double)count_RHS));
       i=0;
       while(i<n){
         if(electronObj[i].x< 0.0){
            k=klast_lhs;
            i_before=i;
            while(amount_to_plopLHS >0 && k<n_ion ){
                if(ionObj[k].x < 0.0)
                {           
                     plop_particles(electronObj[i],ionObj[k],VTE,VTI);
                     amount_to_plopLHS--;
                     i++;
                     no_plopped++;
                     break; //break out of while
                }
                k++;
            } //while k>
            klast_lhs = k+1;
            if(i==i_before){
                i++;
            }
         }//if electronObj.x < 0
         else if(electronObj[i].x>=LENGTH){
            k=klast_rhs;
            i_before=i;
            while(amount_to_plopRHS >0 && k<n_ion ){
               if(ionObj[k].x>= LENGTH)
               {
                     plop_particles(electronObj[i],ionObj[k],VTE,VTI);
                     amount_to_plopRHS--;
                     i++;
                     no_plopped++;
                     break; //break out of while
               }
               k++;
            } //while k> 
            klast_rhs = k+1;
            if(i==i_before){
                i++;
            }
        } //else if electronObj.x > Length
       else if(electronObj[i].x <= LENGTH && electronObj[i].x >= 0.0)
            i++;
       } //while i<n

     cout<<" no plopped= "<<no_plopped<<endl;
}

void boundary::steady_i_boundary(vector<particleVec> &electronObj,vector<particleVec> &ionObj){
    /* This is intentionally blank,everything is done in the electron boundary
     * method.
     */ 
}


