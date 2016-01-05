/*
 * File: fields.cpp
 * -----------------
 * This file implements the fields.h interface.
 *
 * 4/15
 *
 */

#include "fields.h"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>

fields::fields(void)
{
   int j;

   rho = new double[NO_OF_CELLS+1]();
   if (!rho){
      printf("Error allocating memory for rho\n");
      exit(1);
   }

   phi = new double[NO_OF_CELLS+1]();
   if (!phi){
      printf("Error allocating memory for phi\n");
      exit(1);
   }
   
   electric_field = new double[NO_OF_CELLS+1]();
   if (!electric_field){
      printf("Error allocating memory for electric_field\n");
      exit(1);
   }

  
   a = new double[NO_OF_CELLS+1]();
   if (!a){
      printf("Error allocating memory for a\n");
      exit(1);
   }

   b = new double[NO_OF_CELLS+1]();
   if (!b){
      printf("Error allocating memory for b\n");
      exit(1);
   }

   c = new double[NO_OF_CELLS+1]();
   if (!c){
      printf("Error allocating memory for c\n");
      exit(1);
   }

   gam = new double[NO_OF_CELLS + 1]();
   if (!gam){
      printf("Error allocating memory for gam\n");
      exit(1);
   }

  for(j=0;j<(NO_OF_CELLS+1);j++){
     electric_field[j] = 0.0;
     rho[j]=0.0;
     phi[j]=0.0;
     a[j]=0.0;
     b[j]=0.0;
     c[j]=0.0;
     gam[j]=0.0;
  }

  /* Initialize a,b,c for use in calculating the electrostatic potential
   * via a tri-diagonal matrix.
   */

   *(b+1) = -2.0;
   *(c+1) = 1.0;

   *(a + NO_OF_CELLS - 2) = 1.0;
   *(b + NO_OF_CELLS - 2) = -2.0;
   *(c + NO_OF_CELLS - 2) = 0.0;

   for(j=2; j< NO_OF_CELLS -2 ; j++){
      *(a+j) = 1.0;
      *(b+j) = -2.0;
      *(c+j) = 1.0;
   }
   
}

fields::~fields(){
   delete [] rho;
   delete [] phi;
   delete [] electric_field;
}

void fields::grid_rho(particle *ion, double *qe)
{
   int i,j;
   double total_ion_charge[NO_OF_CELLS]={0.0};


   for(j=0; j< NO_OF_CELLS; j++){
      for(i=0;i<NO_ION_SPECIES;i++){
         total_ion_charge[j]+=ion[i].particle::get_q(j);
       }
   }

   for(j=0; j< NO_OF_CELLS ; j++){
      *(rho+j) = (*(total_ion_charge+j) + qe[j] );
   }

   *(rho+NO_OF_CELLS) = 0.0;
} 

/* 
 * Implementation Notes: grid_potential
 * ---------------------------------------
 * Calculate the electrostatic potential at each grid point from the finite-
 * difference Poisson equation using the method outlined in:
 * C K Birdsall and A B Langdon, Plasma Physics via Computer Simulation,
 * Chapter 2 and Appendix D.
 *
 * The resulting tri-diagonal matrix is solved using the method outlined
 * in Numerical Recipes: Tridiagonal and Band-Diagonal Systems of Eqns.
 */
void fields::grid_potential(double phi0, double phiL)
{
   double bet;
   int j;
   double coef = (DX*DX/EPSILON);

   for(j=0; j< NO_OF_CELLS ; j++){
      rho[j] = -coef*rho[j];
   }

   rho[1] -= phi0;
   rho[NO_OF_CELLS -2] -= phiL;

   bet = *(b +1);
   *(phi+1) = *(rho+1)/bet;

   for(j=2; j< NO_OF_CELLS-1; j++){
      *(gam + j) = *(c + (j-1))/bet;
      bet = *(b+j) - *(a+j)*( *(gam+j));
      *(phi+j) = ( *(rho+j) - *(a+j)*( *(phi+(j-1)) ) )/bet;
   }

  *phi = phi0;
  *(phi + NO_OF_CELLS -1 ) = phiL;

  for(j = NO_OF_CELLS -3 ; j>= 1; j--)
     *(phi +j) = *(phi +j) - *(gam + (j+1))*( *(phi + (j+1)) );

  *(phi + NO_OF_CELLS) = 0.0;
}

/* 
 * Implementation Notes: grid_electric_field
 * ------------------------------------------
 * The electric field at each grid point is calculated using the formula in
 * C K Birdsall and A B Langdon, Plasma Physics via Computer Simulation,
 * Chapter 2.
 *
 */
void fields::grid_electric_field(void) 
{
   int j;
   double twoDX=2.0*DX;

   for(j=1; j < NO_OF_CELLS -1; j++)
      *(electric_field + j)=(*(phi+(j-1)) - *(phi+(j+1)) )/twoDX;


   *electric_field = 2.0*(*(electric_field + 1))-*(electric_field+2);
   *(electric_field + NO_OF_CELLS-1) = 2.0*(*(electric_field+NO_OF_CELLS -2)) -                                              *(electric_field + NO_OF_CELLS-3);
   *(electric_field + NO_OF_CELLS) = 0.0;

}
