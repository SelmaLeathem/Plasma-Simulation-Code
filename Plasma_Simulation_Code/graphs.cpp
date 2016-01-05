/*
 * File: graphs.cpp
 * ------------------
 * This file implements the graphs.h interface.
 *
 * 4/15
 *
 */

#include "graphs.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>

graphs::graphs(char *type_inp,int spread_in,double vt_in,double qc_in):type_in(type_inp),spread(spread_in),vt(vt_in),qc(qc_in)
{

   rho_plot = new double[NO_OF_CELLS]();
   if(!rho_plot){
     printf("Error allocating memory for rho_plot\n");
     exit(1);
   }

   q_plot = new double[NO_OF_CELLS]();
   if(!q_plot){
     printf("Error allocating memory for q_plot\n");
     exit(1);
   }

   phi_plot = new double[NO_OF_CELLS]();
   if(!phi_plot){
     printf("Error allocating memory for phi_plot\n");
     exit(1);
   }

   Egrid_plot = new double[NO_OF_CELLS]();
   if(!Egrid_plot){
     printf("Error allocating memory for Egrid_plot\n");
     exit(1);
   }

   d_plot1 = new double[NO_BINS+1]();
   if(!d_plot1){
     printf("Error allocating memory for d_plot1\n");
     exit(1);
   }

   d_plot2 = new double[NO_BINS+1]();
   if(!d_plot2){
     printf("Error allocating memory for d_plot2\n");
     exit(1);
   }

   temp_dist = new double[NO_BINS+1]();
   if(!temp_dist){
     printf("Error allocating memory for temp_dist\n");
     exit(1);
   }

   v_plot = new double[NO_OF_CELLS]();
   if(!v_plot){
     printf("Error allocating memory for v_plot\n");
     exit(1);
   }

   count_particles = new int[NO_OF_CELLS]();
   if(!count_particles){
     printf("Error allocating memory for count_particles\n");
     exit(1);
   }
}

graphs::~graphs(){

   delete [] rho_plot;
   delete [] q_plot;
   delete [] phi_plot;
   delete [] Egrid_plot;
   delete [] d_plot1;
   delete [] d_plot2;
   delete [] v_plot;
   delete [] count_particles;
   delete [] temp_dist;
}

void graphs::set_type_in(char* str){
     type_in=str;
} 

void graphs::set_spread(int val){
     spread=val;
} 

void graphs::set_vt(double vt_in){
     vt=vt_in;
}

void graphs::set_qc(double qc_in){
     qc=qc_in;
}

void graphs::set_dist_location(int location){
     dist_location=location;
}

void graphs::set_dist_location2(int location){
     dist_location2=location;
}

void graphs::sum_fields(double *rho,double *phi,double *Egrid)
{
   int j;
  
   for(j=0;j<NO_OF_CELLS;j++){
      *(rho_plot + j) += *(rho+j);
      *(phi_plot + j) += *(phi+j);
      *(Egrid_plot+j) += *(Egrid + j);
   }
}

void graphs::sum_density(double *q)
{
   int j;
  
   for(j=0;j<NO_OF_CELLS;j++){
      *(q_plot + j) += (q[j]/qc);
   }
}

void graphs::sum_av_velocities(vector<particleVec> &particleObj){
  int i,j;
  int particle_no;
  double v_grid_temp[UPPER_LIMIT] = {0.0};

  memset(count_particles,0,sizeof(*count_particles));
  
  particle_no = particleObj.size();

  for(i=0;i< particle_no; i++){

     j=(int)(round(particleObj[i].particleVec::x/DX));

     if(j>0 && j<NO_OF_CELLS-1){
        *(v_grid_temp+j) += particleObj[i].particleVec::v;
        *(count_particles + j) += 1;
     }
   }

  for(j=1; j< NO_OF_CELLS -1 ; j++){
     if(*(count_particles + j) != 0)
        *(v_plot + j) += *(v_grid_temp+j)/((double)(*(count_particles + j)));
  }
} 

void graphs::sum_av_dist1(vector<particleVec> &particleObj){
   int i,j,vbin,particle_no;

   particle_no = particleObj.size();

   for(i=0; i< NO_BINS; i++)
     *(temp_dist + i) = 0.0;

   for(i=0; i< particle_no; i++){
      j=(int)(round(particleObj[i].particleVec::x/DX));

      if(j>= (dist_location - spread) && j>0 && j<=(dist_location+spread) && j<NO_OF_CELLS){
        vbin =(int)((particleObj[i].particleVec::v/vt)*DIVIDE_BINS_BY); 
 
      if(abs(vbin) <= (NO_BINS-1)/2){
         *(temp_dist+vbin+(NO_BINS-1)/2)+=1.0;
      //   *(temp_dist+vbin+(NO_BINS-1)/2)+=particleObj[i].weight;
      }
     }
  }

   for(j=0; j< NO_BINS; j++)
     *(d_plot1+j) += *(temp_dist + j);
}

void graphs::sum_av_dist2(vector<particleVec> &particleObj){
   int i,j,vbin,particle_no;

   particle_no = particleObj.size();

   for(i=0; i< NO_BINS; i++)
     *(temp_dist + i) = 0.0;

   for(i=0; i< particle_no; i++){
      j=(int)(round(particleObj[i].particleVec::x/DX));

     if(j>= (dist_location - spread) && j>0 && j<=(dist_location+spread) && j<NO_OF_CELLS){

      vbin =(int)((particleObj[i].particleVec::v/vt)*DIVIDE_BINS_BY); 
 
      if(abs(vbin) <= (NO_BINS-1)/2){
         *(temp_dist+vbin+(NO_BINS-1)/2)+= 1.0;
       //  *(temp_dist+vbin+(NO_BINS-1)/2)+=particleObj[i].weight;
      }
     }
   }

   for(j=0; j< NO_BINS; j++)
     *(d_plot1+j) += *(temp_dist + j);
}


void graphs::average_fields(void){
  int j;

  for(j=0; j < NO_OF_CELLS; j++){
     *(rho_plot +j) = *(rho_plot+j)/((double)NO_GRAPH_STEPS);
     *(q_plot +j) = *(q_plot+j)/((double)NO_GRAPH_STEPS);
     *(phi_plot +j) = *(phi_plot+j)/((double)NO_GRAPH_STEPS);
     *(Egrid_plot+j) = *(Egrid_plot+j)/((double)NO_GRAPH_STEPS);
  }
}

void graphs::average_density(void){
  int j;

  for(j=0; j < NO_OF_CELLS; j++){
     *(q_plot +j) = *(q_plot+j)/((double)NO_GRAPH_STEPS);
  }
}

void graphs::average_velocities(void){
  int j;
  //double vte = return_vte();  const

  for(j=0;j < NO_OF_CELLS; j++)
     *(v_plot+j) = *(v_plot+j)/(VTE*((double)NO_GRAPH_STEPS));
}

void graphs::average_dist1(void){
  int j;

  for(j=0; j< NO_BINS; j++)
    *(d_plot1 + j) = *(d_plot1 + j)/((double)NO_GRAPH_STEPS);
}

void graphs::average_dist2(void){
  int j;

  for(j=0; j< NO_BINS; j++)
    *(d_plot1 + j) = *(d_plot1 + j)/((double)NO_GRAPH_STEPS);
}

/*
 * Implementations Notes: average_quantity
 * -----------------------------------------
 * Have not implemented this method yet, since requires several more methods
 * to return references to private variables, which upfront seems more
 * inefficient.
 */
/*
void graphs::average_quantity(int no_items,double factor,double *quantity){
      int j;

      for(j=0;j<no_items;j++)
         quantity[j]=quantity[j]/(factor*((double)NO_GRAPH_STEPS);
}
*/

void graphs::write_fields_to_file(void){
  int j;
  FILE *filePtr;

  filePtr = fopen("rho_phi_Egrid.txt", "w");
  
  for(j=0; j< NO_OF_CELLS; j++)
     fprintf(filePtr, "%e %e %e %e\n",(double)j,*(rho_plot+j),*(phi_plot +j), *(Egrid_plot+j));

  fclose(filePtr);
}

void graphs::write_density_to_file(void){
  int j;
  FILE *filePtr;

  char filename[80] = "density_"; 

  strcat(filename, type_in);
  strcat(filename, ".txt");

  filePtr = fopen(filename, "w");

  for(j=0; j< NO_OF_CELLS; j++)
    fprintf(filePtr,"%e %e \n",(double)j,*(q_plot+j));

  fclose(filePtr);
}

void graphs::write_v_to_file(void){
  int j;
  FILE *filePtr;

  char filename[80] = "x_velocity_";
  strcat(filename, type_in);
  strcat(filename, ".txt");

  filePtr = fopen(filename, "w");

  for(j=0; j< NO_OF_CELLS; j++)
    fprintf(filePtr,"%e %e \n",(double)j,*(v_plot+j));

  fclose(filePtr);
}

void graphs::write_dist1_to_file(void){
  int j;
  FILE *filePtr;

  char filename[80] = "first_v_dist_";
  strcat(filename, type_in);
  strcat(filename, ".txt");
  
  filePtr = fopen(filename, "w");
   
  for(j=0; j< NO_BINS; j++)
   fprintf(filePtr,"%e %e \n",(double)(j-(NO_BINS-1)/2)/((double)DIVIDE_BINS_BY), *(d_plot1+j));

  fclose(filePtr);
}

void graphs::write_dist2_to_file(void){
  int j;
  FILE *filePtr;

  char filename[80] = "second_v_dist_";
  strcat(filename, type_in);
  strcat(filename, ".txt");
  
  filePtr = fopen(filename, "w");
   
  for(j=0; j< NO_BINS; j++)
   fprintf(filePtr,"%e %e \n",(double)(j-(NO_BINS-1)/2)/((double)DIVIDE_BINS_BY), *(d_plot1+j));

  fclose(filePtr);
}

