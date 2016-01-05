/*
 * File: pic.cpp
 * ---------------
 * This file implements the pic.h interface.
 * Eventually, will port to a system where the user defines variables and
 * functions in a single file. In the meantime, all parameters are defined
 * in the parameters.h and input files. User defined collision cross-section
 * methods are declared in species_coll_input_declare and are implemented 
 * in species_coll_input_methods. User defined boundary yield functions are 
 * declared in boundary_input_declare and implemented in boundary_input_methods.
 *
 * 10/15
 *
 */

#include "pic.h"

#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <fstream>
#include <string.h>


/*
 * Constructor: pic
 * ------------------    
 * Initializes a particle-in-cell simulation as described in:
 *
 * -C K Birdsall and A B Langdon, Plasma Physics Via Computer Simulation
 * and the software on www.eecs.berkeley.edu
 *
 * Note: The variable restart can only be initialized via the constructor, and
 * is false by default.
 *
 */
pic::pic(bool restart_p,bool dump_p,int run_time_p,int start_graph_p,int end_graph_p):restart(restart_p),dump(dump_p),run_time(run_time_p),start_graph(start_graph_p),end_graph(end_graph_p){

   int i,ion_particle_no=0;

   ions = new particle[NO_ION_SPECIES]();
   if (!ions){
      printf("Error allocating memory for ions\n");
      exit(1);
   }

   e_collision = new e_collisions[NO_ION_SPECIES]();
   if (!e_collision){
      printf("Error allocating memory for e_collision\n");
      exit(1);
   }

   ion_collision = new ion_collisions[NO_ION_SPECIES]();
   if (!ion_collision){
      printf("Error allocating memory for ion_collision\n");
      exit(1);
   }

   sec_ions = new boundary[NO_ION_SPECIES]();
   if (!sec_ions){
      printf("Error allocating memory for sec_ions\n");
      exit(1);
   }
   
   graphs_ions = new graphs[NO_ION_SPECIES]();
   if (!graphs_ions){
      printf("Error allocating memory for graphs_ions\n");
      exit(1);
   }
   
  /* Parameter definitions */     
     #include "input_files/input"    

/* load particles into simulation region */
   if(restart==true){
      electrons.particle::load_particles();
      for(i=0;i<NO_ION_SPECIES;i++)
         ions[i].particle::load_particles();
   }
   else{
      for(i=0;i< NO_ION_SPECIES;i++){
         ions[i].particle::create_particles();
         ions[i].particle::initiate_positions();
         ions[i].particle::initiate_velocities();
      }
   electrons.particle::create_particles();
   electrons.particle::initiate_positions();
   electrons.particle::initiate_velocities();
   }


/* initialize neutral velocity arrays */
   for(i=0;i<NO_ION_SPECIES;i++){
      ion_collision[i].ion_collisions::calc_neutral_velocity();
      e_collision[i].e_collisions::calc_neutral_velocity();
   }


   for(i=0;i<NO_ION_SPECIES;i++){
      old_size[i] = ions[i].particle::particle_size();
      N_before[i] = old_size[i];
      N_after[i] = old_size[i];
   }

/* Variables used to store data variations in time */

   timedata.open("timeData.txt");
   if(!timedata){
      cout<<"Cannot open file"<<endl;
      exit(1);
   }

   time_plot = new int[SIZE_TIME_DATA]();
   if (!time_plot){
       printf("Error allocating memory for time_plot\n");
       exit(1);
   }

   no_e = new int[SIZE_TIME_DATA]();
   if (!no_e){
       printf("Error allocating memory for no_e\n");
       exit(1);
   }

   no_i = new int*[NO_ION_SPECIES]();
   if (!no_i){
       printf("Error allocating memory for no_i\n");
       exit(1);
   }
   for(i=0;i<NO_ION_SPECIES;i++){
        no_i[i]= new int[SIZE_TIME_DATA];
        if (!no_i[i]){
             printf("Error allocating memory for no_i[i]\n");
             exit(1);
        }
   }

   phi_left = new double[SIZE_TIME_DATA]();
   if (!phi_left){
       printf("Error allocating memory for phi_left\n");
       exit(1);
   }

   phi_mid = new double[SIZE_TIME_DATA]();
   if (!phi_mid){
       printf("Error allocating memory for phi_mid\n");
       exit(1);
   }

  phi_right = new double[SIZE_TIME_DATA]();
  if (!phi_right){
       printf("Error allocating memory for phi_right\n");
       exit(1);
   }

   eden_right = new double[SIZE_TIME_DATA]();
   if (!eden_right){
       printf("Error allocating memory for eden_right\n");
       exit(1);
   }

   iden_right = new double*[NO_ION_SPECIES]();
   if (!iden_right){
       printf("Error allocating memory for iden_right\n");
       exit(1);
   }
   for(i=0;i<NO_ION_SPECIES;i++){
       iden_right[i]= new double[SIZE_TIME_DATA];
       if (!iden_right[i]){
           printf("Error allocating memory for iden_right[i]\n");
           exit(1);
       }
   }

}

void pic::set_dump(bool dump_in){
     dump=dump_in;
}

void pic::set_run_time(int run_time_in){
     run_time=run_time_in;
}

void pic::set_start_graph(int start_graph_in){
     start_graph=start_graph_in;
}

void pic::set_end_graph(int end_graph_in){
     end_graph=end_graph_in;
}

/*
 * Implementation Notes: pic_run
 * -------------------------------
 * Implements the simulation cycle as described in Birdsall and Langdon,
 * C K Birdsall and A B Langdon,Plasma Physics via Computer Simulation, 
 * Chapter 2 and
 * C K Birdsall, IEEE Trans. Plasma Sci.,vol 19, (1991) 65:
 *
 * Integration of equations of motion, moving particles -> Collisions 
 *   ^                                                        |
 *   |                                                        v
 * Wieghting <- Integration of field equations on grid <- Weighting
 * 
 */
void pic::pic_run(void){
  int i,j,time_step,plot_elem=0;

 /* Catch any particles that might accidently get loaded outside of the
  * boundaries due to a poor choice of starting values.
  */  
  //sec_e.boundary::absorb(electrons.particle::get_particle());

  for(i=0;i<NO_ION_SPECIES;i++){
     //sec_ions[i].boundary::absorb(ions[i].particle::get_particle());
    //sec_ions[i].boundary::ion_boundary(electrons.particle::get_particle(),ions[i].particle::get_particle());
    //sec_e.boundary::e_boundary(electrons.particle::get_particle(),ions[i].particle::get_particle());
    
             (sec_ions[i].*sec_ions[i].boundary::implement_i_boundary)(
                                       electrons.particle::get_particle(),
                                       ions[i].particle::get_particle());
             (sec_e.*sec_e.boundary::implement_e_boundary)(
                                       electrons.particle::get_particle(),
                                       ions[i].particle::get_particle());
  }

  for(time_step=0; time_step<run_time; time_step++){

 //cout<<"time is "<<time_step<<endl;

 //Use to randomly reorder particle vector before call collisions which
 //only allows a limited number of ionization events       
   electrons.particle::randomize_vector();

     for(i=0;i<NO_ION_SPECIES;i++)
        e_collision[i].e_collisions::collision(electrons.particle::get_particle(), ions[i].particle::get_particle());

  //Optional particle source during simulation.

     electrons.particle::inject_more_particles();
     for(i=0;i<NO_ION_SPECIES;i++){
        ions[i].particle::inject_more_particles();
     }
 
     electrons.particle::move_particles(grid_fields.fields::get_elec_field());

     //sec_e.boundary::absorb(electrons.particle::get_particle());
    //sec_e.boundary::e_boundary(electrons.particle::get_particle(),ions[i].particle::get_particle());

 //***************************move ions******************************//

     if(time_step%DT_DTI == 0){

         for(i=0;i<NO_ION_SPECIES;i++){
             ion_collision[i].ion_collisions::collision(ions[i].particle::get_particle());

             ions[i].particle::move_particles(grid_fields.fields::get_elec_field());

             //sec_ions[i].boundary::absorb(ions[i].particle::get_particle());
             //sec_ions[i].boundary::ion_boundary(electrons.particle::get_particle(),ions[i].particle::get_particle());
           
              (sec_ions[i].*sec_ions[i].boundary::implement_i_boundary)(
                                       electrons.particle::get_particle(),
                                       ions[i].particle::get_particle());
              (sec_e.*sec_e.boundary::implement_e_boundary)(
                                       electrons.particle::get_particle(),
                                       ions[i].particle::get_particle());

             ions[i].particle::collect_charge();
             old_size[i] = ions[i].particle::particle_size();
         }
     }  //if ion_timestep
     else
         for(i=0;i<NO_ION_SPECIES;i++){
             (sec_e.*sec_e.boundary::implement_e_boundary)(
                                       electrons.particle::get_particle(),
                                       ions[i].particle::get_particle());
  }
 //*******************************************************************//

     electrons.particle::contain_particle_growth();

     for(i=0;i<NO_ION_SPECIES;i++){
          N_before[i] = ions[i].particle::particle_size();
          ions[i].particle::contain_particle_growth();
          N_after[i] = ions[i].particle::particle_size();

          if(N_after[i] < N_before[i]){
             old_size[i] = N_after[i];
             ions[i].particle::collect_charge();
          }

          if(N_after[i] > old_size[i]){
             ions[i].particle::collect_extraCharge(old_size[i]);
             old_size[i] = ions[i].particle::particle_size();
          }
     }

     electrons.particle::collect_charge();

     grid_fields.fields::grid_rho(ions,electrons.particle::get_q());
     grid_fields.fields::grid_potential(PHI_LEFT,PHI_RIGHT);
     grid_fields.fields::grid_electric_field();

   //data vs time

    if (time_step%TIME_DATA_FREQUENCY == 0){
        time_plot[plot_elem]=time_step;
        no_e[plot_elem]=electrons.particle::particle_size();
        phi_left[plot_elem] = grid_fields.get_phi(POINT1);
        phi_mid[plot_elem] = grid_fields.get_phi(POINT2);
        phi_right[plot_elem] = grid_fields.get_phi(POINT3);
        eden_right[plot_elem] = electrons.get_q(POINT3);

        for(i=0;i<NO_ION_SPECIES;i++){
            no_i[i][plot_elem]=ions[i].particle::particle_size();
            iden_right[i][plot_elem] = ions[i].get_q(POINT3);
        }
        plot_elem++;
    }

/* Print particle number to screen during simulation */
    cout<<"time is "<<time_step<<" no of electrons= "<<electrons.particle::particle_size();
    for(i=0;i<NO_ION_SPECIES;i++)
       cout<<" no ions of species"<<i+1<<"  "<<ions[i].particle::particle_size();
    cout<<endl;


    if(time_step >= start_graph && time_step< end_graph){

        graphs_fields.graphs::sum_fields(grid_fields.fields::get_rho(),grid_fields.fields::get_phi(),grid_fields.fields::get_elec_field());

        graphs_electrons.graphs::sum_density(electrons.particle::get_q());
        graphs_electrons.graphs::sum_av_velocities(electrons.particle::get_particle());
        graphs_electrons.graphs::sum_av_dist1(electrons.particle::get_particle());
        graphs_electrons.graphs::sum_av_dist2(electrons.particle::get_particle());

        for(i=0;i<NO_ION_SPECIES;i++){
           graphs_ions[i].graphs::sum_density(ions[i].particle::get_q());
           graphs_ions[i].graphs::sum_av_velocities(ions[i].particle::get_particle());
           graphs_ions[i].graphs::sum_av_dist1(ions[i].particle::get_particle());
           graphs_ions[i].graphs::sum_av_dist2(ions[i].particle::get_particle());
        }
    }

    if(time_step == end_graph){

        graphs_fields.graphs::average_fields();

        graphs_electrons.graphs::average_velocities();
        graphs_electrons.graphs::average_density();
        graphs_electrons.graphs::average_dist1();
        graphs_electrons.graphs::average_dist2();

        graphs_fields.graphs::write_fields_to_file();

        graphs_electrons.graphs::write_density_to_file();
        graphs_electrons.graphs::write_v_to_file();
        graphs_electrons.graphs::write_dist1_to_file();
        graphs_electrons.graphs::write_dist2_to_file();

        for(i=0;i<NO_ION_SPECIES;i++){
           graphs_ions[i].graphs::average_velocities();
           graphs_ions[i].graphs::average_density();
           graphs_ions[i].graphs::average_dist1();
           graphs_ions[i].graphs::average_dist2();

           graphs_ions[i].graphs::write_density_to_file();
           graphs_ions[i].graphs::write_v_to_file();
           graphs_ions[i].graphs::write_dist1_to_file();
           graphs_ions[i].graphs::write_dist2_to_file();
       }
    }


} /*timestep loop*/

  //put time data in file 

   for(j=0;j<SIZE_TIME_DATA ;j++){

      timedata<<time_plot[j]<<" "<<no_e[j];
      for(i=0;i<NO_ION_SPECIES;i++)
          timedata<<" "<<no_i[i][j];
      timedata<<" "<<phi_left[j]<<" "<<phi_mid[j]<<" "<<phi_right[j]<<" "<<eden_right[j];
      for(i=0;i<NO_ION_SPECIES;i++)
          timedata<<" "<<iden_right[i][j];
      timedata<<endl;
   }

 //dump particle data to file

   if(dump==true){
      electrons.particle::dump_data();
      for(i=0;i<NO_ION_SPECIES;i++)
          ions[i].particle::dump_data();
   }
   
   cout<<"no electron reductions= "<<electrons.particle::count_reduce<<endl;
   cout<<"no ion reductions= "<<ions[0].particle::count_reduce<<endl;
}

pic::~pic(){
    int i;
    delete [] ions;
    delete [] e_collision;
    delete [] ion_collision;
    delete [] graphs_ions;
    delete [] time_plot;
    delete [] no_e;
    delete [] phi_mid;
    delete [] phi_right;
    delete [] phi_left;
    delete [] eden_right;
    for(i=NO_ION_SPECIES;i>0;i--){
        delete[] no_i[i-1];
        delete[] iden_right[i-1];
    }
    delete[] no_i;
    delete[] iden_right;
    delete[] sec_ions;
}

