/*
 * File: pic.h
 * ------------
 * This file exports the pic(particle-in-cell code) class. The pic class 
 * initializes and runs a particle-in-cell simulation by declaring objects
 * and calling the methods defined in the boundary, particle, fields,  
 * collisions, e_collisions, ion_collisions, and graphs classes.
 *
 * 10/15
 *
 */


#ifndef pic_h
#define pic_h

#include "input_files/constants.h"
#include "input_files/parameters.h"
#include "boundary.h"
#include "particle.h"
#include "fields.h"
#include "collisions.h"
#include "e_collisions.h"
#include "ion_collisions.h"
#include "graphs.h"

using namespace std;

class pic{

   public:
  
   /*
    * Constructor: pic
    * Usage: pic pic_object;
    *        pic pic_object(true);
    *        pic pic_object(true,true,100000,90000,99999);
    * ------------------------------------------------------  
    */ 
    pic(bool restart=false,bool dump=true,int run_time=TOTAL_RUN_TIME,int start_graph=START_GRAPH_TIME,int end_graph=END_GRAPH_TIME);

   /* Destructor: ~pic */
   ~pic();

   /*
    * Method: pic_run
    * Usage:  pic_run();
    * -------------------
    * Implements the simulation cycle over the number of iterations defined by
    * the run_time variable.
    */
   void pic_run(void);

   /* The following methods set the values of private variables:
    *
    *   void set_dump(bool dump_in);
    *   void set_run_time(int run_time_in);
    *   void set_start_graph(int start_graph_in);
    *   void set_end_graph(int end_graph_in);
    *   void set_time_point1(int time_point);
    *   void set_time_point2(int time_point);
    *   void set_time_point3(int time_point);
    *
    */

   /*
    * Method: set_dump
    * Usage:  set_dump(false);
    * --------------------------
    * Pass false to turn off dumping simulation data into text files at the 
    * end of the simulation. The default is true. These files can be used
    * to start a new simulation where the previous one left off.
    */
   void set_dump(bool dump_in);

   /*
    * Method: set_run_time
    * Usage:  set_run_time(number_of_iterations); 
    * --------------------------------------------
    * An alternative way to set the number of iterations that the simulation 
    * steps through. 
    */
   void set_run_time(int run_time_in);

   /*
    * Method: set_start_graph
    * Usage:  set_start_graph(time_to_start_summing_av_data);
    * --------------------------------------------------------
    * An alternative way to set the iteration number from which to start
    * summing average data for the graphs class.
    */
   void set_start_graph(int start_graph_in);

   /*
    * Method: set_end_graph
    * Usage:  set_end_graph(time_to_stop_summing_av_data);
    * -----------------------------------------------------
    * An alternative way to set the iteration number to stop summing
    * average data for the graphs class.
    */
   void set_end_graph(int end_graph_in);

   /*
    * Method: set_time_point1, set_time_point2, set_time_point3
    * Usage:  set_time_point1(cell_position1),set_time_point2(cell_position2),
    *         set_time_point3(cell_position3)
    * -------------------------------------------------------------------------
    * The three cells at which to dump the variation in time of several
    * variables including the potential and charge density.
    */
   void set_time_point1(int time_point);
   void set_time_point2(int time_point);
   void set_time_point3(int time_point);

   private:

   particle electrons,*ions;
   fields grid_fields;
   e_collisions *e_collision;
   ion_collisions *ion_collision;
   graphs graphs_fields,graphs_electrons,*graphs_ions;

   int run_time,start_graph,end_graph;
   int no_ions_inject[NO_ION_SPECIES];
   int old_size[NO_ION_SPECIES],N_before[NO_ION_SPECIES],N_after[NO_ION_SPECIES];
  /*
   * Variables used to store the variation of select data in time, including
   * the total particle_number for each species, electric_potential at three
   * poistions, and the particle charge density of each species at position3.
   */
   int *time_plot,*no_e,**no_i;
   double *phi_mid,*phi_left,*phi_right,*eden_right,**iden_right;
  
   ofstream timedata; 

   bool dump,restart;
   boundary sec_e, *sec_ions; /* might need a sec_electron for each atomic 
                                 species */

};

#endif /*pic*/
