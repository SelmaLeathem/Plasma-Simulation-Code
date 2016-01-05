/*
 * File: graphs.h
 * ---------------
 * This file exports the graphs class which stores variables and methods
 * that make text files of diagnostic data including the average fields, 
 * particle velocites, densities, and particle velocity distributions over
 * a period of time selected by the user.
 *
 * 4/15
 */

#include "input_files/constants.h"
#include "input_files/parameters.h"
#include "particleVec.h" 
#include <vector>

using namespace std;

class graphs{
 
  public:

   /*
    * Constructor: graphs
    * Usage: graphs graph_object;
    *        graphs graph_object(ION_TYPE,SPREAD,VT,QCI);
    */       
    graphs(char *type_inp=ION_TYPE,int spread_in=SPREAD,double vt_in=VTI,double qc_in=QCI);
    
   /* Destructor: ~graphs */
    ~graphs();

   /* These methods sum the indicated quantities over a number of time-steps
    * defined by the user:
    *
    * void sum_fields(double *rho,double *phi,double *Egrid);
    * void sum_density(double *q);
    * void sum_av_velocities(vector<particleVec> &particleObj);
    * void sum_av_dist1(vector<particleVec> &particleObj);
    * void sum_av_dist2(vector<particleVec> &particleObj);
    */

   /*
    * Method: sum_fields
    * Usage:  sum_fields(net_charge_density,electric_potential,electric_field);
    * -------------------------------------------------------------------------
    * Sums field quantities, including the net_charge density, electric 
    * potential and the electric field.
    */
    void sum_fields(double *rho,double *phi,double *Egrid);

   /*
    * Method: sum_density
    * Usage:  sume_density(particle_charge_density);
    * -----------------------------------------------
    * Sums the charge density of each species.
    */
    void sum_density(double *q);

   /* 
    * Method: sum_av_velocities
    * Usage:  sum_av_velocities(particle_vector);
    * --------------------------------------------
    * Sums the average velocity of each grid cell for the species passed to 
    * the arguement.
    */
    void sum_av_velocities(vector<particleVec> &particleObj);

   /*
    * Methods: sum_av_dist1, sum_av_dist2
    * Usage:   sum_av_dist1(particle_vector), sum_av_dist2(particle_vector);
    * -----------------------------------------------------------------------
    * Sums the particle velocity distribution at cell locations defined by the
    * dist_location and dist_location2 variables. Particle velocities are 
    * included from several cells around each location indicated by the 
    * spread variable.
    */ 
    void sum_av_dist1(vector<particleVec> &particleObj);
    void sum_av_dist2(vector<particleVec> &particleObj);

   /*
    * The following functions take averages of the indicated quantity by
    * dividing it over the number of steps each item was summed up over.
    * 
    * void average_fields(void);
    * void average_density(void);
    * void average_velocities(void);
    * void average_dist1(void);
    * void average_dist2(void);
    *
    */

    void average_fields(void);
    void average_density(void);
    void average_velocities(void);
    void average_dist1(void);
    void average_dist2(void);

   /*
    * No general average function has been implemented at this time because
    * would need to return reference to private graph variables requiring
    * several more get methods.
    * 
    * void average_quantity(int no_items,double factor,double *quantity);
    *
    */

   /* 
    * The following methods write the average of the indicated quantity to
    * a text file.
    *
    * void write_fields_to_file(void);
    * void write_density_to_file(void);
    * void write_v_to_file(void);
    * void write_dist1_to_file(void);
    * void write_dist2_to_file(void);
    *
    */

    void write_fields_to_file(void);
    void write_density_to_file(void);
    void write_v_to_file(void);
    void write_dist1_to_file(void);
    void write_dist2_to_file(void);

   /*
    * A group of methods that sets the values of private variables.
    *
    * void set_dist_location(int location);
    * void set_dist_location2(int location);
    * void set_type_in(char* str);
    * void set_spread(int val);
    * void set_vt(double vt_in);
    * void set_qc(double qc_in);
    *
    */

   /* 
    * Method: set_dist_location
    * Usage:  set_dist_location(cell_number);
    * ----------------------------------------   
    * Pick the central cell from which particle velocities are picked to
    * determine the particle velocity distribution.
    */
    void set_dist_location(int location);

   /* 
    * Method: set_dist_location2
    * Usage:  set_dist_location2(cell_number);
    * ----------------------------------------   
    * Pick a second central cell from which particle velocities are picked to
    * determine the particle velocity distribution.
    */
    void set_dist_location2(int location);

   /*
    * Method: set_spread
    * Usage:  set_spread(number_of_cells);
    * -------------------------------------
    * Set the number of cells around the central one from which to pick
    * particle velocities from. These are used to calculate the particle
    * velocity distribution.
    */
    void set_spread(int val);
 
   /*
    * Method: set_type_in
    * Usage:  set_type_in("i1");
    * ----------------------------
    * Define a string that is used to distinguish which species each diagnostic
    * text file is reffering to.
    */
    void set_type_in(char* str);
   
   /*
    * Method: set_vt
    * Usage:  set_vt(thermal_velocity);
    * ----------------------------------
    * Set the value of the thermal velocity used to normalize the particle
    * velocity distribution.
    */
    void set_vt(double vt_in);

   /*
    * Method: set_qc
    * Usage:  set_qc(charge_of_macro_particle);
    * -------------------------------------------
    * Use to set the the value of the charge of each macro-particle that is
    * multiplied by the particle number density of each cell. To plot the number
    * density set this value to 1.0.
    */
    void set_qc(double qc_in);

    private:

    int spread;
    int dist_location,dist_location2;

    double *rho_plot;
    double *q_plot;
    double *phi_plot;  
    double *Egrid_plot;
    double *d_plot1;
    double *d_plot2;
    double *v_plot;
    
    int *count_particles; 
    double *temp_dist;
    char *type_in;

   /* cell size */
    double dx;
    double vt;
    double qc;

};
