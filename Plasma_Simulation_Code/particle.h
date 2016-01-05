/*
 * File: particle.h
 * ------------------
 * This file exports the particle class which stores particle vectors and
 * associated functions that operate on particles, such as acceleration by
 * the grid electric field, and the collection of the total particle charge    
 * in each cell.
 *
 * 4/15
 *
 */

#ifndef particle_h
#define particle_h

#include "input_files/constants.h"
#include "input_files/parameters.h"
#include "particleVec.h"
#include <vector>

using namespace std;

class particle
{

   public:
           int count_reduce;
     /*
      * Constructor: particle
      * Usage: particle;
      *        particle(dt,vt,q_over_m,qc,mass,dxp);
      *
      */ 
      particle(double dt_p=DTI, double vt_p=VTI, double q_over_m_p=Q_OVER_MI, double qc_p=QCI,double mass=ION_MASS,double dxp=DXPI);

    /* Destructor: ~particle */
    ~particle();

    /*
     * Method: load_particles
     * Usage:  load_particles();
     * --------------------------
     * Starts the simulation by loading particle data dumped from a previous
     * run if the pic::restart is initialized to true via the pic constructor.
     */
      void load_particles(void); 

    /*
     * Method: randomize_vector
     * Usage:  randomize_vector();
     * ----------------------------
     * Shuffles the elements of a vector randomly.
     */
      void randomize_vector(void);

    /*
     * Method: dump_data
     * Usage:  dump_data();
     * ---------------------
     * Dumps particle data to a text file for each species at the end of a run
     * if pic::dump == true. Note, that this is set to true by default.
     */
      void dump_data(void);

    /*
     * Method: set_type
     * Usage:  set_type("i1");
     * -------------------------
     * Sets a string expression that is used to identify each species type
     * being dumped in a text file or read in for particle loading.
     */
      void set_type(char *type_in);


     /*
      * Method: get_init_no
      * Usage:  int init_particle_no = get_init_no();
      * ------------------------
      * Returns the initial number of particles to be loaded if not loading  
      * particles from a previous run.
      */
      int get_init_no(void);

     /*
      * Method: create_particles
      * Usage:  create_particles();
      * -----------------------------
      * Adds init_particle_no vector elements to the particle vectors before
      * the simulation starts, if pic::restart=false.
      */
     void create_particles(void);

     /*
      * Method: initiate_positions
      * Usage:  initiate_positions();
      * -------------------------------
      * Initiates particles positions along the x-axis before the simulations
      * starts if pic::restart=false.
      */
     void initiate_positions(void);

     /*
      * Method: initiate_velocities
      * Usage:  initiate_velocities();
      * --------------------------------
      * Initiates particle velocites from a Maxwellian distribution before the  
      * simulation starts if pic::restart=false.
      */
     void initiate_velocities(void);
    
     /*
      * Method: inject_more_particles
      * Usage:  inject_more_particles();
      * ----------------------------------
      * Inject no_inject number of particles each time-step randomly along the
      * x-axis.
      */  
     void inject_more_particles(void);

     /*
      * Method: contain_particle_growth
      * Usage:  contain_particle_growth()
      * -----------------------------------
      * Implements the methods that limit the size of the particle vectors, such
      * as initialize_map,and reduce_population.
      */
     void contain_particle_growth(void);

     /*
      * Method: initialize_map
      * Usage:  initialize_map();
      * ---------------------------
      * Initializes values used to reduce the particle population and increase
      * the particle wieghting. This function is called by 
      * contain_particle_growth.
      */
     void initialize_map(void);

     /*
      * Method: reduce_population
      * Usage:  reduce_population();
      * ------------------------------
      * Reduces the population of each velocity group down to the value
      * GROUP_SIZE, a constant determined by the user.
      */
     void reduce_population(void);

     /*
      * Method: group_v_ranges 
      * Usage:  group_v_ranges();
      * -----------------------------
      * Defines the lower and upper velocities of each velocity group using
      * the constants NO_GROUPS (the number of velocity groups) and 
      * VT_RANGE (the size of each velocity group, eg. 5*vt).
      */
     void group_v_ranges(void); 

     /*
      * Method: get_particle_velocity
      * Usage:  double v = get_particle_velocity();
      * --------------------------------------------
      * Returns a random particle velocity taken off a Maxwellian distribution.
      * This method is a temporary place holder for a more accurate one.
      */
     double get_particle_velocity(void);

     /*
      * Method: collect_charge
      * Usage:  collect_charge();
      * --------------------------
      * Calculates charge on the grid.    
      */
     void collect_charge(void);

     /*
      * Method: collect_extraCharge
      * Usage:  collect_extraCharge(old_ion_particle_number);
      * ------------------------------------------------------
      * Collects any extra ion charge in between ion time-steps if the ions
      * have a larger time-step, and the old ion population size is less than
      * the current one. 
      */
     void collect_extraCharge(int old_size);

     /*
      * Method: move_particles
      * Usage:  move_particle(grid_electric_field_array);
      * --------------------------------------------------
      * Calculates the new particle velocities and positions due to the
      * electric and optional magnetic fields.
      */
     void move_particles(double *Egrid); 

     /*
      * Method: is_equall_to
      * Usage:  is_equall_to(vector_to_be_changed,vector_to_copy,
      *                        first_vec_element_no,second_vec_element_no)
      * --------------------------------------------------------------------
      * Sets an element of the first vector equall to an element in the second.
      * This might not be a necessary method.
      */
     void is_equall_to(vector<particleVec> &obj1, vector<particleVec> &obj2,int i, int j);
    
    /* A set of methods that return private values for use in other classes,
     * such as the boundary, particle-collisions, and graph classes: 
     *
     *       vector<particleVec>& get_particle(void){return particleObj;}
     *       double get_q(int elem){return q[elem];}
     *       double* get_q(void){return q;}
     *       int particle_size(void);
     */

    /*
     * Method: get_particle 
     * Usage:  vector<particleVec> &particles =object.particle::get_particle();      * ------------------------------------------------------------------------
     * Returns a reference to the object's particle vector.    
     */ 
     vector<particleVec>& get_particle(void){return particleObj;}

    /*
     * Method: get_q
     * Usage:  double q_plot = object.particle::get_q(element_number);
     * -----------------------------------------------
     * Returns the particle charge on the grid at the specified cell number.
     */ 
     double get_q(int elem){return q[elem];}

    /*
     * Method: get_q
     * Usage:  double *q_plot = object.particle::get_q();
     * -----------------------------------
     * Returns a pointer to the particle charge on the grid array.
     */ 
     double* get_q(void){return q;}

    /*
     * Method: particle_size
     * Usage:  int number_particles = object.particle::particle_size();
     * -------------------------------------------------
     * Returns the length of the particle vector.
     */
     int particle_size(void){return particleObj.size();}

     /* A group of methods that sets the values of private variables:
      * 
      *     void set_vt(double vt);
      *     void set_q_over_m(double q_over_m);
      *     void set_mass(double mass);
      *     void set_dxp(double dxp);
      *     void set_init_no(int init_no);
      *     void set_dt(double dt);
      *     void set_qc(double qc);
      *     void set_no_inject(int no_to_inject);
      *     void set_no_groups(int val);
      *     void set_v_limit(double val,int element);
      *     void set_group_limit(int val,int element);
      *     void set_new_group_size(int val,int element);
      *     void set_size_to_check(int val);
      */

     /*
      * Method: set_vt
      * Usage: set_vt(thermal_velocity)
      * ---------------------------------    
      * Set the particle thermal velocity, sqrt( (kT/e)*e/m).
      */
      void set_vt(double vt);

     /*
      * Method: set_q_over_m
      * Usage:  set_q_over_m(mass_over_charge);
      * ----------------------------------------
      * Set the particle charge to mass ratio.                 
      */
      void set_q_over_m(double q_over_m);

     /*
      * Method: set_mass
      * Usage:  set_mass(particle_mass);
      * ---------------------------------
      * Set the particle mass.
      */
      void set_mass(double mass);

     /*
      * Method: set_dxp
      * Usage:  set_dxp(distance_between_particles);
      * ---------------------------------------------
      * Particles are loaded uniformly, at a distance dxp apart,
      * at the start of the simulation if pic::restart=false.
      */
      void set_dxp(double dxp);

     /*
      * Method: set_init_no
      * Usage:  set_init_no(initial_no_of_particles)
      * ---------------------------------------------
      * Set the initial number of particles to load into the simulation if
      * pic::restart = false.
      */ 
      void set_init_no(int init_no);

     /*
      * Method: set_dt
      * Usage:  set_dt(time_step);
      * ----------------------------
      * Set the time-step. The ion time-step can be optionally larger than
      * the electron.
      */
      void set_dt(double dt);

     /*
      * Method: set_qc
      * Usage:  set_qc(charge_of_each_macroParticle); 
      * ----------------------------------------------
      * Set the charge of each macro-particle. Typically this is equall to
      * the size of each macro-particle multiplied by the charge.
      */
      void set_qc(double qc);

     /*
      * Method: set_no_inject
      * Usage:  set_no_inject(no_particles_to_inject);
      * -----------------------------------------------
      * Set the number of particles to inject each time-step. Currently,
      * particles are inserted randomly along the x-axis. This variable is
      * currently mandatory. For no particles, set it to zero.
      */
      void set_no_inject(int no_to_inject);

     /*
      * Method: set_no_groups
      * Usage:  set_no_groups(no_groups);
      * -----------------------------------------------
      * Set the number of velocity groups to be used when decreasing
      * the particle number.
      */
      void set_no_groups(int val);

     /*
      * Method: set_v_limit
      * Usage:  set_v_limit(v_limit,group_element);
      * -----------------------------------------------
      * Set the upper velocity of a velocity group. Used to increase particle
      * weighting of particles that are first placed in user defined
      * velocity groups.
      */
      void set_v_limit(double val,int element);

     /*
      * Method: set_group_limit
      * Usage:  set_group_limit(group_limit,group_element);
      * -----------------------------------------------
      * Set the particle number limit of a velocity group. Used to increase 
      * particle weighting of particles that are first placed in user defined
      * velocity groups.
      */
      void set_group_limit(int val,int element);

     /*
      * Method: set_new_group_size
      * Usage:  set_new_group_size(group_limit,group_element);
      * -----------------------------------------------
      * Set the particle desired particle number of a velocity group. Used to  
      * increase particle weighting of particles that are first placed in user 
      * defined velocity groups.
      */
      void set_new_group_size(int val,int element);

     /*
      * Method: set_size_to_check
      * Usage:  set_size_to_check(val);
      * -----------------------------------------------
      * Putting particles into groups to check each group size is time
      * consuming so this is only done when the number of particles of a given
      * species reaches the value given by size_to_check.
      */  
      void set_size_to_check(int val);

    private:
  
      /*instance variables*/

     vector<particleVec> particleObj;

     int number_of_cells, particle_number,no_groups,size_to_check;
     double dx, dt, vt, q_over_m, qc,mass,dxp;

     double *q;
     double *v_limit;
     int *group_limit,*new_group_size;
     double *group_weight;
     int *j_group;
     int *Ngroup;
     int *Nremove;
     int init_particle_no;
     int no_inject;
     char *type;

};

#endif /*particle*/
