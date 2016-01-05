/*
 * File: e_collisions.h
 * ---------------------
 * This file exports the e_collisions class which stores data and defines
 * methods that calculate a variety of electron impact collisions.
 *
 * 6/15
 */

#ifndef e_collision_h
#define e_collision_h

#include "collisions.h"

class e_collisions:public collisions
{

 public:

  /*
   * Constuctor: e_collisions
   * Usage: e_collisions e_coll;
   *        e_collisions e_coll(IS_ELECTRON,dt);
   */
   e_collisions(int if_electron_in=IS_ELECTRON,double dt_in=DT); 

  /* Destructor: ~e_collisions */
   ~e_collisions();


  /*
   * Method: calc_energy
   * Usage: calc_energy(x_velocity,y_velocity,z_velocity)
   * ------------------------------------------------------
   * Calculates the particle energy for use in cross-section and
   * collision calculations.
   */
   void calc_energy(double vx,double vy,double vz);

  /*
   * Method: collision
   * Usage: collision(electron_vector,ion_vector)
   * ---------------------------------------------
   * Determines if a collision will happen and if it does calls the method
   * to implement it.
   */
   void collision(vector<particleVec> &particleObj_e,vector<particleVec> &particleObj_i);

  /*
   * Method: ionization
   * Usage: ionization(element_number,electron_vector,ion_vector,
   *                                          electron_energy,coll_factor)
   * ----------------------------------------------------------------------
   * Implements electron-atom impact ionization collisions.          
   */
   void ionization(int i, vector<particleVec> &particleObj_e, vector<particleVec> &particleObj_i,double energy, double coll_factor);

  /*A group of methods that sets the values of private variables:
   *   
   *   void set_neutral_density(double density);
   *   void set_vtn(double neutral_velocity);
   *   void set_ion_mass(double mass);
   *   void set_constant_elastic(double val);
   *   void set_constant_ionization(double val);
   *   void set_crossSec_elastic(double (e_collisions::*func1)(void));
   *   void set_crossSec_ionization(double (e_collisions::*func1)(void));
   *   void set_ionization_threshold(double val);
   */

  /*
   * Method: set_neutral_density
   * Usage:  set_neutral_density(NEUTRAL_DENSITY);
   * ----------------------------------------------
   * Set the neutral density. This value is used to calculate the collision
   * probability.
   */
   void set_neutral_density(double density);

  /*
   * Method: set_vtn
   * Usage:  set_vtn(neutral_thermal_velocity);
   * -------------------------------------------
   * Set the neutral thermal velocity.
   */
   void set_vtn(double neutral_velocity);

  /*
   * Method: set_ion_mass
   * Usage:  set_ion_mass(ION_MASS);
   * --------------------------------
   * Set the ion_mass of the species the electrons are interacting with.
   */
   void set_ion_mass(double mass);
 
  /*
   * Method: set_constant_elastic
   * Usage:  set_constant_elastic(constant_cross_section_val);
   * -----------------------------------------------------------  
   * Sets the constant elastic cross-section value, if it is constant. This
   * variable does not need to be set. If the cross-section is not constant
   * do not set this value, otherwise, this must be set -in addition- to
   * setting the cross-section function pointer to point to the
   * sigma_const_elastic method defined in this file.
   * 
   */
   void set_constant_elastic(double val);

  /*
   * Method: set_constant_ionization
   * Usage:  set_constant_ionization(constant_cross_section_val);
   * -----------------------------------------------------------  
   * Sets the constant ionization cross-section value, if it is constant. This
   * variable does not need to be set. If the cross-section is not constant
   * do not set this value, otherwise, this must be set -in addition- to setting
   * the cross-section function pointer to point to the 
   * sigma_const_ionization method defined in this file.
   */
   void set_constant_ionization(double val);

  /*
   * Method: set_crossSec_elastic
   * Usage:  set_crossSec_elastic(&e_collisions::elastic_crossSec_function);
   * ------------------------------------------------------------------------
   * Sets the 'double (e_collisions::*sigma_elastic)(void)' function pointer to 
   * point to the function passed in the arguement. User defined functions
   * can be passed. This function returns energy-dependent elastic
   * cross sections.
   */
   void set_crossSec_elastic(double (e_collisions::*func1)(void));

  /*
   * Method: set_crossSec_ionization
   * Usage: set_crossSec_ionization(&e_collisions::ionization_crossSec_func);
   * ------------------------------------------------------------------------
   * Sets the 'double (e_collisions::*sigma_ionization)(void)' function pointer 
   * to point to the function passed in the arguement. User defined functions
   * can be passed. This function returns energy-dependent ionization
   * cross sections.
   */
   void set_crossSec_ionization(double (e_collisions::*func1)(void));

  /*
   * Method: set_ionization_threshold
   * Usage:  set_ionization_threshold(threshold_value);
   * ----------------------------------------------------  
   * Sets the ionization threshold in eV.
   */
   void set_ionization_threshold(double val);

  /*
   * Method: sigma_const_elastic
   * Usage: object.set_crossSec_elastic(&e_collisions::sigma_const_elastic); 
   * -------------------------------------------------------------------------
   * Passed as an arguement to the set_crossSec_elastic method to implement
   * the constant elastic cross-section defined using the
   * set_constant_elastic method. Note: set_constant_elastic(value) must
   * be called first. 
   */
   double sigma_const_elastic(void);

  /*
   * Method: sigma_const_ionization
   * Usage: object.set_crossSec_ionization(
                                       &e_collisions::sigma_const_ionization); 
   * -------------------------------------------------------------------------
   * Passed as an arguement to the set_crossSec_ionization method to implement
   * the constant ionization cross-section defined using the
   * set_constant_ionization method. Note: the 
   * set_constant_ionization(value) method must be called first.
   */
   double sigma_const_ionization(void);

   /* A set of sample methods that return energy dependent cross-sections
    * to be passed as an arguement to the set_crossSec_elastic or
    * set_crossSec_ionization methods. In general the user is expected to 
    * supply functions.
    */
     double sigma_Ar_elastic(void);
     double sigma_Ar_ionize(void);

   /* User defined cross-section method declarations are defined in
    * "coll_input_declare"
    */
   #include "input_files/e_coll_input_declare"

   private:
     double energy;
     double threshold;
     double constant_elastic;
     double constant_ionization;
    
   /*
    * Method: (e_collisions::*sigma_elastic)(void)
    * Usage:  double sigma_elastic = (this->*sigma_elastic)();
    * ---------------------------------------------------
    * Points to a function selected by the user. This function returns
    * the elastic cross-section.
    */
     double (e_collisions::*sigma_elastic)(void);

   /*
    * Method: (e_collisions::*sigma_ionization)(void)
    * Usage:  double sigma_ionization = (this->*sigma_ionization)();
    * ----------------------------------------------------------------
    * Points to a function selected by the user. This function returns
    * the ionization cross-section.
    */
     double (e_collisions::*sigma_ionization)(void);
};

#endif
