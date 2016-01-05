/*
 * File: ion_collisions.h
 * ------------------------
 * This file exports the ion_collisions class which stores data and defines
 * methods that calculate a variety of ion collisions
 *
 * 6/15
 *
 */

#ifndef ion_collision_h
#define ion_collision_h

#include "collisions.h"

class ion_collisions: public collisions
{

 public:

  /*
   * Constuctor: ion_collisions
   * Usage: ion_collisions ion_coll;
   *        ion_collisions ion_coll(IS_ION,dt);
   */
   ion_collisions(int if_electron_in=IS_ION, double dt_in=DTI, bool CM_elastic=false, bool CM_CX=false);

  /* Destructor: ~ion_collisions */
   ~ion_collisions();

  /*
   * Method: calc_energy
   * Usage: calc_energy(x_velocity,y_velocity,z_velocity)
   * ------------------------------------------------------
   * Calculates the particle energy for use in cross-section and
   * collision calculations.
   */
   void calc_energy(double vx,double vy,double vz);

  /*
   * Method: calc_CM_energy
   * Usage: calc_CM_energy(x_velocity,y_velocity,z_velocity)
   * ------------------------------------------------------
   * Calculates the particle center-of-mass energy for use in cross-section
   * calculations.
   */
   void calc_CM_energy(double vx,double vy,double vz);

  /*
   * Method: collision
   * Usage: collision(ion_vector)
   * ---------------------------------------------
   * Determines if a collision will happen and if it does calls the method
   * to implement it.
   */
   void collision(vector<particleVec> &particleObj);

  /*
   * Method: charge_exchange_velocity
   * Usage:  charge_exchange_velocity( element_number,ion_vector);
   * --------------------------------------------------------------
   * Implements charge-exchange collisions.
   */
   void charge_exchange_velocity(int i,vector<particleVec> &particleObj);

  /* A group of methods that sets the values of private variables:
   *
   *   void set_neutral_density(double density);
   *   void set_vtn(double neutral_velocity);
   *   void set_ion_mass(double mass);
   *   void set_CMtrue_elastic(bool val);
   *   void set_CMtrue_CX(bool val);
   *   void set_constant_elastic(double val);
   *   void set_constant_CX(double val);
   *   void set_crossSec_elastic(double (ion_collisions::*func1)(void));
   *   void set_crossSec_CX(double (ion_collisions::*func1)(void));
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
   * Method: set_CMtrue_elastic
   * Usage:  set_CMtrue_elastic(true);
   * ----------------------------------
   * Set center_mass_El to true if the ion-atom center of mass velocity is
   * required to calculate the elastic cross-section. The default is false.
   */ 
   void set_CMtrue_elastic(bool val);

  /*
   * Method: set_CMtrue_CX
   * Usage:  set_CMtrue_CX(true);
   * ----------------------------------
   * Set center_mass_CX to true if the ion-atom center of mass velocity is
   * required to calculate the charge-exchange cross-section. The default 
   * is false.
   */ 
   void set_CMtrue_CX(bool val);

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
   * Method: set_constant_CX
   * Usage:  set_constant_CX(constant_cross_section_val);
   * -----------------------------------------------------------
   * Sets the constant charge-exchange cross-section value, if it is constant.
   * This variable does not need to be set. If the cross-section is not constant
   * do not set this value, otherwise, this must be set -in addition- to setting
   * the cross-section function pointer to point to the
   * sigma_const_cx method defined in this file.
   */
   void set_constant_CX(double val);

  /*
   * Method: set_crossSec_elastic
   * Usage:  set_crossSec_elastic(&ion_collisions::elastic_crossSec_function);
   * ------------------------------------------------------------------------
   * Sets the 'double (ion_collisions::*sigma_elastic)(void)' function pointer
   * to point to the function passed in the arguement. User defined functions
   * can be passed. This function returns energy-dependent elastic cross 
   * sections.
   */
   void set_crossSec_elastic(double (ion_collisions::*func1)(void));

  /*
   * Method: set_crossSec_CX
   * Usage: set_crossSec_CX(&ion_collisions::charge_exchange_crossSec_func);
   * ------------------------------------------------------------------------
   * Sets the 'double (ion_collisions::*sigma_CX)(void)' function pointer
   * to point to the function passed in the arguement. User defined functions
   * can be passed. This function returns energy-dependent charge-exchange
   * cross sections.
   */
   void set_crossSec_CX(double (ion_collisions::*func1)(void));

  /*
   * Method: sigma_const_elastic
   * Usage: object.set_crossSec_elastic(&ion_collisions::sigma_const_elastic);
   * -------------------------------------------------------------------------
   * Passed as an arguement to the set_crossSec_elastic method to implement
   * the constant elastic cross-section defined using the
   * set_constant_elastic method. Note: the set_constant_elastic(value) must
   * be called first.
   */
   double sigma_const_elastic(void);

  /*
   * Method: sigma_const_cx
   * Usage: object.set_crossSec_CX(&ion_collisions::sigma_const_cx);
   * -------------------------------------------------------------------------
   * Passed as an arguement to the set_crossSec_CX method to implement
   * the constant charge-exchange cross-section defined using the
   * set_constant_CX method. Note: the set_constant_CX(value) method must
   * be called first.
   */
   double sigma_const_cx(void);

   /* A set of sample methods that return energy dependent cross-sections
    * to be passed as an arguement to the set_crossSec_elastic or
    * set_crossSec_CX methods. In general the user is expected to
    * supply functions.
    */
   double sigma_Ar_elastic(void);
   double sigma_Ar_cx(void);

   /*User defined ion cross-section method declarations*/
   #include "input_files/ion_coll_input_declare"

   private:
     double energy, CMenergy;
     bool center_mass_El,center_mass_CX;
     double constant_elastic;
     double constant_cx;

   /*
    * Method: (ion_collisions::*sigma_elastic)(void)
    * Usage:  double sigma_elastic = (this->*sigma_elastic)();
    * ---------------------------------------------------
    * Points to a function selected by the user. This function returns
    * the elastic cross-section.
    */
     double (ion_collisions::*sigma_elastic)(void);

   /*
    * Method: (ion_collisions::*sigma_CX)(void)
    * Usage:  double sigma_CX = (this->*sigma_CX)();
    * ---------------------------------------------------
    * Points to a function selected by the user. This function returns
    * the charge-exchange cross-section.
    */
     double (ion_collisions::*sigma_CX)(void);
     
};

#endif
