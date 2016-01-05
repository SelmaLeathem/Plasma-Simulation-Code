/*
 * File: fields.h
 * ----------------
 * This file exports the fields class which stores data and methods that store
 * and calculate field values on the grid such as the charge density, electric-
 * potential, and the electric field.
 *
 * 4/15
 *
 */

#ifndef fields_h
#define fields_h

#include "input_files/constants.h"
#include "input_files/parameters.h"
#include "particle.h"

using namespace std;

class fields
{
   public:
    
    /*
     * Constructor: fields
     * Usage: fields field_object;
     *
     */
     fields();

    /* Destructor: ~fields */
     ~fields();

    /*
     * Method: grid_rho
     * Usage:  grid_rho(ion_objects,e_grid_charge);
     * ----------------------------------------------
     * Calculates the net charge density on each grid point.
     *
     */
     void grid_rho(particle *ion,double *qe);

    /*
     * Method: grid_potential
     * Usage:  grid_potential(lhs_boundary_potential,rhs_boundary_potential);
     * ------------------------------------------------------------------------
     * Calculates the the electrostatic potential at each grid point due to 
     * the net charge density as defined by Poisson's equation.
     *
     */
     void grid_potential(double phi0, double phiL);

    /*
     * Method: grid_electric_field
     * Usage:  grid_electric_field();
     * -------------------------------     
     * Calculates the electric field on the grid from the electrostatic  
     * potential.
     */
     void grid_electric_field(void);

    /* Return pointers to private field variables. These are called from 
     * within the graph  and particle classes.
     *
     *     double* get_elec_field(void){return electric_field;}
     *     double* get_rho(void){return rho;}
     *     double* get_phi(void){return phi;}
     *     double  get_phi(int elem){return phi[elem];}
     */

    /*
     * Method: get_elec_field
     * Usage: double *E_field = object.fields::get_elec_field();
     * ---------------------------------------------------------
     * Returns a pointer to the electric field array. This is called
     * from the particle and graph class methods.
     */ 
     double* get_elec_field(void){return electric_field;}
 
    /*
     * Method: get_rho
     * Usage: double *rho= object.fields::get_rho();
     * ----------------------------------------------
     * Returns a pointer to the net charge density array. This is called
     * by methods in the graph class.
     */
     double* get_rho(void){return rho;}

    /*
     * Method: get_phi
     * Usage: double *phi= object.fields::get_phi();
     * -----------------------------------------------
     * Returns a pointer to the electrostatic potential array. This is 
     * called by methods in the graph class.
     */
     double* get_phi(void){return phi;}

    /*
     * Method: get_phi
     * Usage: double phi_element= object.fields::get_phi(element_number);
     * -------------------------------------------------------------------
     * Returns the electrostatic potential at the specified grid point. This
     * is used to graph the variation of the potential at a specific point in
     * time.
     */
     double  get_phi(int elem){return phi[elem];}

     private:

     /* Instance variables */

     double *rho;
     double *phi;
     double *electric_field;

     /* a,b,c,gam arrays are used to solve a tridiagonal matrix */
     double *a;
     double *b;
     double *c;
     double *gam;

};

#endif
