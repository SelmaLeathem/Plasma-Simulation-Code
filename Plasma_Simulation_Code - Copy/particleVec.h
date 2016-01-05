/*
 * File: particleVec.h
 * -------------------
 * This file defines a class that holds each particle's position, velocities,
 * and wieght. 
 *
 * 4/15
 */

using namespace std;

#ifndef particleVec_h
#define particleVec_h

class particleVec
{
    public:
      double x;        /*position along x axis*/
      double v;        /*velocity along x axis*/
      double vy;       /*y axis velocity*/
      double vz;       /*z axis velocity*/
      double weight;   /*wieght of particle*/

      /*
       * Default constructor: particleVec 
       * Usage: particleVec;
       */ 
      particleVec(); 

      /*
       * Constructor: particleVec
       * Usage: particleVec particle
                particleVec particle(x,x_velocity,y_velocity,z_velocity,wieght);
       */
      particleVec(double x_in,double v_in,double vy_in,double vz_in,double weight_in);

      /*  
       * Destructor
       */
      ~particleVec();

};
#endif /*particleVec*/

 
