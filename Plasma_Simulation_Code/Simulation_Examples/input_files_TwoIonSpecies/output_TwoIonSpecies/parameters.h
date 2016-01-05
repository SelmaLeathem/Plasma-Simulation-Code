/*
 * File: parameters.h
 * ------------------  
 * This file temporarily houses all input parameters to the simulation.
 *
 *
 */

/*
 * This is a sample parameter file used to model a radio-frequency 
 * discharge over 600 RF cycles with two species: Argon and Hydrogen.  
 *
 * The simulation takes about 4-6 hours to run on a computer with 4GB RAM,
 * Intel core duo E8400 @3GHz processor, and a 64bit operating system.
 *
 */


#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <fstream>
#include <string.h>

#include "constants.h"

using namespace std;

#ifndef parameters_h
#define parameters_h


//constants used directly by program:
//Random seed generator. Supply an integer to repeat a simultion.
#define RAND_SEED 1 // time(NULL)

//At the end of each simulation particle data is dumped to textfiles.
//These may be used to start a new simulation where an older one left off
//by setting RESTART to true
#define RESTART false

//Create an array of neutral velocities to randomly pick from during a
//collision
#define INIT_NO_NEUTRALS (int)10000

//Number of ions and electrons to inject uniformly throughout the domain, per
//time-step. In the future will have option to choose point source instead.
#define NO_INJECT (int)0  //inject per timestep     

//Limit the total number of ionization events per timestep, otherwise
//sometimes have issues of runaway ionization.
#define IONIZE_LIMIT (int)15

//Electron particle density is used to calculate electron plasma frequency
//and the Debye length. This density is the order of magnitude that I expect
//the simulation to attain after ionization. 
#define ELECTRON_PARTICLE_DENSITY (double)1.0e15
#define ION_PARTICLE_DENSITY (double)1.0e15

//The macroparticle sizes are input parameters to the program. In this case
//they so happen to equall the density reduced to two-dimensions. In otherwords
//they represent charge sheets.
#define LOADING_DENSITY (double)1.0e12
#define ELECTRON_1D_DENSITY pow(LOADING_DENSITY,(1.0/3.0))
#define ION_1D_DENSITY pow(LOADING_DENSITY,(1.0/3.0))
#define ELECTRON_2D_DENSITY (double)1.0e5 // pow(LOADING_DENSITY,(2.0/3.0))
#define ION_2D_DENSITY (double)1.0e5  //pow(LOADING_DENSITY,(2.0/3.0))
#define E_MACROPARTICLE_SIZE ELECTRON_2D_DENSITY
#define ION_MACROPARTICLE_SIZE ION_2D_DENSITY

//Gas pressure in pascals.
#define GAS_PRESSURE ((double)0.5*((double)0.2)*((double)133.32))
#define NEUTRAL_DENSITY (GAS_PRESSURE/((double)0.025*CHARGE))
#define NEUTRAL_DENSITY2  NEUTRAL_DENSITY

//Thermal velocities, VTE,VTI,VTN are inputs to the program. 
#define TEMPERATURE  (double)1.5
#define VTE sqrt(TEMPERATURE*CHARGE/ELECTRON_MASS)
#define ION_MASS ((double)39.948*AMU) 
#define ION_CHARGE CHARGE
#define VTI sqrt(TEMPERATURE*CHARGE/ION_MASS)
#define VTN sqrt((double)0.025*CHARGE/ION_MASS)
#define ION_MASS2 AMU
#define VTI2 sqrt(TEMPERATURE*CHARGE/ION_MASS2)
#define VTN2 sqrt((double)0.025*CHARGE/ION_MASS2)

//Used for secondary electrons and ionization 
#define THRESHOLD (double)15.75  //Argon ionization threshold
#define WORK_FUNCTION (double)4.6 //Work function of surface
#define THRESHOLD2 (double)13.6  //Hydrogen ionization threshold

//LENGTH is the input parameter to the program giving the length of the system
//in meters.
#define NO_OF_DEBYE_LENGTHS (int)105
#define CELLS_PER_DEBYE_LENGTH (int)4
#define DEBYE_LENGTH sqrt((EPSILON*TEMPERATURE)/(ELECTRON_PARTICLE_DENSITY*CHARGE))
#define LENGTH ((double)NO_OF_DEBYE_LENGTHS*DEBYE_LENGTH)

//NO_OF_CELLS is an input parameter to the program, giving the number of cells
//in length LENGTH
#define NO_OF_CELLS (NO_OF_DEBYE_LENGTHS*CELLS_PER_DEBYE_LENGTH)

//DT, DTI,DT_DTI,TOTAL_RUN_TIME,START_GRAPH_TIME,END_GRAPH_TIME, and SIZE_TIME_DATA 
//are input parameters to the program

#define PLASMA_FREQUENCY sqrt((ELECTRON_PARTICLE_DENSITY*ELECTRON_CHARGE2)/(EPSILON*ELECTRON_MASS))

//Electron time-step
#define DT ((double)0.1/PLASMA_FREQUENCY)
//Ion time-step
//A different DTI should not be used when implementing a steady_state_boundary.
#define DTI ((double)1.0*DT)
//Ratio of ion time-step over electron
#define DT_DTI (int)(DTI/DT)


#define POTENTIAL_FREQUENCY (double)13.6e6
#define TWO_PI_F_DT (2.0*DT*PI*POTENTIAL_FREQUENCY) 
#define PERIOD (1.0/TWO_PI_F_DT)

//number of steps to take average of data over
#define NO_GRAPH_STEPS ((int)100*(int)PERIOD)
//total run-time of simulation
#define TOTAL_RUN_TIME ((int)600*(int)PERIOD + 5)
//when to start collecting graph data points
#define START_GRAPH_TIME ((int)500*(int)PERIOD)

//when to stop collecting graph data points
#define END_GRAPH_TIME (START_GRAPH_TIME+NO_GRAPH_STEPS)

//SIZE_TIME_DATA is input parameter that defines the length of the arrays  
//storing the variety in data at a particular location in time
#define TIME_DATA_FREQUENCY (int)100
#define SIZE_TIME_DATA TOTAL_RUN_TIME/TIME_DATA_FREQUENCY +(int)1

//initial particle number loaded into the domain
#define INIT_PARTICLE_NOE (int)(ELECTRON_1D_DENSITY*LENGTH) 
#define INIT_PARTICLE_NOI (int)(0.5*ELECTRON_1D_DENSITY*LENGTH) 
#define INIT_PARTICLE_NOI2 INIT_PARTICLE_NOI  

//DX,QCE,QCI,Q_OVER_ME,Q_OVER_MI,DXPE DXPI are used directly in the program

//DX is the cell length
#define DX (LENGTH/(double)(NO_OF_CELLS-1))

//QC is the charge of each macroparticle
#define QCI (ION_MACROPARTICLE_SIZE*ION_CHARGE)
#define QCE (E_MACROPARTICLE_SIZE*ELECTRON_CHARGE)
#define QCI2 QCI

//q_over_m is independent of the size of the macroparticle and is used to
//in the equations of motion
#define Q_OVER_ME (ELECTRON_CHARGE/ELECTRON_MASS)
#define Q_OVER_MI (ION_CHARGE/ION_MASS)
#define Q_OVER_MI2 (ION_CHARGE/ION_MASS2)

//Dxp is the distance between particles during the initial load at the start
//of the simulation
#define DXPE  (LENGTH/(INIT_PARTICLE_NOE +(double)1.0))
#define DXPI  (LENGTH/(INIT_PARTICLE_NOI +(double)1.0))
#define DXPI2  (LENGTH/(INIT_PARTICLE_NOI2 +(double)1.0))

//LENGTH_MINUS_DX is used directly by the program
#define LENGTH_MINUS_DX (LENGTH - DX)

//optional magnetic field on x-y plane, Bx = Bcos(angle),By=Bsin(angle)
//set magnetic field =0 for no field
#define  MAGNETIC_FIELD (double)0.0
#define  ANGLE ((double)30.0)*PI/((double)180.0)

//boundary voltage
#define PHI_LEFT (double)0.0
//#define PHI_RIGHT ((double)100.0*cos(TWO_PI_F_DT*(double)time_step))
#define PHI_RIGHT (double)100.0*sin(TWO_PI_F_DT*(double)time_step)

//NO_ION_SPECIES is an input parameter describing the number of species
#define NO_ION_SPECIES (int)2

//Parameters used to make text files of particle velocity distributions
//Number of velocity bins currently corresponds to the max number of rows in
//excel workbooks
#define NO_BINS 30001
//Bin scaling. This will become an internal variable 
#define DIVIDE_BINS_BY ((NO_BINS-1)/((int)100))
//Locations of where the particle velocity distributions are obtained
#define DIST_LOCATION  (int)((float)0.9*(float)NO_OF_CELLS)
#define DIST_LOCATION2 (int)((float)0.1*(float)NO_OF_CELLS)
//The number of cells around the dist_location from which velocity data is 
//pulled
#define SPREAD (int)((float)0.05*(float)NO_OF_CELLS)

//Which cell numbers to take time-variable data from
#define POINT1  (int)3
#define POINT2  (int)((double)NO_OF_CELLS/2.0)
#define POINT3  (int)(NO_OF_CELLS -3)

//Names used to distinguish text files of plot data
#define NO_TYPE (char *)"n"
#define ION_TYPE (char *)"i1"
#define ION_TYPE2 (char *)"i2"
#define ELECTRON_TYPE (char *)"e"


//parameters used for reducing the particle number when the particle vectors
//get too large, one group example:

//Number of velocity groups to account for when reducing the particle number.  
//Choose more than one velocity group to preserve high velocity tails in the
//random reduction of the number of particles. This is a mandatory constant 
//that is also used to construct group arrays.
//#define NO_GROUPS (int)1  //(int)2  

//Maximum allowable particle number of each species in each velocity group.
//Note: Different groups can have different sizes.
//#define UPPER_GROUP_LIMIT (int)1000000
//Desired group size after reduction
//#define NEW_GROUP_SIZE (int)500000
//The largest velocity value to include in the reduction. Note: can exclude tail//tails.
//#define UPPER_V_LIMIT ((double)200.0*VTE)
//#define UPPER_V_LIMIT_I ((double)200.0*VTI)


//more than one group
#define NO_GROUPS (int)3   
#define UPPER_GROUP_LIMIT (int)400000
#define UPPER_GROUP_LIMIT2 (int)100000
#define NEW_GROUP_SIZE (int)200000
#define NEW_GROUP_SIZE2 (int)50000
#define UPPER_V_LIMIT ((double)1000.0*VTE)
#define UPPER_V_LIMIT_I ((double)1000.0*VTI)
#define UPPER_V_LIMIT_I2 ((double)1000.0*VTI2)

#define VI1 (double)2.0*VTI
#define VI2 (double)3.0*VTI
#define VI12 (double)2.0*VTI2
#define VI22 (double)3.0*VTI2
#define VE1 (double)2.0*VTE
#define VE2 (double)3.0*VTE


//Putting particles into groups to check each group size is time consuming
//so this is only done when the number of particles reaches the below limit.
//Set to one to check each group size every time-step.
#define SIZE_TO_CHECK (int)600001
#define SIZE_TO_CHECK_E (int)1000001

//For a steady boundary define how a particle who reaches the wall is 
//placed back into the simulation.
//random_number = a random number between 0 and 1
//LENGTH = the domain length
//DX = cell length
#define PLACE_AT (DX+random_number*(LENGTH-2.0*DX)) //0.5*LENGTH

//Number of ions and electrons to inject uniformly throughout the domain, per
//time-step. In the future will have option to choose point source instead.
#define NO_INJECT (int)0  //inject per timestep
//Define where to place the injected particle.
//random_number= a random number between 0 and 1
//Be aware that placing particles at a single point can result in large
//localized electric fields that might accelerate electrons to relativistic
//values.

#define PUT_PARTICLE_AT (DX+random_number*(LENGTH-2.0*DX)) //0.5*LENGTH

#endif
