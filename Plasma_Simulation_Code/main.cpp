/* 
 * File: main.cpp
 * ----------------
 * This file implements the main function, which initializes a simulation via a
 * pic(particle-in-cell) object of the pic class, and runs it using the
 * the pic_run() method.
 *
 */

#include "main.h" 
#include <time.h>
#include "input_files/constants.h"
#include "input_files/parameters.h"

int main(){

  /* Set the random seed to a constant to repeat a simulation. */
   srand(RAND_SEED);


  /* Need to initialize restart to true via the constructor in order to 
   * start a new simulation from dump files.         
   *
   * e.g. pic pic1(true);  
   */

   pic pic1(RESTART);
   pic1.pic::pic_run();

   return 0;
}
