/*
 * File: particleVec.cpp
 * ---------------------
 * This file implements the particleVec.h interface
 *
 */


#include "particleVec.h"


particleVec::particleVec(void){
   x=v=vy=vz= 0.0;
   weight=1.0;
}

particleVec::particleVec(double x_in,double v_in,double vy_in,double vz_in,double weight_in):x(x_in),v(v_in),vy(vy_in),vz(vz_in),weight(weight_in)
{
    
}

particleVec::~particleVec(void)
{

}

