/*
 * Constants.h
 *
 *  Created on: Sep 16, 2012
 *      Author: dhf
 *      NOTE ALL CONSTANTS ARE IN SI UNITS
 */

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#include <stdlib.h>

const double pi = 3.141592658285;
const double Gsi = 6.67e-11;
const double Gmau = 1.0;  // Value of G for solar mass-AU units (time units...2pi units= 1 year)
const double Gmkpc = 4.49e-12; // G in kpc^3 Msol^{-1} kpc^{-2}
const double AU = 1.496e11;
const double msol = 1.99e30;
const double rsol = 6.95e8;
const double kpc = 3.08e19;
const double AUtokpc = AU/kpc;
const double rsoltokpc = rsol/kpc;




#endif /* CONSTANTS_H_ */
