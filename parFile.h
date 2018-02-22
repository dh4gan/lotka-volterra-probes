/*
 * parFile.h
 *
 *  Created on: Sep 23, 2013
 *      Author: davidharvey
 */

#ifndef PARFILE_H_
#define PARFILE_H_

#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include <string.h>

#include <fstream>
#include <sstream>
using namespace std;

class parFile {
public:
	parFile( );
	parFile(string name);

	string parFileName;

	int nVertices;
	int iseed;
	double initialPrey;
	double initialPred;
	double tmax;
	double dt;

	string icChoice; // "GHZ" or "cluster"
	string parChoice; // "constant" = constant parameters, "uniform" = random sampling, "gaussian" =Gaussian sampling

	double rmin;
	double rmax;
	double scaleLength;

	double range;


	// Variable 1 is the constant or minimum value selected, or the mean for gaussian sampling
	// Variable 2 is the maximum value for uniform sampling, or the standard deviation for Gaussian sampling

	double predGrow1;
	double predGrow2;

	double predDeath1;
	double predDeath2;

	double preyGrow1;
	double preyGrow2;

	double preyDeath1;
	double preyDeath2;

	double mutationRate1;
	double mutationRate2;

	double outflowRate1;
	double outflowRate2;

	double velocity1;
	double velocity2;


	void readParams();
	void writeParamsToFile();
	void writeParamsToScreen();

};


#endif /* PARFILE_H_ */
