/*
 *
 *
 *  Created on: Sep 28, 2017
 *      Author: dh4gan
 *
 *  tests single LKSystem object
 *
 */

#include "LKVertex.h"
#include <iostream>
using namespace std;

int main()
    {

    // Set up single LK System

    double preyGrow = 0.6666;
    double preyDeath = 1.333;

    double predGrow = 1.0;
    double predDeath = 1.0;

    double mutate = 0.0;
    double outflow = 0.0;
    double velocity = 0.0;
    double t0 = 0.0;

    double initialPrey = 1.8;
    double initialPred = 1.8;

    double dt =0.0001;

    double tmax = 5.0;
    double t = 0.0;
    int ID = 1;

    LKVertex test(ID,initialPrey, initialPred, preyGrow,
	    preyDeath,predGrow, predDeath, mutate,
	    outflow, velocity, t0);

    // do test integrations

    test.initialiseSystem(t0, dt, initialPrey, initialPred);

    while(t < tmax)
	{

	test.updateLKSystem(t);
	cout << t << endl;
	test.writeToFile(t);
	t = t+dt;
	}


    }
